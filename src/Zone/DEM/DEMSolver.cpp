#include "DEMSolver.h"
#include "LinearSpringDashpot.h"
#include "HertzMindlin.h"
#include "KDTree.h"
#include "Log.h"
#include "ZaranError.h"
#include "CommonPara.h"
#include <fstream>
#include <algorithm>
#include <cmath>
#include <set>
#include <vtkCellArray.h>
#include <vtkDelaunay2D.h>
#include <vtkIdList.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>

namespace
{
using Polygon2D = std::vector<Eigen::Vector2d>;

Polygon2D ClipVoronoiHalfPlane(const Polygon2D& polygon,
                               const Eigen::Vector2d& normal,
                               double offset)
{
    Polygon2D result;
    if (polygon.empty()) return result;
    const double tolerance = 1.0e-13 * std::max(1.0, std::abs(offset));
    for (std::size_t i = 0; i < polygon.size(); ++i)
    {
        const Eigen::Vector2d& current = polygon[i];
        const Eigen::Vector2d& previous = polygon[(i + polygon.size() - 1) % polygon.size()];
        const double current_value = normal.dot(current) - offset;
        const double previous_value = normal.dot(previous) - offset;
        const bool current_inside = current_value <= tolerance;
        const bool previous_inside = previous_value <= tolerance;
        if (current_inside != previous_inside)
        {
            const double denominator = previous_value - current_value;
            if (std::abs(denominator) > 1.0e-30)
                result.push_back(previous + (previous_value / denominator) * (current - previous));
        }
        if (current_inside) result.push_back(current);
    }
    return result;
}

double PolygonArea(const Polygon2D& polygon)
{
    double twice_area = 0.0;
    for (std::size_t i = 0; i < polygon.size(); ++i)
    {
        const auto& a = polygon[i];
        const auto& b = polygon[(i + 1) % polygon.size()];
        twice_area += a.x() * b.y() - a.y() * b.x();
    }
    return 0.5 * std::abs(twice_area);
}
}

namespace zaran
{

    DEMSolver::DEMSolver(index_type index,
        const std::string& name,
        shared_ptr<DEMSolverParam> para,
        shared_ptr<DEMFieldData>   dem_data)
        : FieldSolver(index, name, para, nullptr, nullptr), // DEM 无背景网格，无传统 DataManager
        m_dem_data(std::move(dem_data))
    {
    }

    DEMSolverParam* DEMSolver::GetDEMParam() const
    {
        return static_cast<DEMSolverParam*>(const_cast<DEMSolver*>(this)->GetPara());
    }

    void DEMSolver::Init()
    {
        InitSolver();
        InitField();
    }

    void DEMSolver::InitSolver()
    {
        auto* para = GetDEMParam();
        para->Init();

        const std::string& model_name = para->GetContactModel();
        if (model_name == "LinearSpringDashpot" || model_name == "Linear")
        {
            m_contact_model = make_unique<LinearSpringDashpot>();
            Log::info("DEMSolver: contact model = LinearSpringDashpot");
        }
        else if (model_name == "HertzMindlin" || model_name == "Hertz")
        {
            m_contact_model = make_unique<HertzMindlin>();
            Log::info("DEMSolver: contact model = HertzMindlin");
        }
        else
        {
            throw ZaranError("Unknown DEM contact model: " + model_name);
        }
    }

    void DEMSolver::InitField()
    {
        auto* para = GetDEMParam();
        // 将默认材料参数应用到每个粒子（若粒子未显式指定则用全局值）
        for (auto& p : m_dem_data->GetParticles())
        {
            if (p.radius <= 0.0 || p.inertia <= 0.0)
                throw ZaranError("DEM particle radius and inertia factor must be positive");
            if (p.mass <= 0.0)
            {
                // 默认球体：m = 4/3 π r³ ρ
                double rho = para->GetDensity();
                p.mass = (4.0 / 3.0) * PI * std::pow(p.radius, 3) * rho;
                p.inertia = 0.4; // I = 2/5 m r²
            }
            if (!p.pos.allFinite() || !p.vel.allFinite() || !p.omega.allFinite())
                throw ZaranError("DEM particle state contains a non-finite value");
            if (p.temperature <= 0.0 || p.reaction_progress < 0.0
                || p.reaction_progress > 1.0 || p.specific_heat <= 0.0
                || p.thermal_conductivity < 0.0 || p.reaction_heat < 0.0
                || p.arrhenius_prefactor < 0.0 || p.activation_temperature < 0.0)
                throw ZaranError("DEM particle thermal/reaction data are invalid");
            if (p.reaction_progress >= para->GetBurnoutThreshold()) p.phase = 2;
            else if (p.reaction_progress >= para->GetMixedPhaseThreshold()) p.phase = 1;
            else p.phase = 0;
            // 仅在粒子未从文件显式设置材料参数时，才使用全局参数覆盖
            if (!p.material_from_file)
            {
                p.young_modulus = para->GetYoungModulus();
                p.poisson_ratio = para->GetPoissonRatio();
                p.friction_coeff = para->GetFrictionCoeff();
                p.restitution_coeff = para->GetRestitutionCoeff();
            }
            if (p.reference_volume <= 0.0)
                p.reference_volume = p.mass / para->GetDensity();
            if (p.energetic && para->GetIgnitionRadius() > 0.0
                && (p.pos - para->GetIgnitionCenter()).norm() <= para->GetIgnitionRadius())
                p.temperature = std::max(p.temperature, para->GetIgnitionCenterTemperature());
            p.gas_temperature = p.temperature;
            if (p.reaction_progress > 0.0 && p.gas_internal_energy <= 0.0)
                p.gas_internal_energy = p.mass * p.reaction_progress
                * (p.specific_heat * p.temperature + p.reaction_heat);
            UpdateGasSolidState(p);
        }
        double rest_length_sum = 0.0;
        std::size_t rest_length_count = 0;
        for (auto& bond : m_dem_data->GetBonds())
        {
            if (bond.idx_a >= m_dem_data->GetParticleNum()
                || bond.idx_b >= m_dem_data->GetParticleNum()
                || bond.idx_a == bond.idx_b || bond.rest_length <= 0.0
                || bond.normal_stiffness < 0.0 || bond.tangential_stiffness < 0.0
                || bond.normal_damping < 0.0 || bond.tangential_damping < 0.0
                || bond.conduction_area < 0.0)
                throw ZaranError("DEM bond data are invalid");
            if (para->GetSurfaceEnergy() > 0.0)
            {
                const double boundary_ratio = (para->GetSurfaceEnergy()
                    - para->GetGrainBoundaryEnergy()) / para->GetSurfaceEnergy();
                bond.fracture_energy = para->GetSurfaceEnergy() * bond.rest_length
                    * para->GetLatticeThickness() / std::sqrt(3.0) * boundary_ratio;
            }
            const double interface_area = bond.conduction_area > 0.0
                ? bond.conduction_area : bond.rest_length * para->GetLatticeThickness();
            if (para->GetNormalViscosity() > 0.0)
                bond.normal_damping = para->GetNormalViscosity()
                    * interface_area / bond.rest_length;
            if (para->GetTangentialViscosity() > 0.0)
                bond.tangential_damping = para->GetTangentialViscosity()
                    * interface_area / bond.rest_length;
            rest_length_sum += bond.rest_length;
            ++rest_length_count;
        }
        if (rest_length_count > 0)
            m_lattice_spacing = rest_length_sum / static_cast<double>(rest_length_count);
        UpdateGasVoronoiMesh(true);
        Log::info("DEMSolver InitField: {} particles, {} bonds",
            m_dem_data->GetParticleNum(), m_dem_data->GetBondNum());
    }

    void DEMSolver::Preprocess()
    {
        // 当前版本无前处理工作（扩展时可加：MPI 通信等）
    }

    void DEMSolver::Postprocess()
    {
        // 当前版本无后处理工作
    }

    void DEMSolver::Solve()
    {
        m_contacts_this_step.clear();
        ZeroForce();
        CalcBondForce();
        ContactDetection();
        CalcContactForce();
        CalcWallForce();
        PruneContactHistory();
        CalcGravity();
        CalcThermalReaction();
        CalcGasPressureForce();
        Integrate();
        ++m_dem_step;
    }

    void DEMSolver::UpdateGasVoronoiMesh(bool force)
    {
        if (!GetDEMParam()->GetReactionEnabled() || !GetDEMParam()->GetVoronoiEnabled())
            return;
        if (!force && m_voronoi_valid
            && (m_dem_step == 0
                || (m_dem_step % static_cast<std::uint64_t>(GetDEMParam()->GetVoronoiUpdateInterval())) != 0))
            return;

        const auto& particles = m_dem_data->GetParticles();
        const index_type count = particles.size();
        if (count < 3) return;

        vtkNew<vtkPoints> points;
        points->SetNumberOfPoints(static_cast<vtkIdType>(count));
        for (index_type i = 0; i < count; ++i)
            points->SetPoint(static_cast<vtkIdType>(i), particles[i].pos.x(), particles[i].pos.y(), 0.0);
        vtkNew<vtkPolyData> input;
        input->SetPoints(points);
        vtkNew<vtkDelaunay2D> triangulation;
        triangulation->SetInputData(input);
        triangulation->SetTolerance(1.0e-10);
        triangulation->Update();

        std::vector<std::set<index_type>> neighbors(count);
        vtkCellArray* triangles = triangulation->GetOutput()->GetPolys();
        vtkNew<vtkIdList> ids;
        triangles->InitTraversal();
        while (triangles->GetNextCell(ids))
        {
            if (ids->GetNumberOfIds() != 3) continue;
            for (vtkIdType edge = 0; edge < 3; ++edge)
            {
                const index_type a = static_cast<index_type>(ids->GetId(edge));
                const index_type b = static_cast<index_type>(ids->GetId((edge + 1) % 3));
                if (a >= count || b >= count || a == b) continue;
                neighbors[a].insert(b);
                neighbors[b].insert(a);
            }
        }

        m_gas_voronoi_cells.assign(count, GasVoronoiCell{});
        const double xmin = GetDEMParam()->GetBoxXMin();
        const double xmax = GetDEMParam()->GetBoxXMax();
        const double ymin = GetDEMParam()->GetBoxYMin();
        const double ymax = GetDEMParam()->GetBoxYMax();
        const Polygon2D box = {{xmin, ymin}, {xmax, ymin}, {xmax, ymax}, {xmin, ymax}};
        for (index_type i = 0; i < count; ++i)
        {
            Polygon2D polygon = box;
            const Eigen::Vector2d pi(particles[i].pos.x(), particles[i].pos.y());
            for (index_type j : neighbors[i])
            {
                const Eigen::Vector2d pj(particles[j].pos.x(), particles[j].pos.y());
                const Eigen::Vector2d normal = pj - pi;
                if (normal.squaredNorm() <= 1.0e-30) continue;
                const double offset = 0.5 * (pj.squaredNorm() - pi.squaredNorm());
                polygon = ClipVoronoiHalfPlane(polygon, normal, offset);
                if (polygon.empty()) break;
            }
            m_gas_voronoi_cells[i].vertices = std::move(polygon);
            m_gas_voronoi_cells[i].area = PolygonArea(m_gas_voronoi_cells[i].vertices);
        }

        m_gas_voronoi_faces.clear();
        const double domain_diagonal = std::hypot(xmax - xmin, ymax - ymin);
        const double face_tolerance = std::max(1.0e-12, 1.0e-8 * domain_diagonal);
        for (index_type i = 0; i < count; ++i)
        {
            const Eigen::Vector2d pi(particles[i].pos.x(), particles[i].pos.y());
            for (index_type j : neighbors[i])
            {
                if (j <= i) continue;
                const Eigen::Vector2d pj(particles[j].pos.x(), particles[j].pos.y());
                const Eigen::Vector2d normal = pj - pi;
                const double normal_length = normal.norm();
                if (normal_length <= 1.0e-15) continue;
                const double offset = 0.5 * (pj.squaredNorm() - pi.squaredNorm());
                double face_length = 0.0;
                const auto& polygon = m_gas_voronoi_cells[i].vertices;
                for (std::size_t edge = 0; edge < polygon.size(); ++edge)
                {
                    const auto& a = polygon[edge];
                    const auto& b = polygon[(edge + 1) % polygon.size()];
                    const double da = std::abs(normal.dot(a) - offset) / normal_length;
                    const double db = std::abs(normal.dot(b) - offset) / normal_length;
                    if (da <= face_tolerance && db <= face_tolerance)
                        face_length = std::max(face_length, (b - a).norm());
                }
                if (face_length > face_tolerance)
                    m_gas_voronoi_faces.push_back({i, j, face_length,
                        face_length * GetDEMParam()->GetLatticeThickness()});
            }
        }
        m_voronoi_valid = !m_gas_voronoi_faces.empty();
        if (force)
            Log::info("DEM gas Voronoi: {} cells, {} internal faces",
                      m_gas_voronoi_cells.size(), m_gas_voronoi_faces.size());
    }

    double DEMSolver::JwlPressure(const DEMParticle& p, double gas_volume) const
    {
        if (p.reaction_progress <= 0.0 || p.reference_volume <= 0.0 || gas_volume <= 0.0)
            return 0.0;
        const double reacted_reference_volume = p.reference_volume * p.reaction_progress;
        const double V = std::max(0.25, gas_volume / reacted_reference_volume);
        const double energy_density = p.gas_internal_energy / reacted_reference_volume;
        const double pressure = GetDEMParam()->GetJwlA()
            * (1.0 - GetDEMParam()->GetJwlOmega() / (GetDEMParam()->GetJwlR1() * V))
            * std::exp(-GetDEMParam()->GetJwlR1() * V)
            + GetDEMParam()->GetJwlB()
            * (1.0 - GetDEMParam()->GetJwlOmega() / (GetDEMParam()->GetJwlR2() * V))
            * std::exp(-GetDEMParam()->GetJwlR2() * V)
            + GetDEMParam()->GetJwlOmega() * energy_density / V;
        return std::min(GetDEMParam()->GetPressureCap(), std::max(0.0, pressure));
    }

    void DEMSolver::UpdateGasSolidState(DEMParticle& p) const
    {
        const double alpha = std::min(1.0, std::max(0.0, p.reaction_progress));
        p.total_volume = p.reference_volume * p.volume_ratio;
        const double thickness = GetDEMParam()->GetLatticeThickness();
        p.gas_radius = std::sqrt(std::max(0.0, p.total_volume) / (PI * thickness));
        if (alpha <= 0.0)
        {
            p.solid_volume = p.total_volume;
            p.gas_volume = 0.0;
            p.solid_core_radius = p.gas_radius;
            p.gas_temperature = p.temperature;
            p.gas_pressure = 0.0;
            return;
        }

        const double gas_mass = p.mass * alpha;
        if (p.gas_internal_energy <= 0.0)
            p.gas_internal_energy = gas_mass * p.specific_heat * p.temperature;
        p.gas_temperature = std::max(1.0,
            p.gas_internal_energy / (gas_mass * p.specific_heat));

        if (alpha >= GetDEMParam()->GetBurnoutThreshold())
        {
            p.solid_volume = 0.0;
            p.gas_volume = p.total_volume;
        }
        else
        {
            const double Vg0 = p.reference_volume * alpha;
            const double Vs0 = p.reference_volume * (1.0 - alpha);
            const double bulk_modulus = p.young_modulus
                / (3.0 * (1.0 - 2.0 * p.poisson_ratio));
            double lo = std::max(0.25 * Vg0, p.total_volume - Vs0);
            double hi = std::max(lo, p.total_volume * (1.0 - 1.0e-12));
            auto residual = [&](double gas_volume)
            {
                    const double solid_volume = p.total_volume - gas_volume;
                    const double strain = std::min(0.75,
                        std::max(0.0, 1.0 - solid_volume / Vs0));
                    double solid_pressure = bulk_modulus * strain;
                    if (GetDEMParam()->GetHighPressureKnn() > 0.0
                        && GetDEMParam()->GetHighPressureRnn() > 0.0)
                    {
                        // LSM gas_phase_mod 的二维非线性体积压缩闭合。
                        const double exponent = GetDEMParam()->GetHighPressureExponent();
                        const double radial_compression = 1.0 - std::sqrt(1.0 - strain);
                        const double factor = bulk_modulus * GetDEMParam()->GetHighPressureKnn()
                            / std::pow(1.0 - GetDEMParam()->GetHighPressureRnn(), exponent);
                        solid_pressure = 2.0 * factor
                            * std::pow(radial_compression, exponent + 1.0)
                            * (1.0 / (exponent + 1.0)
                               - radial_compression / (exponent + 2.0));
                    }
                    return JwlPressure(p, gas_volume) - solid_pressure;
            };
            double flo = residual(lo);
            double fhi = residual(hi);
            if (flo * fhi <= 0.0)
            {
                for (int iteration = 0; iteration < 60; ++iteration)
                {
                    const double mid = 0.5 * (lo + hi);
                    const double fmid = residual(mid);
                    if (flo * fmid <= 0.0) { hi = mid; fhi = fmid; }
                    else { lo = mid; flo = fmid; }
                }
                p.gas_volume = 0.5 * (lo + hi);
            }
            else
            {
                p.gas_volume = std::abs(flo) < std::abs(fhi) ? lo : hi;
            }
            p.solid_volume = std::max(0.0, p.total_volume - p.gas_volume);
        }
        p.solid_core_radius = std::sqrt(p.solid_volume / (PI * thickness));
        p.gas_pressure = JwlPressure(p, p.gas_volume);
    }

    void DEMSolver::CalcThermalReaction()
    {
        if (!GetDEMParam()->GetReactionEnabled()) return;

        UpdateGasVoronoiMesh();
        auto& particles = m_dem_data->GetParticles();
        const double dt = GetDEMParam()->GetTimeStep();
        std::vector<double> energy_delta(particles.size(), 0.0);
        std::vector<double> surface_delta(particles.size(), 0.0);
        std::vector<double> core_delta(particles.size(), 0.0);
        std::vector<double> length_ratio_sum(particles.size(), 0.0);
        std::vector<int> length_ratio_count(particles.size(), 0);
        std::vector<double> gas_pressure_sum(particles.size(), 0.0);
        std::vector<int> gas_pressure_count(particles.size(), 1);

        if (!m_voronoi_valid)
        {
            for (const auto& bond : m_dem_data->GetBonds())
            {
                if (bond.rest_length <= 0.0) continue;
                const double ratio = (particles[bond.idx_b].pos - particles[bond.idx_a].pos).norm()
                    / bond.rest_length;
                length_ratio_sum[bond.idx_a] += ratio * ratio;
                length_ratio_sum[bond.idx_b] += ratio * ratio;
                ++length_ratio_count[bond.idx_a];
                ++length_ratio_count[bond.idx_b];
            }
        }
        for (index_type i = 0; i < particles.size(); ++i)
        {
            DEMParticle& p = particles[i];
            if (m_voronoi_valid && i < m_gas_voronoi_cells.size()
                && m_gas_voronoi_cells[i].area > 0.0 && p.reference_volume > 0.0)
            {
                const double voronoi_volume = m_gas_voronoi_cells[i].area
                                            * GetDEMParam()->GetLatticeThickness();
                p.volume_ratio = std::min(GetDEMParam()->GetVoronoiMaxVolumeRatio(),
                    std::max(0.25, voronoi_volume / p.reference_volume));
            }
            else if (length_ratio_count[i] > 0)
                p.volume_ratio = std::max(0.25, length_ratio_sum[i] / length_ratio_count[i]);

            UpdateGasSolidState(p);
            gas_pressure_sum[i] = p.gas_pressure;
        }

        if (m_voronoi_valid)
        {
            for (const auto& face : m_gas_voronoi_faces)
            {
                const DEMParticle& pa = particles[face.idx_a];
                const DEMParticle& pb = particles[face.idx_b];
                if (pb.reaction_progress > 0.0)
                {
                    gas_pressure_sum[face.idx_a] += pb.gas_pressure;
                    ++gas_pressure_count[face.idx_a];
                }
                if (pa.reaction_progress > 0.0)
                {
                    gas_pressure_sum[face.idx_b] += pa.gas_pressure;
                    ++gas_pressure_count[face.idx_b];
                }
            }
        }
        else
        {
            for (const auto& bond : m_dem_data->GetBonds())
            {
                const DEMParticle& pa = particles[bond.idx_a];
                const DEMParticle& pb = particles[bond.idx_b];
                if (pb.reaction_progress > 0.0)
                {
                    gas_pressure_sum[bond.idx_a] += pb.gas_pressure;
                    ++gas_pressure_count[bond.idx_a];
                }
                if (pa.reaction_progress > 0.0)
                {
                    gas_pressure_sum[bond.idx_b] += pa.gas_pressure;
                    ++gas_pressure_count[bond.idx_b];
                }
            }
        }

        // 气固共存点内部传热和固相核表面燃烧。LSM 使用二维固相核周长乘厚度
        // 作为界面面积，并用本点与相邻气相点的平均压力计算燃速。
        for (index_type i = 0; i < particles.size(); ++i)
        {
            DEMParticle& p = particles[i];
            p.internal_heat_transfer = 0.0;
            if (p.phase != 1 || p.solid_core_radius <= 0.0 || p.gas_volume <= 0.0)
                continue;
            const double interface_area = 2.0 * PI * p.solid_core_radius
                * GetDEMParam()->GetLatticeThickness();
            const double transfer_distance = 0.5 * (p.solid_core_radius + p.gas_radius);
            if (transfer_distance > 0.0 && p.thermal_conductivity > 0.0)
            {
                double transfer = p.thermal_conductivity
                    * (p.gas_temperature - p.temperature) * interface_area
                    / transfer_distance * dt;
                const double gas_limit = GetDEMParam()->GetInternalHeatLimit()
                    * p.gas_internal_energy;
                const double solid_energy = p.mass * (1.0 - p.reaction_progress)
                    * p.specific_heat * p.temperature;
                const double solid_limit = GetDEMParam()->GetInternalHeatLimit() * solid_energy;
                if (transfer >= 0.0) transfer = std::min(transfer, gas_limit);
                else transfer = std::max(transfer, -solid_limit);
                p.internal_heat_transfer = transfer;
                p.gas_internal_energy -= transfer;
                energy_delta[i] += transfer;
            }

            const double average_pressure = gas_pressure_sum[i] / gas_pressure_count[i];
            if (average_pressure > GetDEMParam()->GetBurnPressureThreshold()
                && p.gas_temperature > GetDEMParam()->GetBurnTemperatureThreshold())
            {
                double burn_rate = 1.0e-3 * GetDEMParam()->GetBurnA()
                    * std::pow(average_pressure * 1.0e-6, GetDEMParam()->GetBurnB());
                if (p.gas_temperature < GetDEMParam()->GetGasTemperatureReference())
                    burn_rate *= p.gas_temperature / GetDEMParam()->GetGasTemperatureReference();
                core_delta[i] = burn_rate * interface_area * dt / p.reference_volume;
            }
        }

        // LSM 的固相连接传热：q = k * A * (Tj - Ti) / L。
        for (auto& bond : m_dem_data->GetBonds())
        {
            bond.heat_flow_a = 0.0;
            if (!bond.active || bond.conduction_area <= 0.0) continue;
            const DEMParticle& pa = particles[bond.idx_a];
            const DEMParticle& pb = particles[bond.idx_b];
            const double length = (pb.pos - pa.pos).norm();
            if (length <= 1.0e-15) continue;
            double conductivity = 0.0;
            if (pa.thermal_conductivity > 0.0 && pb.thermal_conductivity > 0.0)
                conductivity = 2.0 * pa.thermal_conductivity * pb.thermal_conductivity
                / (pa.thermal_conductivity + pb.thermal_conductivity);
            if (conductivity <= 0.0) continue;
            bond.heat_flow_a = conductivity * bond.conduction_area
                * (pb.temperature - pa.temperature) / length;
            energy_delta[bond.idx_a] += bond.heat_flow_a * dt;
            energy_delta[bond.idx_b] -= bond.heat_flow_a * dt;
        }

        // LSM 的颗粒间表面燃烧：Voronoi 气固共边是实际燃烧面积。
        auto accumulate_surface_burn = [&](index_type source_index, index_type target_index,
                                           double interface_area)
        {
            const DEMParticle& source = particles[source_index];
            const DEMParticle& target = particles[target_index];
            if (!target.energetic || target.reaction_progress >= 1.0
                || source.reaction_progress <= 0.19
                || source.gas_pressure <= GetDEMParam()->GetBurnPressureThreshold()
                || source.gas_temperature <= GetDEMParam()->GetBurnTemperatureThreshold()
                || target.reference_volume <= 0.0 || interface_area <= 0.0)
                return;
            double burn_rate = 1.0e-3 * GetDEMParam()->GetBurnA()
                * std::pow(source.gas_pressure * 1.0e-6, GetDEMParam()->GetBurnB());
            if (source.gas_temperature < GetDEMParam()->GetGasTemperatureReference())
                burn_rate *= source.gas_temperature / GetDEMParam()->GetGasTemperatureReference();
            surface_delta[target_index] += burn_rate * interface_area * dt
                                         / target.reference_volume;
        };
        if (m_voronoi_valid)
        {
            for (const auto& face : m_gas_voronoi_faces)
            {
                accumulate_surface_burn(face.idx_a, face.idx_b, face.area);
                accumulate_surface_burn(face.idx_b, face.idx_a, face.area);
            }
        }
        else
        {
            for (const auto& bond : m_dem_data->GetBonds())
            {
                accumulate_surface_burn(bond.idx_a, bond.idx_b, bond.conduction_area);
                accumulate_surface_burn(bond.idx_b, bond.idx_a, bond.conduction_area);
            }
        }

        for (index_type i = 0; i < particles.size(); ++i)
        {
            DEMParticle& p = particles[i];
            p.reaction_rate = 0.0;
            p.body_reaction_increment = 0.0;
            p.core_burn_increment = 0.0;
            p.neighbor_burn_increment = 0.0;
            energy_delta[i] += p.heat_source * dt;

            // 与 LSM calculation.f90 一致的一级 Arrhenius 基体热分解项：
            // d_alpha = Z (1-alpha) exp(-Ta/T) dt，并在本步截断到 [0,1]。
            if (p.energetic
                && (p.temperature >= GetDEMParam()->GetIgnitionTemperature()
                    || p.reaction_progress >= GetDEMParam()->GetMixedPhaseThreshold())
                && p.reaction_progress < 1.0 && p.arrhenius_prefactor > 0.0
                && p.activation_temperature >= 0.0)
            {
                p.reaction_rate = p.arrhenius_prefactor * (1.0 - p.reaction_progress)
                    * std::exp(-p.activation_temperature / p.temperature);
                p.body_reaction_increment = std::max(0.0, p.reaction_rate * dt);
            }
            if (p.energetic && p.reaction_progress < 1.0)
            {
                p.core_burn_increment = std::max(0.0, core_delta[i]);
                p.neighbor_burn_increment = std::max(0.0, surface_delta[i]);
                double d_alpha = p.body_reaction_increment + p.core_burn_increment
                    + p.neighbor_burn_increment;
                const double remaining = 1.0 - p.reaction_progress;
                if (d_alpha > remaining && d_alpha > 0.0)
                {
                    const double scale = remaining / d_alpha;
                    p.body_reaction_increment *= scale;
                    p.core_burn_increment *= scale;
                    p.neighbor_burn_increment *= scale;
                    d_alpha = remaining;
                }
                const double old_alpha = p.reaction_progress;
                p.reaction_progress += d_alpha;
                p.reaction_rate = d_alpha / dt;
                // 新生成气体带入原固相显热，反应热只加入气相内能。
                p.gas_internal_energy += p.mass * d_alpha
                    * (p.specific_heat * p.temperature + p.reaction_heat);
                if (old_alpha <= 0.0 && d_alpha > 0.0)
                    p.phase = 1;
            }

            const double solid_fraction = std::max(0.0, 1.0 - p.reaction_progress);
            if (solid_fraction > 1.0e-12)
            {
                const double heat_capacity = p.mass * p.specific_heat * solid_fraction;
                p.temperature += energy_delta[i] / heat_capacity;
            }
            else
            {
                // 燃尽点已经没有固相热容；连接传热/外热源应进入气相内能，
                // 避免用一个人为的极小固相热容造成温度发散。
                p.gas_internal_energy += energy_delta[i];
            }
            if (!std::isfinite(p.temperature))
                throw ZaranError("DEM thermal/reaction update produced a non-finite temperature");
            p.temperature = std::max(1.0, p.temperature);
            if (GetDEMParam()->GetMaxTemperature() > 0.0)
                p.temperature = std::min(p.temperature, GetDEMParam()->GetMaxTemperature());

            if (p.reaction_progress >= GetDEMParam()->GetBurnoutThreshold()) p.phase = 2;
            else if (p.reaction_progress >= GetDEMParam()->GetMixedPhaseThreshold()) p.phase = 1;
            else p.phase = 0;

            UpdateGasSolidState(p);
        }
    }

    void DEMSolver::CalcGasPressureForce()
    {
        if (!GetDEMParam()->GetReactionEnabled()) return;
        auto& particles = m_dem_data->GetParticles();
    auto apply_pressure_face = [&](index_type index_a, index_type index_b,
                                   double face_area, double reference_length)
    {
        if (face_area <= 0.0) return;
        DEMParticle& pa = particles[index_a];
        DEMParticle& pb = particles[index_b];
        const double gas_fraction = std::max(pa.reaction_progress, pb.reaction_progress);
        if (gas_fraction <= 0.0) return;
        const Eigen::Vector3d delta = pb.pos - pa.pos;
        const double length = delta.norm();
        if (length <= 1.0e-15) return;

        // LSM gas_phase_mod 使用 Voronoi 面扣除两固相核遮挡后的气相暴露面，
        // 并保证一个随反应度逐渐张开的最小面积。Zaran 当前没有动态 Voronoi 网，
        // 因而在原始连接面上采用同一最小开口比例作局部近似。
        const double exposed_fraction = std::min(1.0, std::max(0.0,
            2.0 * (1.0 - std::sqrt(std::max(0.0, 1.0 - gas_fraction)))));
        const double pressure_area = face_area * exposed_fraction;
        const double pressure = 0.5 * (pa.gas_pressure + pb.gas_pressure);
        // LSM 还以参考格距/当前格距修正气相排斥力。
        const Eigen::Vector3d force = pressure * pressure_area
                                    * (reference_length / length) * delta / length;
        pa.force -= force;
        pb.force += force;
    };
    if (m_voronoi_valid)
    {
        for (const auto& face : m_gas_voronoi_faces)
        {
            const double distance = (particles[face.idx_b].pos - particles[face.idx_a].pos).norm();
            apply_pressure_face(face.idx_a, face.idx_b, face.area,
                m_lattice_spacing > 0.0 ? m_lattice_spacing : distance);
        }
    }
    else
    {
        for (const auto& bond : m_dem_data->GetBonds())
            apply_pressure_face(bond.idx_a, bond.idx_b, bond.conduction_area, bond.rest_length);
        }
    }

    void DEMSolver::CalcBondForce()
    {
        auto& particles = m_dem_data->GetParticles();
        const double dt = GetDEMParam()->GetTimeStep();
        for (auto& bond : m_dem_data->GetBonds())
        {
            if (!bond.active) continue;
            DEMParticle& pa = particles[bond.idx_a];
            DEMParticle& pb = particles[bond.idx_b];
            const double solid_a = std::max(0.0, 1.0 - pa.reaction_progress);
            const double solid_b = std::max(0.0, 1.0 - pb.reaction_progress);
            if (solid_a <= 1.0 - GetDEMParam()->GetBurnoutThreshold()
                || solid_b <= 1.0 - GetDEMParam()->GetBurnoutThreshold())
            {
                // 与 LSM 的固相格点燃尽删除一致：任一端燃尽后，该连接不再承载固相力。
                // 固定连接仍保留为当前近似气相拓扑，供压力和燃烧界面使用。
                bond.active = false;
                bond.force_a.setZero();
                continue;
            }
            const Eigen::Vector3d separation = pb.pos - pa.pos;
            const double length = separation.norm();
            if (length <= 1.0e-15)
                throw ZaranError("A DEM bond has coincident endpoints");

            // LSM 通过缩小固相核几何尺寸退化承载截面，而不是直接按反应度降低
            // 材料的断裂应变。几何平均保证端点交换对称。
            const double solid_area_scale = std::sqrt(solid_a * solid_b);
            const double limit = GetDEMParam()->GetBondBreakStrain();
            const double peak = GetDEMParam()->GetBondPeakStrain();
            const double tensile_strain = std::max(0.0, (length - bond.rest_length) / bond.rest_length);
            bond.maximum_tensile_strain = std::max(bond.maximum_tensile_strain, tensile_strain);

            const Eigen::Vector3d normal = separation / length;
            const Eigen::Vector3d arm_a = 0.5 * length * normal;
            const Eigen::Vector3d arm_b = -arm_a;
            Eigen::Vector3d relative_velocity = pa.vel + pa.omega.cross(arm_a)
                - pb.vel - pb.omega.cross(arm_b);
            const double velocity_n = relative_velocity.dot(normal);
            const Eigen::Vector3d velocity_t = relative_velocity - velocity_n * normal;

            bond.delta_t += velocity_t * dt;
            bond.delta_t -= bond.delta_t.dot(normal) * normal;
            bond.extension = length - bond.rest_length;

            bond.elastic_energy = 0.5 * solid_area_scale
                * (bond.normal_stiffness * bond.extension * bond.extension
                   + bond.tangential_stiffness * bond.delta_t.squaredNorm());
            if (bond.fracture_energy > 0.0 && bond.extension >= 0.0
                && bond.elastic_energy >= bond.fracture_energy * solid_area_scale)
            {
                // LSM LEP：法向与切向储能之和超过表面能控制的单键断裂能。
                bond.dissipated_fracture_energy += bond.elastic_energy;
                bond.damage = 1.0;
                bond.active = false;
                bond.force_a.setZero();
                continue;
            }
            if (bond.fracture_energy <= 0.0)
            {
                if (bond.maximum_tensile_strain > peak && limit > peak)
                    bond.damage = std::max(bond.damage,
                        std::min(1.0, (bond.maximum_tensile_strain - peak) / (limit - peak)));
                if (tensile_strain > limit)
                {
                    bond.dissipated_fracture_energy += bond.elastic_energy;
                    bond.active = false;
                    bond.force_a.setZero();
                    continue;
                }
            }

            const double effective_stiffness = (1.0 - bond.damage) * solid_area_scale;
            double elastic_normal_force = bond.normal_stiffness * bond.extension;
            if (bond.extension < 0.0 && GetDEMParam()->GetHighPressureKnn() > 0.0
                && GetDEMParam()->GetHighPressureRnn() > 0.0)
            {
                // LSM LEP/REP 高压压缩支路：在 r/r0=rnn 处切线刚度为 knn*Kn0。
                const double exponent = GetDEMParam()->GetHighPressureExponent();
                const double compression_scale = (GetDEMParam()->GetHighPressureRnn() - 1.0)
                                               * bond.rest_length;
                const double normalized_compression = bond.extension / compression_scale;
                const double coefficient = GetDEMParam()->GetHighPressureKnn()
                    * bond.normal_stiffness * compression_scale / (exponent + 1.0);
                elastic_normal_force = coefficient
                    * std::pow(std::max(0.0, normalized_compression), exponent + 1.0);
            }
            const Eigen::Vector3d force_n = (effective_stiffness * elastic_normal_force
                - bond.normal_damping * velocity_n) * normal;
            const Eigen::Vector3d force_t = -effective_stiffness * bond.tangential_stiffness * bond.delta_t
                - bond.tangential_damping * velocity_t;
            bond.force_a = force_n + force_t;

            pa.force += bond.force_a;
            pa.torque += arm_a.cross(force_t);
            pb.force -= bond.force_a;
            pb.torque += arm_b.cross(-force_t);
        }
    }

    DEMSolver::ContactKey DEMSolver::MakeContactKey(const DEMContact& contact) const
    {
        return std::make_tuple(static_cast<int>(contact.type), contact.idx_a, contact.idx_b);
    }

    void DEMSolver::RestoreContactHistory(DEMContact& contact) const
    {
        const auto it = m_tangential_history.find(MakeContactKey(contact));
        if (it != m_tangential_history.end())
            contact.delta_t = it->second;
    }

    void DEMSolver::SaveContactHistory(const DEMContact& contact)
    {
        const ContactKey key = MakeContactKey(contact);
        m_tangential_history[key] = contact.delta_t;
        m_contacts_this_step.insert(key);
    }

    void DEMSolver::PruneContactHistory()
    {
        for (auto it = m_tangential_history.begin(); it != m_tangential_history.end();)
        {
            if (m_contacts_this_step.find(it->first) == m_contacts_this_step.end())
                it = m_tangential_history.erase(it);
            else
                ++it;
        }
    }

    void DEMSolver::ZeroForce()
    {
        for (auto& p : m_dem_data->GetParticles())
        {
            p.force.setZero();
            p.torque.setZero();
        }
    }

    void DEMSolver::ContactDetection()
    {
        auto& particles = m_dem_data->GetParticles();
        const index_type N = particles.size();

        m_dem_data->ClearContacts();

        if (N == 0) return;

        // 固定粒子也进入邻域搜索；它们不积分，但仍可向活动粒子施加接触力。
        point_vec pts;
        pts.reserve(N);
        for (index_type i = 0; i < N; ++i)
        {
            pts.push_back({ particles[i].pos.x(), particles[i].pos.y(), particles[i].pos.z() });
        }

        if (N < 2) return;

        KDTree kd(pts);

        // 对每个粒子 i，查询半径 = r_i + r_max_neighbor 内的邻居
        double r_max = 0.0;
        for (const auto& particle : particles)
            r_max = std::max(r_max, particle.radius);

        for (index_type i = 0; i < N; ++i)
        {
            const DEMParticle& pa = particles[i];
            coord_vec pt_i = { pa.pos.x(), pa.pos.y(), pa.pos.z() };

            double search_r = pa.radius + r_max;
            auto neighbors = kd.NeighborhoodIndices(pt_i, search_r);

            for (index_type kj : neighbors)
            {
                index_type j = kj;
                if (j <= i) continue; // 每对只算一次

                const DEMParticle& pb = particles[j];
                if (!pa.IsDynamic() && !pb.IsDynamic()) continue;
                Eigen::Vector3d d = pb.pos - pa.pos;
                double dist = d.norm();
                double sum_r = pa.radius + pb.radius;
                if (dist >= sum_r || dist < 1.0e-15) continue;

                DEMContact c;
                c.type = ContactType::ParticleParticle;
                c.idx_a = i;
                c.idx_b = j;
                c.overlap_n = sum_r - dist;
                c.normal = d / dist; // n 从 A 指向 B
                c.contact_point = pa.pos + (pa.radius - 0.5 * c.overlap_n) * c.normal;
                RestoreContactHistory(c);
                m_dem_data->AddContact(c);
            }
        }
    }

    void DEMSolver::CalcContactForce()
    {
        auto& particles = m_dem_data->GetParticles();
        double dt = GetDEMParam()->GetTimeStep();

        for (auto& c : m_dem_data->GetContacts())
        {
            if (c.type != ContactType::ParticleParticle) continue;
            DEMParticle& pa = particles[c.idx_a];
            DEMParticle& pb = particles[c.idx_b];

            m_contact_model->CalcNormalForce(pa, pb, c, dt);
            m_contact_model->CalcTangentialForce(pa, pb, c, dt);
            SaveContactHistory(c);

            Eigen::Vector3d F = c.force_n + c.force_t;

            // 作用-反作用
            pa.force += F;
            pa.torque += c.contact_point.cross(c.force_t) - pa.pos.cross(c.force_t);
            pb.force -= F;
            pb.torque -= c.contact_point.cross(c.force_t) - pb.pos.cross(c.force_t);
        }
    }

    void DEMSolver::CalcWallForce()
    {
        auto& particles = m_dem_data->GetParticles();
        auto& walls = m_dem_data->GetWalls();
        double dt = GetDEMParam()->GetTimeStep();

        for (auto& wall : walls)
        {
            for (index_type pi = 0; pi < particles.size(); ++pi)
            {
                DEMParticle& pa = particles[pi];
                if (!pa.IsDynamic()) continue;
                double d = wall.SignedDist(pa.pos);
                double overlap = pa.radius - d;
                if (overlap <= 0.0) continue;

                // 构造虚拟墙粒子用于接触力计算
                DEMParticle pb_wall;
                pb_wall.radius = 1.0e10; // 极大半径 → 平面
                pb_wall.mass = 1.0e30;
                pb_wall.vel.setZero();
                pb_wall.omega.setZero();
                pb_wall.young_modulus = wall.young_modulus;
                pb_wall.poisson_ratio = wall.poisson_ratio;
                pb_wall.friction_coeff = wall.friction_coeff;
                pb_wall.restitution_coeff = wall.restitution_coeff;
                // 墙粒子质心沿法向偏移极大半径
                pb_wall.pos = pa.pos - (d - pa.radius + 1.0e10) * wall.normal;

                DEMContact c;
                c.type = ContactType::ParticleWall;
                c.idx_a = pi;
                c.idx_b = wall.id;
                c.overlap_n = overlap;
                c.normal = -wall.normal; // 法向从 A 指向墙
                c.contact_point = pa.pos - pa.radius * wall.normal;
                RestoreContactHistory(c);

                m_contact_model->CalcNormalForce(pa, pb_wall, c, dt);
                m_contact_model->CalcTangentialForce(pa, pb_wall, c, dt);
                SaveContactHistory(c);

                pa.force += c.force_n + c.force_t;
                pa.torque += c.contact_point.cross(c.force_t) - pa.pos.cross(c.force_t);
            }
        }
    }

    void DEMSolver::CalcGravity()
    {
        const Eigen::Vector3d& g = GetDEMParam()->GetGravity();
        for (auto& p : m_dem_data->GetParticles())
        {
            if (p.IsDynamic())
                p.force += p.mass * g;
        }
    }

    void DEMSolver::Integrate()
    {
        double dt = GetDEMParam()->GetTimeStep();
        for (auto& p : m_dem_data->GetParticles())
        {
            if (!p.active) continue;
            if (p.kinematic)
            {
                p.rotation += p.omega * dt;
                p.pos += p.vel * dt;
                continue;
            }
            // 半隐式 Euler：先更新速度，再用新速度更新位置。
            Eigen::Vector3d acc = p.force / p.mass;
            Eigen::Vector3d alpha = p.torque / (p.inertia * p.mass * p.radius * p.radius);

            p.vel += acc * dt;
            p.omega += alpha * dt;
            p.rotation += p.omega * dt;
            p.pos += p.vel * dt;

        }
    }

    void DEMSolver::BackupField(std::string& back_folder)
    {
        BackupField(static_cast<const std::string&>(back_folder));
    }

    void DEMSolver::BackupField(const std::string& back_folder) const
    {
        std::string fname = back_folder + "/particles.dat";
        std::ofstream fout(fname);
        if (!fout.is_open())
        {
            Log::warn("DEMSolver::BackupField: cannot open {}", fname);
            return;
        }
        fout << "variables=id,group,radius,mass,px,py,pz,vx,vy,vz,rx,ry,rz,ox,oy,oz,fx,fy,fz,motion_type,temperature,reaction_progress,reaction_rate,phase,energetic,gas_temperature,gas_pressure,volume_ratio,total_volume,solid_volume,gas_volume,solid_core_radius,gas_radius,gas_internal_energy,internal_heat_transfer,body_reaction_increment,core_burn_increment,neighbor_burn_increment\n";
        for (const auto& p : m_dem_data->GetParticles())
        {
            fout << p.id << ","
                << p.group << ","
                << p.radius << ","
                << p.mass << ","
                << p.pos.x() << "," << p.pos.y() << "," << p.pos.z() << ","
                << p.vel.x() << "," << p.vel.y() << "," << p.vel.z() << ","
                << p.rotation.x() << "," << p.rotation.y() << "," << p.rotation.z() << ","
                << p.omega.x() << "," << p.omega.y() << "," << p.omega.z() << ","
                << p.force.x() << "," << p.force.y() << "," << p.force.z() << ","
                << (p.kinematic ? 2 : (p.active ? 1 : 0)) << ","
                << p.temperature << "," << p.reaction_progress << ","
                << p.reaction_rate << "," << p.phase << "," << (p.energetic ? 1 : 0) << ","
                << p.gas_temperature << "," << p.gas_pressure << "," << p.volume_ratio << ","
                << p.total_volume << "," << p.solid_volume << "," << p.gas_volume << ","
                << p.solid_core_radius << "," << p.gas_radius << ","
                << p.gas_internal_energy << "," << p.internal_heat_transfer << ","
                << p.body_reaction_increment << "," << p.core_burn_increment << ","
                << p.neighbor_burn_increment << "\n";
        }
        fout.close();

        const std::string bond_name = back_folder + "/bonds.dat";
        std::ofstream bond_out(bond_name);
        if (!bond_out.is_open())
        {
            Log::warn("DEMSolver::BackupField: cannot open {}", bond_name);
            return;
        }
        bond_out << "id,particle_a_id,particle_b_id,rest_length,extension,fx,fy,fz,active,heat_flow_a,damage,maximum_tensile_strain,elastic_energy,fracture_energy,dissipated_fracture_energy\n";
        const auto& particles = m_dem_data->GetParticles();
        for (const auto& bond : m_dem_data->GetBonds())
        {
            bond_out << bond.id << "," << particles[bond.idx_a].id << ","
                << particles[bond.idx_b].id << "," << bond.rest_length << ","
                << bond.extension << "," << bond.force_a.x() << ","
                << bond.force_a.y() << "," << bond.force_a.z() << ","
                << (bond.active ? 1 : 0) << "," << bond.heat_flow_a << ","
                << bond.damage << "," << bond.maximum_tensile_strain << ","
                << bond.elastic_energy << "," << bond.fracture_energy << ","
                << bond.dissipated_fracture_energy << "\n";
        }
        Log::info("DEMSolver::BackupField: {} particles and {} bonds written to {}",
            m_dem_data->GetParticleNum(), m_dem_data->GetBondNum(), back_folder);

        if (m_voronoi_valid)
        {
            const std::string face_name = back_folder + "/gas_voronoi_faces.dat";
            std::ofstream face_out(face_name);
            if (face_out.is_open())
            {
                face_out << "particle_a_id,particle_b_id,length,area,phase_a,phase_b,pressure_a,pressure_b\n";
                for (const auto& face : m_gas_voronoi_faces)
                {
                    const auto& pa = particles[face.idx_a];
                    const auto& pb = particles[face.idx_b];
                    face_out << pa.id << "," << pb.id << "," << face.length << ","
                        << face.area << "," << pa.phase << "," << pb.phase << ","
                        << pa.gas_pressure << "," << pb.gas_pressure << "\n";
                }
            }
        }
    }

} // namespace zaran
