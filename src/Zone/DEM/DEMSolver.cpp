#include "DEMSolver.h"
#include "LinearSpringDashpot.h"
#include "HertzMindlin.h"
#include "FastNumberFormat.h"
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

/// @brief splitmix64：用于逐键 Weibull 抽样的确定性哈希。
/// 键与粒子 id 绑定（而不是数组下标），因此同一物理连接无论输入顺序如何都得到同一强度。
std::uint64_t SplitMix64(std::uint64_t x)
{
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

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
        // --- 可选：仿真开始时（t=0）按几何邻近关系建立弹簧连接网络 ---
        // 与逐步接触检测相互独立：仅在初始化阶段执行一次，之后每个时间步
        // 由 CalcBondForce 统一计算弹簧力与断裂判据。
        if (para->GetSpringNetworkEnabled())
            BuildSpringNetwork();

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
        // --- 可选：逐键强度异质性（Weibull）---
        // 放在键合参数（含 fracture_energy）确定之后、且覆盖文件给定与自动建链两种来源。
        AssignBondStrengthScales();
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
        ZeroForce();
        CalcBondForce();
        ContactDetection();
        CalcContactForce();
        CalcWallForce();
        AdvanceContactHistory();
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
        // 温度更新与反应解耦：
        //   反应 / 气相 / Voronoi / 燃烧 / Arrhenius 化学只在 reaction_enabled 时执行；
        //   而「键合导热 + 外部体热源 + 机械耗散（接触摩擦与阻尼、键合阻尼）」以及
        //   温度更新本身，对所有算例（含惰性算例）都执行。
        const bool reaction = GetDEMParam()->GetReactionEnabled();
        if (reaction) UpdateGasVoronoiMesh();
        auto& particles = m_dem_data->GetParticles();
        const double dt = GetDEMParam()->GetTimeStep();
        const index_type np = particles.size();
        // 复用成员缓冲（容量跨步保留），避免每步分配 7 个 N 大小数组
        auto& energy_delta        = m_tmp_energy_delta;
        auto& surface_delta       = m_tmp_surface_delta;
        auto& core_delta          = m_tmp_core_delta;
        auto& length_ratio_sum    = m_tmp_length_ratio_sum;
        auto& length_ratio_count  = m_tmp_length_ratio_count;
        auto& gas_pressure_sum    = m_tmp_gas_pressure_sum;
        auto& gas_pressure_count  = m_tmp_gas_pressure_count;
        energy_delta.assign(np, 0.0);
        surface_delta.assign(np, 0.0);
        core_delta.assign(np, 0.0);
        length_ratio_sum.assign(np, 0.0);
        length_ratio_count.assign(np, 0);
        gas_pressure_sum.assign(np, 0.0);
        gas_pressure_count.assign(np, 1);

        if (reaction)
        {
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

        if (reaction)
        {
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
        }

        for (index_type i = 0; i < particles.size(); ++i)
        {
            DEMParticle& p = particles[i];
            if (reaction)
            {
                p.reaction_rate = 0.0;
                p.body_reaction_increment = 0.0;
                p.core_burn_increment = 0.0;
                p.neighbor_burn_increment = 0.0;
            }
            // ==== 以下对所有算例（含惰性）都执行 ====
            // 外部体热源（来自粒子文件）
            energy_delta[i] += p.heat_source * dt;
            // 机械耗散生热：接触法向阻尼 + 接触切向摩擦 + 键合阻尼
            if (GetDEMParam()->GetMechanicalHeatingEnabled())
                energy_delta[i] += m_tmp_dissipation[i];

            if (reaction)
            {
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

            if (reaction) UpdateGasSolidState(p);
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

    void DEMSolver::AssignBondStrengthScales()
    {
        const double modulus = GetDEMParam()->GetBondWeibullModulus();
        if (!(modulus > 0.0))
            return; // 均质：所有键保持 strength_scale = 1

        // Weibull(k = modulus, λ) 的形状参数不为 1 时均值不等于 λ，
        // 取 λ = 1/Γ(1+1/k) 把均值归一到 1 —— 这样"名义阈值"仍是期望值，
        // 引入的只是离散度，不改变上一节按宏观量标定出的整体强度。
        const double lambda = 1.0 / std::tgamma(1.0 + 1.0 / modulus);
        const std::uint64_t seed = GetDEMParam()->GetBondWeibullSeed();
        const auto& particles = m_dem_data->GetParticles();

        double sum = 0.0;
        double minimum = 1.0e30;
        double maximum = -1.0e30;
        for (auto& bond : m_dem_data->GetBonds())
        {
            const std::uint64_t ia = static_cast<std::uint64_t>(particles[bond.idx_a].id);
            const std::uint64_t ib = static_cast<std::uint64_t>(particles[bond.idx_b].id);
            const std::uint64_t key = SplitMix64(seed
                ^ (ia * 0x9E3779B97F4A7C15ULL)
                ^ ((ib + 0x165667B19E3779F9ULL) * 0xC2B2AE3D27D4EB4FULL));
            // u ∈ (0,1)，避免 u = 0 时抽样退化为 0
            const double u = (static_cast<double>(key >> 11) + 0.5) * (1.0 / 9007199254740992.0);
            const double sample = lambda * std::pow(-std::log(1.0 - u), 1.0 / modulus);
            // 保护：极端尾部不应产生 0 或爆炸式阈值
            bond.strength_scale = std::min(std::max(sample, 1.0e-2), 1.0e2);
            sum += bond.strength_scale;
            minimum = std::min(minimum, bond.strength_scale);
            maximum = std::max(maximum, bond.strength_scale);
        }

        const double count = static_cast<double>(std::max<std::size_t>(1, m_dem_data->GetBondNum()));
        Log::info("DEMSolver: bond strength heterogeneity ON "
                  "(Weibull m={:g}, seed={}, {} bonds, mean scale={:.4f}, min={:.4f}, max={:.4f})",
                  modulus, seed, m_dem_data->GetBondNum(), sum / count, minimum, maximum);
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
                bond.elastic_energy = 0.0; // 不承载也不储能
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
            // 逐键强度折减（Weibull 异质性）：阈值按 strength_scale 缩放，
            // strength_scale = 1 时与均质情形完全一致。
            const double strength_scale = bond.strength_scale;
            const double limit = GetDEMParam()->GetBondBreakStrain() * strength_scale;
            const double peak = GetDEMParam()->GetBondPeakStrain() * strength_scale;
            const double tensile_strain = std::max(0.0, (length - bond.rest_length) / bond.rest_length);
            const double compressive_strain = std::max(0.0, (bond.rest_length - length) / bond.rest_length);
            bond.maximum_tensile_strain = std::max(bond.maximum_tensile_strain, tensile_strain);
            bond.maximum_compressive_strain = std::max(bond.maximum_compressive_strain, compressive_strain);

            const Eigen::Vector3d normal = separation / length;
            const Eigen::Vector3d arm_a = 0.5 * length * normal;
            const Eigen::Vector3d arm_b = -arm_a;
            Eigen::Vector3d relative_velocity = pa.vel + pa.omega.cross(arm_a)
                - pb.vel - pb.omega.cross(arm_b);
            const double velocity_n = relative_velocity.dot(normal);
            const Eigen::Vector3d velocity_t = relative_velocity - velocity_n * normal;

            // 键合阻尼耗散：阻尼力 -nd*vn*n - td*vt 对相对运动做负功，
            // 该对被耗散的功率为 nd*vn² + td*|vt|²（恒 ≥ 0），两端各记一半。
            {
                const double damping_energy = (bond.normal_damping * velocity_n * velocity_n
                    + bond.tangential_damping * velocity_t.squaredNorm()) * dt;
                if (damping_energy > 0.0)
                {
                    m_tmp_dissipation[bond.idx_a] += 0.5 * damping_energy;
                    m_tmp_dissipation[bond.idx_b] += 0.5 * damping_energy;
                }
            }

            bond.delta_t += velocity_t * dt;
            bond.delta_t -= bond.delta_t.dot(normal) * normal;
            bond.extension = length - bond.rest_length;

            bond.elastic_energy = 0.5 * solid_area_scale
                * (bond.normal_stiffness * bond.extension * bond.extension
                   + bond.tangential_stiffness * bond.delta_t.squaredNorm());

            // --- 压缩/屈曲失效（脆性材料在压缩侧同样有强度上限）---
            // 线性弹性下压缩强度 σ_c 对应应变 ε_c = σ_c/E，而抗拉阈值 ε_t = σ_t/E，
            // 故 ε_c = (σ_c/σ_t)·ε_t：默认直接由抗拉阈值按强度比派生，
            // 也可用 dem.bond_break_strain_compression 直接给定。
            // 与拉伸的"先软化后断"不同，压缩破坏取脆性判据（达阈值即断），
            // 因为粉碎区的渐进失效正来自应变场的非均匀性而非软化段。
            // 必须独立于断裂能准则：能量准则只在 extension ≥ 0 时生效（LSM LEP），
            // 若把压缩判据也放进那个分支，能量口径的算例会失去压缩强度。
            const double compression_limit = GetDEMParam()->GetBondBreakStrainCompression() > 0.0
                ? GetDEMParam()->GetBondBreakStrainCompression() * strength_scale
                : GetDEMParam()->GetBondCompressionStrengthRatio() * limit;
            if (compression_limit > 0.0 && compressive_strain > compression_limit)
            {
                bond.dissipated_fracture_energy += bond.elastic_energy;
                // 断开的键不再储能：清零后 Σelastic_energy 只统计完好键，
                // 不会与 dissipated_fracture_energy 重复计数（否则能量账会算错）。
                bond.elastic_energy = 0.0;
                bond.damage = 1.0;
                bond.active = false;
                bond.force_a.setZero();
                continue;
            }

            if (bond.fracture_energy > 0.0 && bond.extension >= 0.0
                && bond.elastic_energy >= bond.fracture_energy * solid_area_scale
                    * strength_scale * strength_scale)
            {
                // LSM LEP：法向与切向储能之和超过表面能控制的单键断裂能。
                // 强度按 s 缩放时应力应变同乘 s，故断裂能按 s² 缩放，保持口径一致。
                bond.dissipated_fracture_energy += bond.elastic_energy;
                // 断开的键不再储能：清零后 Σelastic_energy 只统计完好键，
                // 不会与 dissipated_fracture_energy 重复计数（否则能量账会算错）。
                bond.elastic_energy = 0.0;
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
                    bond.elastic_energy = 0.0; // 断开的键不再储能
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

    DEMSolver::ContactKey DEMSolver::MakeContactKey(const DEMContact& contact)
    {
        // 打包为 64 位：粒子-粒子用 (min<<32)|max；粒子-墙用最高位置 1。
        if (contact.type == ContactType::ParticleWall)
        {
            return (std::uint64_t(1) << 63)
                 | (std::uint64_t(contact.idx_a) << 32)
                 | std::uint64_t(contact.idx_b & 0xffffffffu);
        }
        index_type a = contact.idx_a;
        index_type b = contact.idx_b;
        if (a > b) std::swap(a, b);
        return (std::uint64_t(a) << 32) | std::uint64_t(b & 0xffffffffu);
    }

    void DEMSolver::RestoreContactHistory(DEMContact& contact) const
    {
        if (const Eigen::Vector3d* v = m_history_prev.Find(MakeContactKey(contact)))
            contact.delta_t = *v;
    }

    void DEMSolver::SaveContactHistory(const DEMContact& contact)
    {
        m_history_cur.Insert(MakeContactKey(contact), contact.delta_t);
    }

    void DEMSolver::AdvanceContactHistory()
    {
        // 当前表成为下一步的历史；另一张表 O(1) 清空后作为新表。
        std::swap(m_history_prev, m_history_cur);
        m_history_cur.Reset();
    }

    void DEMSolver::ZeroForce()
    {
        auto& particles = m_dem_data->GetParticles();
        // 机械耗散缓冲：容量跨步保留，每步清零（由键合/接触/墙面三处累加）
        m_tmp_dissipation.assign(particles.size(), 0.0);
        for (auto& p : particles)
        {
            p.force.setZero();
            p.torque.setZero();
        }
    }

    double DEMSolver::BuildUniformGrid(const dynamic_array<DEMParticle>& particles, double cell)
    {
        const index_type N = particles.size();

        // 包围盒
        Eigen::Vector3d lo = particles[0].pos;
        Eigen::Vector3d hi = lo;
        for (index_type i = 0; i < N; ++i)
        {
            lo = lo.cwiseMin(particles[i].pos);
            hi = hi.cwiseMax(particles[i].pos);
        }
        m_grid_lo = lo;

        if (!(cell > 0.0)) cell = 1.0; // 退化保护（半径为 0）
        const Eigen::Vector3d ext = (hi - lo).cwiseMax(Eigen::Vector3d::Zero());

        // 单元总数上限：稀疏/超大域时放大 cell（放大只会增加候选，不影响正确性）
        const double k_max_cells = 1.0e6;
        index_type nx = 1, ny = 1, nz = 1;
        for (int guard = 0; guard < 64; ++guard)
        {
            nx = static_cast<index_type>(std::floor(ext.x() / cell)) + 1;
            ny = static_cast<index_type>(std::floor(ext.y() / cell)) + 1;
            nz = static_cast<index_type>(std::floor(ext.z() / cell)) + 1;
            if (double(nx) * double(ny) * double(nz) <= k_max_cells) break;
            cell *= 2.0;
        }
        m_grid_nx  = nx;
        m_grid_ny  = ny;
        m_grid_nz  = nz;
        m_grid_nyz = ny * nz;
        m_grid_cell = cell;

        // 链表式单元表（O(N)，缓冲跨步复用，无每步分配）
        const index_type ncell = nx * m_grid_nyz;
        m_grid_head.assign(ncell, static_cast<index_type>(-1));
        m_grid_next.resize(N);
        for (index_type i = 0; i < N; ++i)
        {
            const index_type c = GridCellOf(particles[i].pos);
            m_grid_next[i] = m_grid_head[c];
            m_grid_head[c] = i;
        }
        return cell;
    }

    index_type DEMSolver::GridCellOf(const Eigen::Vector3d& q) const
    {
        index_type ix = static_cast<index_type>((q.x() - m_grid_lo.x()) / m_grid_cell);
        index_type iy = static_cast<index_type>((q.y() - m_grid_lo.y()) / m_grid_cell);
        index_type iz = static_cast<index_type>((q.z() - m_grid_lo.z()) / m_grid_cell);
        if (ix >= m_grid_nx) ix = m_grid_nx - 1;
        if (iy >= m_grid_ny) iy = m_grid_ny - 1;
        if (iz >= m_grid_nz) iz = m_grid_nz - 1;
        return ix * m_grid_nyz + iy * m_grid_nz + iz;
    }

    void DEMSolver::ContactDetection()
    {
        auto& particles = m_dem_data->GetParticles();
        const index_type N = particles.size();

        m_dem_data->ClearContacts();

        if (N < 2) return;

        // ---- 包围盒与最大半径 ----
        // 接触判据：dist < r_i + r_j ≤ 2·r_max。取单元边长 ≥ 2·r_max，
        // 则任一接触对必落在相邻（含自身）的 3×3×3 单元内。
        double r_max = 0.0;
        for (index_type i = 0; i < N; ++i)
            r_max = std::max(r_max, particles[i].radius);
        BuildUniformGrid(particles, 2.0 * r_max);

        // ---- 3×3×3 邻域配对（每个单元只访问一次，天然无重复对）----
        for (index_type i = 0; i < N; ++i)
        {
            const DEMParticle& pa = particles[i];
            const index_type ci = GridCellOf(pa.pos);
            const index_type ix0 = ci / m_grid_nyz;
            const index_type iy0 = (ci / m_grid_nz) % m_grid_ny;
            const index_type iz0 = ci % m_grid_nz;

            const index_type xb = (ix0 == 0) ? 0 : ix0 - 1;
            const index_type xe = std::min(ix0 + 1, m_grid_nx - 1);
            const index_type yb = (iy0 == 0) ? 0 : iy0 - 1;
            const index_type ye = std::min(iy0 + 1, m_grid_ny - 1);
            const index_type zb = (iz0 == 0) ? 0 : iz0 - 1;
            const index_type ze = std::min(iz0 + 1, m_grid_nz - 1);

            for (index_type ix = xb; ix <= xe; ++ix)
            for (index_type iy = yb; iy <= ye; ++iy)
            for (index_type iz = zb; iz <= ze; ++iz)
            {
                const index_type base = ix * m_grid_nyz + iy * m_grid_nz + iz;
                for (index_type j = m_grid_head[base]; j != static_cast<index_type>(-1); j = m_grid_next[j])
                {
                    if (j <= i) continue; // 每对只算一次
                    const DEMParticle& pb = particles[j];
                    const bool a_dyn = pa.IsDynamic();
                    const bool b_dyn = pb.IsDynamic();
                    if (!a_dyn && !b_dyn) continue;

                    const double sum_r = pa.radius + pb.radius;
                    const Eigen::Vector3d delta = pb.pos - pa.pos;

                    // --- 刚性边界（可选）：把 prescribed-motion 粒子当作刚体平面 ---
                    // 球-球法向由两球心连线给出，一旦粒子中心越过压板球心，法向翻转，
                    // 接触力由"压入"变成"抛出"，粒子会被弹出边界（表现为穿透）。
                    // 规定运动体的推进方向是确定的，故取其作为固定法向，并用投影长度
                    // 计算重叠量：这样即使深度侵入，力仍始终把粒子推回试件内侧。
                    Eigen::Vector3d normal;
                    double gap = 0.0;
                    bool use_rigid = false;
                    if (GetDEMParam()->GetRigidBoundaryEnabled() && (a_dyn != b_dyn))
                    {
                        const DEMParticle& kine = a_dyn ? pb : pa;
                        const double speed2 = kine.vel.squaredNorm();
                        if (speed2 > 1.0e-24)
                        {
                            const Eigen::Vector3d vhat = kine.vel / std::sqrt(speed2);
                            // normal 由 A 指向 B：B 为刚体时取 -v̂，A 为刚体时取 +v̂
                            normal = b_dyn ? vhat : -vhat;
                            gap = delta.dot(normal);
                            use_rigid = true;
                        }
                    }
                    if (!use_rigid)
                    {
                        const double dist2 = delta.squaredNorm();
                        if (dist2 >= sum_r * sum_r) continue;
                        const double dist = std::sqrt(dist2);
                        if (dist < 1.0e-15) continue;
                        normal = delta / dist; // n 从 A 指向 B
                        gap = dist;
                    }

                    const double overlap = sum_r - gap;
                    if (overlap <= 0.0) continue;

                    DEMContact c;
                    c.type = ContactType::ParticleParticle;
                    c.idx_a = i;
                    c.idx_b = j;
                    c.overlap_n = overlap;
                    c.normal = normal;
                    c.contact_point = pa.pos + (pa.radius - 0.5 * overlap) * normal;
                    RestoreContactHistory(c);
                    m_dem_data->AddContact(c);
                }
            }
        }

        // 预留接触历史容量（含墙面接触余量），避免写入过程中的再散列
        m_history_cur.Reserve(m_dem_data->GetContacts().size() * 2 + 64);
    }

    void DEMSolver::BuildSpringNetwork()
    {
        auto* para = GetDEMParam();
        const auto& particles = m_dem_data->GetParticles();
        const index_type N = particles.size();
        if (N < 2) return;

        // 已存在的连接对（一般来自 bonds.csv），避免重复建链
        std::set<std::pair<index_type, index_type>> existing;
        index_type next_bond_id = 0;
        for (const auto& b : m_dem_data->GetBonds())
        {
            index_type a = b.idx_a;
            index_type c = b.idx_b;
            if (a > c) std::swap(a, c);
            existing.insert({a, c});
            next_bond_id = std::max(next_bond_id, b.id + 1);
        }

        // 建链判据：d ≤ (r_a + r_b) · (1 + gap)
        const double gap = para->GetSpringNetworkGap();
        const double factor = 1.0 + gap;

        double r_max = 0.0;
        for (index_type i = 0; i < N; ++i)
            r_max = std::max(r_max, particles[i].radius);

        // 单元边长需 ≥ 最大连接距离 2·r_max·(1+gap)，保证 3×3×3 邻域即可覆盖
        BuildUniformGrid(particles, 2.0 * r_max * factor);

        const double user_kn = para->GetSpringNetworkStiffness();
        const double user_kt = para->GetSpringNetworkTangentialStiffness();
        const double fracture_strain = para->GetSpringNetworkFractureStrain();
        const bool energy_from_surface = para->GetSurfaceEnergy() > 0.0;

        index_type created = 0;
        const index_type npos = static_cast<index_type>(-1);

        for (index_type i = 0; i < N; ++i)
        {
            const DEMParticle& pa = particles[i];
            if (!pa.IsDynamic() && !pa.kinematic) continue;

            const index_type ci = GridCellOf(pa.pos);
            const index_type ix0 = ci / m_grid_nyz;
            const index_type iy0 = (ci / m_grid_nz) % m_grid_ny;
            const index_type iz0 = ci % m_grid_nz;
            const index_type xb = (ix0 == 0) ? 0 : ix0 - 1;
            const index_type xe = std::min(ix0 + 1, m_grid_nx - 1);
            const index_type yb = (iy0 == 0) ? 0 : iy0 - 1;
            const index_type ye = std::min(iy0 + 1, m_grid_ny - 1);
            const index_type zb = (iz0 == 0) ? 0 : iz0 - 1;
            const index_type ze = std::min(iz0 + 1, m_grid_nz - 1);

            for (index_type ix = xb; ix <= xe; ++ix)
            for (index_type iy = yb; iy <= ye; ++iy)
            for (index_type iz = zb; iz <= ze; ++iz)
            {
                const index_type base = ix * m_grid_nyz + iy * m_grid_nz + iz;
                for (index_type j = m_grid_head[base]; j != npos; j = m_grid_next[j])
                {
                    if (j <= i) continue; // 每对只建一次
                    const DEMParticle& pb = particles[j];
                    if (!pb.IsDynamic() && !pb.kinematic) continue;

                    const double dist = (pb.pos - pa.pos).norm();
                    if (dist <= 1.0e-15) continue;
                    if (dist > (pa.radius + pb.radius) * factor) continue;
                    if (!existing.insert({i, j}).second) continue; // 已有连接

                    DEMBond bond;
                    bond.id = next_bond_id++;
                    bond.idx_a = i;
                    bond.idx_b = j;
                    bond.rest_length = dist; // 初始无应力长度 = 初始中心距
                    bond.active = true;
                    bond.source = 1; // 标记为"初始时刻自动生成的弹簧连接网络"

                    if (user_kn > 0.0)
                    {
                        bond.normal_stiffness = user_kn;
                    }
                    else
                    {
                        // 与线性接触模型一致的等效刚度：k_n = 2·E*·R*
                        const double inv_E = (1.0 - pa.poisson_ratio * pa.poisson_ratio) / pa.young_modulus
                                           + (1.0 - pb.poisson_ratio * pb.poisson_ratio) / pb.young_modulus;
                        const double E_star = 1.0 / inv_E;
                        const double R_star = (pa.radius * pb.radius) / (pa.radius + pb.radius);
                        bond.normal_stiffness = 2.0 * E_star * R_star;
                    }
                    bond.tangential_stiffness = user_kt > 0.0 ? user_kt : 0.5 * bond.normal_stiffness;

                    // 断裂能：优先由表面能（InitField 中统一覆盖）；否则由断裂应变换算：
                    // 轴向弹性能 0.5·k_n·(ε·L0)² 达到该阈值即断。
                    if (!energy_from_surface && fracture_strain > 0.0)
                    {
                        const double critical_extension = fracture_strain * bond.rest_length;
                        bond.fracture_energy = 0.5 * bond.normal_stiffness
                                             * critical_extension * critical_extension;
                    }

                    m_dem_data->AddBond(bond);
                    ++created;
                }
            }
        }

        Log::info("DEMSolver: spring network built at t=0 -> {} connections "
                  "({} pre-existing pairs, {} bonds total, gap={:g}, fracture_strain={:g})",
                  created, existing.size() - created, m_dem_data->GetBondNum(), gap, fracture_strain);
        if (created == 0)
            Log::warn("DEMSolver: spring network enabled but no connection was generated "
                      "(check dem.spring_network_gap / particle spacing)");
    }

    bool DEMSolver::ApplyContactRebound(DEMParticle& pa, DEMParticle& pb,
                                        const DEMContact& contact) const
    {
        const double ratio = GetDEMParam()->GetContactReboundRatio();
        if (ratio <= 0.0) return false;

        const bool a_dyn = pa.IsDynamic();
        const bool b_dyn = pb.IsDynamic();
        if (!a_dyn && !b_dyn) return false;

        const double sum_r = pa.radius + pb.radius;
        if (sum_r <= 0.0) return false;

        // 接触距离。几何法向下 dist = (r_a+r_b) − δ；刚性边界模式（rigid_boundary）下
        // δ 是沿规定运动方向的投影重叠量，此处同样按该口径折算，用作统一的接近度指标。
        const double dist = sum_r - contact.overlap_n;
        if (dist >= ratio * sum_r) return false; // 尚未过深，不干预

        const Eigen::Vector3d& n = contact.normal; // A → B
        // 法向相对速度：< 0 表示正在相向接近，才需要回弹
        const double v_n = (pb.vel - pa.vel).dot(n);
        if (v_n >= 0.0) return false;

        // 目标：把法向相对速度由 v_n 改为 −e·v_n（即按恢复系数 e 反弹），
        // 再按线动量守恒把该增量分配给两端。
        double e = std::min(pa.restitution_coeff, pb.restitution_coeff);
        e = std::max(e, 1.0e-3);
        const double dv_n = -e * v_n - v_n; // > 0

        // 由 m_a·α_a + m_b·α_b = 0 与 (α_b − α_a)·n = Δv_n 解得：
        //   α_a = −Δv_n·m_b/(m_a+m_b)，α_b = +Δv_n·m_a/(m_a+m_b)
        if (!b_dyn)
        {
            // B 为规定运动体（无限质量）：A 独立承担全部修正，B 速度不变
            pa.vel += n * (-dv_n);
        }
        else if (!a_dyn)
        {
            pb.vel += n * dv_n;
        }
        else
        {
            const double msum = pa.mass + pb.mass;
            if (msum <= 0.0) return false;
            pa.vel += n * (-dv_n * pb.mass / msum);
            pb.vel += n * ( dv_n * pa.mass / msum);
        }
        return true;
    }

    void DEMSolver::AmplifyDeepOverlapForce(DEMContact& contact, double sum_radius) const
    {
        const double ratio = GetDEMParam()->GetContactStiffenRatio();
        if (ratio <= 0.0) return;
        if (sum_radius <= 0.0) return;

        const double limit = ratio * sum_radius;
        const double dist = sum_radius - contact.overlap_n;
        if (dist >= limit) return;

        // 参考 LSM：fn = (rn_limit/rn)^20 · fn。
        // 配合回弹阈值（默认 0.5），放大倍数实际有界（≤(0.6/0.5)^20≈38）；
        // 仍加上限保护，避免两端几乎完全重合时出现 inf/NaN。
        double factor = std::pow(limit / std::max(dist, 1.0e-30), 20.0);
        if (!(factor > 1.0)) return; // 含 NaN/未放大
        if (factor > 1.0e6) factor = 1.0e6;
        contact.force_n *= factor;
        // 法向阻尼耗散同步放大（此时 dissipation 只含法向阻尼项，切向尚未计算）
        contact.dissipation *= factor;
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

            // --- 重叠限制（① 刚性回弹）---
            // 深压缩且仍相向接近时，先把法向速度按弹性碰撞改写。
            // 放在法向力之前，使阻尼项使用回弹后的相对速度（与 LSM 的处理顺序一致）。
            ApplyContactRebound(pa, pb, c);

            m_contact_model->CalcNormalForce(pa, pb, c, dt);

            // --- 重叠限制（② 过深压缩时放大法向排斥力）---
            // 必须在切向之前：切向库仑摩擦上限取 contact.force_n 的模，
            // 放大后再算切向，可使摩擦上限同步放大（与 LSM 一致）。
            AmplifyDeepOverlapForce(c, pa.radius + pb.radius);

            m_contact_model->CalcTangentialForce(pa, pb, c, dt);
            SaveContactHistory(c);

            // 机械耗散（法向阻尼 + 切向摩擦）对半分给两端 → 后续计入温度
            const double dissipated = 0.5 * c.dissipation;
            m_tmp_dissipation[c.idx_a] += dissipated;
            m_tmp_dissipation[c.idx_b] += dissipated;

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

                // 墙面不可动（无限质量），该对的耗散能量全部计入粒子
                m_tmp_dissipation[pi] += c.dissipation;

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

    namespace
    {
        /// @brief 一次性写入文本文件（内容已在内存拼好），避免逐字段 << 的开销
        void WriteTextFile(const std::string& path, const std::string& content)
        {
            std::ofstream f(path, std::ios::binary);
            if (!f.is_open())
            {
                Log::warn("DEMSolver::BackupField: cannot open {}", path);
                return;
            }
            f.write(content.data(), static_cast<std::streamsize>(content.size()));
        }
    } // namespace

    void DEMSolver::BackupField(const std::string& back_folder) const
    {
        const auto& particles = m_dem_data->GetParticles();
        const auto& bonds     = m_dem_data->GetBonds();
        std::string& buf = m_out_buffer;

        // --- particles.dat ---
        // 缺陷修复：此前的表头以 "variables=" 开头（沿用了残差文件的 Tecplot 风格），
        // 与同目录的 bonds.dat、gas_voronoi_faces.dat 以及输入文件 particles.csv
        // 的纯 CSV 表头不一致，且会破坏下游按 csv.DictReader 解析 "id" 的脚本。
        // 现统一为纯 CSV 表头。
        //
        // 性能：改为"整块字符串缓冲 + std::to_chars + 单次 write"。
        // 旧实现逐字段 `fout <<`，万级粒子下每帧需数十万次流格式化（实测 ~200 ms/帧）。
        buf.clear();
        buf.reserve(particles.size() * 320 + 512);
        buf.append("id,group,radius,mass,px,py,pz,vx,vy,vz,rx,ry,rz,ox,oy,oz,fx,fy,fz,motion_type,temperature,reaction_progress,reaction_rate,phase,energetic,gas_temperature,gas_pressure,volume_ratio,total_volume,solid_volume,gas_volume,solid_core_radius,gas_radius,gas_internal_energy,internal_heat_transfer,body_reaction_increment,core_burn_increment,neighbor_burn_increment\n");
        for (const auto& p : particles)
        {
            AppendNumber(buf, p.id);                                    buf.push_back(',');
            AppendNumber(buf, p.group);                                 buf.push_back(',');
            AppendNumber(buf, p.radius);                                buf.push_back(',');
            AppendNumber(buf, p.mass);                                  buf.push_back(',');
            AppendNumber(buf, p.pos.x()); buf.push_back(',');
            AppendNumber(buf, p.pos.y()); buf.push_back(',');
            AppendNumber(buf, p.pos.z());                               buf.push_back(',');
            AppendNumber(buf, p.vel.x()); buf.push_back(',');
            AppendNumber(buf, p.vel.y()); buf.push_back(',');
            AppendNumber(buf, p.vel.z());                               buf.push_back(',');
            AppendNumber(buf, p.rotation.x()); buf.push_back(',');
            AppendNumber(buf, p.rotation.y()); buf.push_back(',');
            AppendNumber(buf, p.rotation.z());                          buf.push_back(',');
            AppendNumber(buf, p.omega.x()); buf.push_back(',');
            AppendNumber(buf, p.omega.y()); buf.push_back(',');
            AppendNumber(buf, p.omega.z());                             buf.push_back(',');
            AppendNumber(buf, p.force.x()); buf.push_back(',');
            AppendNumber(buf, p.force.y()); buf.push_back(',');
            AppendNumber(buf, p.force.z());                             buf.push_back(',');
            AppendNumber(buf, p.kinematic ? 2 : (p.active ? 1 : 0));     buf.push_back(',');
            AppendNumber(buf, p.temperature);                           buf.push_back(',');
            AppendNumber(buf, p.reaction_progress);                     buf.push_back(',');
            AppendNumber(buf, p.reaction_rate);                         buf.push_back(',');
            AppendNumber(buf, p.phase);                                 buf.push_back(',');
            AppendNumber(buf, p.energetic ? 1 : 0);                     buf.push_back(',');
            AppendNumber(buf, p.gas_temperature);                       buf.push_back(',');
            AppendNumber(buf, p.gas_pressure);                          buf.push_back(',');
            AppendNumber(buf, p.volume_ratio);                          buf.push_back(',');
            AppendNumber(buf, p.total_volume);                          buf.push_back(',');
            AppendNumber(buf, p.solid_volume);                          buf.push_back(',');
            AppendNumber(buf, p.gas_volume);                            buf.push_back(',');
            AppendNumber(buf, p.solid_core_radius);                     buf.push_back(',');
            AppendNumber(buf, p.gas_radius);                            buf.push_back(',');
            AppendNumber(buf, p.gas_internal_energy);                   buf.push_back(',');
            AppendNumber(buf, p.internal_heat_transfer);                buf.push_back(',');
            AppendNumber(buf, p.body_reaction_increment);               buf.push_back(',');
            AppendNumber(buf, p.core_burn_increment);                   buf.push_back(',');
            AppendNumber(buf, p.neighbor_burn_increment);               buf.push_back('\n');
        }
        WriteTextFile(back_folder + "/particles.dat", buf);

        // --- bonds.dat ---
        buf.clear();
        buf.reserve(bonds.size() * 176 + 256);
        buf.append("id,particle_a_id,particle_b_id,rest_length,extension,fx,fy,fz,active,heat_flow_a,damage,maximum_tensile_strain,elastic_energy,fracture_energy,dissipated_fracture_energy,source,maximum_compressive_strain,strength_scale\n");
        for (const auto& b : bonds)
        {
            AppendNumber(buf, b.id);                                    buf.push_back(',');
            AppendNumber(buf, particles[b.idx_a].id);                   buf.push_back(',');
            AppendNumber(buf, particles[b.idx_b].id);                   buf.push_back(',');
            AppendNumber(buf, b.rest_length);                           buf.push_back(',');
            AppendNumber(buf, b.extension);                             buf.push_back(',');
            AppendNumber(buf, b.force_a.x()); buf.push_back(',');
            AppendNumber(buf, b.force_a.y()); buf.push_back(',');
            AppendNumber(buf, b.force_a.z());                           buf.push_back(',');
            AppendNumber(buf, b.active ? 1 : 0);                        buf.push_back(',');
            AppendNumber(buf, b.heat_flow_a);                           buf.push_back(',');
            AppendNumber(buf, b.damage);                                buf.push_back(',');
            AppendNumber(buf, b.maximum_tensile_strain);                buf.push_back(',');
            AppendNumber(buf, b.elastic_energy);                        buf.push_back(',');
            AppendNumber(buf, b.fracture_energy);                       buf.push_back(',');
            AppendNumber(buf, b.dissipated_fracture_energy);            buf.push_back(',');
            AppendNumber(buf, b.source);                                buf.push_back(',');
            AppendNumber(buf, b.maximum_compressive_strain);            buf.push_back(',');
            AppendNumber(buf, b.strength_scale);                        buf.push_back('\n');
        }
        WriteTextFile(back_folder + "/bonds.dat", buf);

        Log::info("DEMSolver::BackupField: {} particles and {} bonds written to {}",
            m_dem_data->GetParticleNum(), m_dem_data->GetBondNum(), back_folder);

        if (m_voronoi_valid)
        {
            const auto& faces = m_gas_voronoi_faces;
            buf.clear();
            buf.reserve(faces.size() * 104 + 256);
            buf.append("particle_a_id,particle_b_id,length,area,phase_a,phase_b,pressure_a,pressure_b\n");
            for (const auto& face : faces)
            {
                const auto& pa = particles[face.idx_a];
                const auto& pb = particles[face.idx_b];
                AppendNumber(buf, pa.id);          buf.push_back(',');
                AppendNumber(buf, pb.id);          buf.push_back(',');
                AppendNumber(buf, face.length);    buf.push_back(',');
                AppendNumber(buf, face.area);      buf.push_back(',');
                AppendNumber(buf, pa.phase);       buf.push_back(',');
                AppendNumber(buf, pb.phase);       buf.push_back(',');
                AppendNumber(buf, pa.gas_pressure); buf.push_back(',');
                AppendNumber(buf, pb.gas_pressure); buf.push_back('\n');
            }
            WriteTextFile(back_folder + "/gas_voronoi_faces.dat", buf);
        }
    }

} // namespace zaran
