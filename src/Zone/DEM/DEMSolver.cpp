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
        // 切向接触阻尼（c_t = λ·c_n）。λ ≤ 0 时接触模型内部完全跳过该分支。
        m_contact_model->SetTangentialDampingScale(para->GetTangentialDampingScale());
        if (para->GetTangentialDampingScale() > 0.0)
            Log::info("DEMSolver: 切向接触阻尼已启用，c_t = {:g}·c_n",
                para->GetTangentialDampingScale());
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
        // --- 多通道失效：切向独立断裂能阈值 UtIII = α·UnIII ---
        // 必须在 fracture_energy 全部确定之后（含文件给定与 surface_energy 派生两种来源）。
        // α ≤ 0（默认）时全部为 0 ⇒ 切向断裂/压剪两条通道完全不参与计算。
        {
            const double shear_ratio = para->GetBondShearEnergyRatio();
            index_type shear_bonds = 0;
            double shear_ratio_min = 1.0e30, shear_ratio_max = 0.0;
            for (auto& bond : m_dem_data->GetBonds())
            {
                bond.fracture_energy_shear = (shear_ratio > 0.0)
                    ? shear_ratio * bond.fracture_energy : 0.0;
                if (bond.fracture_energy_shear > 0.0 && bond.fracture_energy > 0.0)
                {
                    ++shear_bonds;
                    const double r = bond.fracture_energy_shear / bond.fracture_energy;
                    shear_ratio_min = std::min(shear_ratio_min, r);
                    shear_ratio_max = std::max(shear_ratio_max, r);
                }
            }
            if (shear_bonds > 0)
            {
                Log::info("DEM 多通道失效: {} 条键启用切向独立断裂通道，"
                    "UtIII/UnIII ∈ [{:.4f}, {:.4f}]（切向先断时按 LSM GBratio 折减法向阈值）",
                    shear_bonds, shear_ratio_min, shear_ratio_max);
            }
            else if (shear_ratio > 0.0)
            {
                Log::warn("DEM 多通道失效: dem.bond_shear_energy_ratio = {:g} > 0，"
                    "但没有任何键有正的断裂能（fracture_energy ≤ 0 ⇒ 需要 dem.surface_energy），"
                    "切向独立断裂通道实际未生效。", shear_ratio);
            }
        }
        // γ 控制的内聚律参数（δ_p、δ_f）必须在 fracture_energy 与 strength_scale
        // 都确定之后、且覆盖"文件给定 + 自动建链"两种来源。
        PrepareCohesiveLaw();
        // 每个键的初始半径（径向损伤剖面的拉格朗日坐标），建键后赋值一次
        {
            const Eigen::Vector3d center = para->GetGasCavityCenter();
            for (auto& b : m_dem_data->GetBonds())
            {
                if (b.idx_a >= m_dem_data->GetParticleNum()
                    || b.idx_b >= m_dem_data->GetParticleNum())
                {
                    continue;
                }
                const auto& pa = m_dem_data->GetParticles()[b.idx_a];
                const auto& pb = m_dem_data->GetParticles()[b.idx_b];
                b.initial_radius = std::hypot(
                    0.5 * (pa.pos.x() + pb.pos.x()) - center.x(),
                    0.5 * (pa.pos.y() + pb.pos.y()) - center.y());
            }
        }
        UpdateGasVoronoiMesh(true);

        // --- 可选：高压气腔加载（孔洞内充压气体把脆性材料撑碎）---
        // 必须放在键合参数（含 m_lattice_spacing）确定之后：边界环识别用格距作容差。
        if (para->GetGasCavityEnabled())
        {
            CavityGasOptions gas;
            gas.enabled = true;
            gas.center = para->GetGasCavityCenter();
            gas.cavity_radius = para->GetGasCavityRadius();
            gas.pressure_initial = para->GetGasPressureInitial();
            gas.polytropic_index = para->GetGasPolytropicIndex();
            gas.ramp_time = para->GetGasRampTime();
            gas.thickness = para->GetLatticeThickness();
            gas.shell_tolerance = para->GetGasCavityShellTolerance();
            gas.disk_radius = para->GetSampleRadius();
            gas.vent_area_ratio = para->GetGasCavityVentAreaRatio();
            gas.pressure_cap = para->GetGasCavityPressureCap();
            m_cavity_gas.Initialize(m_dem_data->GetParticles(), gas, m_lattice_spacing);
        }

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

        // 高压气腔加载：受力必须在积分**之前**（用步初位置确定边界环几何），
        // 体积/内能的更新必须在积分**之后**；两者共用同一套边界离散 ⇒ W ≡ pΔV。
        const double dt = GetDEMParam()->GetTimeStep();
        const double time = static_cast<double>(m_dem_step) * dt;
        CalcCavityGasForce(time);

        Integrate();
        UpdateCavityGasState(dt, time);

        // 机械耗散累计：m_tmp_dissipation 在 ZeroForce 里清零、本步值即 J，
        // 单独累计一份供能量账本闭合检查（与写入温度的用法互不干扰）。
        for (index_type i = 0; i < m_tmp_dissipation.size(); ++i)
        {
            m_dissipation_total += m_tmp_dissipation[i];
        }

        ++m_dem_step;
    }

    void DEMSolver::CalcCavityGasForce(double time)    {
        if (!m_cavity_gas.IsEnabled()) return;
        m_cavity_gas.ApplyPressureForce(m_dem_data->GetParticles(), time);
    }

    void DEMSolver::UpdateCavityGasState(double dt, double time)
    {
        if (!m_cavity_gas.IsEnabled()) return;
        m_cavity_gas.UpdateState(m_dem_data->GetParticles(), dt, time);
    }

    DEMSolver::EnergyBudget DEMSolver::GetEnergyBudget() const
    {
        EnergyBudget budget;
        const auto& particles = m_dem_data->GetParticles();
        for (const auto& p : particles)
        {
            budget.kinetic_translation += 0.5 * p.mass * p.vel.squaredNorm();
            budget.kinetic_rotation += 0.5 * p.inertia * p.mass * p.radius * p.radius
                * p.omega.squaredNorm();
            budget.max_speed = std::max(budget.max_speed, p.vel.norm());
            const double dx = p.pos.x() - GetDEMParam()->GetGasCavityCenter().x();
            const double dy = p.pos.y() - GetDEMParam()->GetGasCavityCenter().y();
            budget.max_fragment_distance = std::max(budget.max_fragment_distance,
                std::sqrt(dx * dx + dy * dy));
        }
        for (const auto& b : m_dem_data->GetBonds())
        {
            budget.bond_elastic += b.elastic_energy;
            budget.bond_fracture += b.dissipated_fracture_energy;
            // 切向独立断裂的耗散单列统计，并并入断裂耗散总量
            // （未启用切向通道时恒为 0，历史账本逐位不变）。
            if (b.dissipated_shear_energy > 0.0)
            {
                budget.bond_fracture += b.dissipated_shear_energy;
                budget.bond_shear_fracture += b.dissipated_shear_energy;
            }
            if (b.active) ++budget.active_bonds; else ++budget.broken_bonds;
        }
        budget.dissipation = m_dissipation_total;
        budget.contact_elastic = m_contact_elastic;
        budget.rebound_loss = m_rebound_loss;
        return budget;
    }

    // ==================================================================
    // γ 控制的双线性内聚律：逐键参数（δ_p、δ_f）与自洽性检查
    // ==================================================================
    // 关键关系（推导见 DEMBond.h 的注释）：
    //   曲线下面积 = ½·F_p·δ_f，令其等于 Efract = γ·L0·t/√3·(1−Egb/γ)·s²
    //   ⇒ δ_f = 2·Efract/(k_n·δ_p)        （solid_area_scale 自动约掉）
    // 物理含义：**δ_f ≥ δ_p 等价于 L0 ≤ γ_eff·E/σ_t²** —— γ 与 σ_t 是两个独立材料量，
    // 网格越粗就越无法同时表示它们；这是"裂纹能被分辨开"的先决条件，必须显式检查。
    void DEMSolver::PrepareCohesiveLaw()
    {
        auto* para = GetDEMParam();
        const bool want = para->GetBondCohesiveEnabled();
        const double eps_peak_input = para->GetBondPeakStrain();
        const bool have_peak = eps_peak_input > 0.0 && eps_peak_input < 1.0e29;

        index_type cohesive_bonds = 0;
        index_type degenerate_bonds = 0;
        double min_ratio = 1.0e30;
        double max_ratio = 0.0;
        double min_zone_bonds = 1.0e30;
        const double thickness = para->GetLatticeThickness();
        const auto& particles = m_dem_data->GetParticles();

        for (auto& b : m_dem_data->GetBonds())
        {
            b.cohesive = false;
            b.cohesive_peak_separation = 0.0;
            b.cohesive_failure_separation = 0.0;
            if (!want || !have_peak || b.fracture_energy <= 0.0
                || b.normal_stiffness <= 0.0 || b.rest_length <= 0.0)
            {
                continue;
            }
            ++cohesive_bonds;
            const DEMParticle& bp = particles[static_cast<size_t>(b.idx_a)];
            const double s = b.strength_scale;
            // δ_p 用**未放大**的峰值应变乘 L0（含逐键强度折减 s）
            const double dp = eps_peak_input * s * b.rest_length;
            // Efract 已按 s² 缩放（强度按 s 缩放时应力应变同乘 s）
            const double efract_scaled = b.fracture_energy * s * s;
            const double df = 2.0 * efract_scaled / (b.normal_stiffness * dp);
            b.cohesive_peak_separation = dp;
            b.cohesive_failure_separation = df;

            if (!(df > dp))
            {
                // 储能还没到峰值就已经用完了全部断裂能 ⇒ 无法表示软化段。
                // 默认直接报错（静默退化是最难发现的一类错误），可用
                // dem.bond_cohesive_strict = false 降级为"脆断"。
                b.cohesive = false;
                ++degenerate_bonds;
                continue;
            }
            b.cohesive = true;
            const double ratio = df / dp;
            min_ratio = std::min(min_ratio, ratio);
            max_ratio = std::max(max_ratio, ratio);

            // 内聚区长度（平面应力）：l_cz = E·G_f/σ_p²，G_f = Efract/A_bond，σ_p = F_p/A_bond。
            // 注意 **δ_f 本身不是内聚区长度** —— 它是材料分离度（~10 μm），而内聚区
            // 是裂尖附近牵引力起作用的那一段（~mm 量级），必须按 E·G_f/σ_t² 估。
            const double a_bond = (thickness > 0.0)
                ? b.rest_length * thickness / 1.7320508075688772 : 0.0;
            if (a_bond > 0.0 && bp.young_modulus > 0.0)
            {
                const double sigma_p = b.normal_stiffness * dp / a_bond;
                if (sigma_p > 0.0)
                {
                    const double zone = bp.young_modulus * (efract_scaled / a_bond)
                        / (sigma_p * sigma_p);
                    min_zone_bonds = std::min(min_zone_bonds, zone / b.rest_length);
                }
            }
        }

        if (cohesive_bonds == 0)
        {
            return;
        }

        // 网格分辨力的判据：内聚区至少要覆盖 min_bonds 个键长
        const double min_bonds_required = para->GetBondCohesiveMinSofteningBonds();
        if (min_bonds_required > 0.0
            && min_zone_bonds < min_bonds_required)
        {
            Log::warn("DEM cohesive law: 内聚区只有 {:.2f} 个键长（要求 ≥ {:g}）—— "
                "牵引力起作用的区域落在网格尺度以内，裂纹无法沿路径扩展而会弥散粉碎。"
                "增大断裂功（dem.surface_energy）或加密网格：需要 "
                "l_cz = E·γ/σ_t² ≥ {:g}·L0。", min_zone_bonds, min_bonds_required,
                min_bonds_required);
        }

        if (degenerate_bonds > 0)
        {
            const std::string msg = "DEM cohesive law: " + std::to_string(degenerate_bonds)
                + " / " + std::to_string(cohesive_bonds)
                + " 条键的 δ_f ≤ δ_p（峰值储能已超过断裂能）⇒ 无法表示软化段。"
                "请提高 dem.surface_energy 或降低 dem.bond_peak_strain；"
                "自洽条件为 γ_eff ≥ E·ε_p²·L0（等价 L0 ≤ γ_eff·E/σ_t²）。";
            if (para->GetBondCohesiveStrict())
            {
                throw ZaranError("DEM cohesive law 参数不自洽: " + msg);
            }
            Log::warn(msg + " 已按 dem.bond_cohesive_strict = false 退化为脆断口径。");
        }

        Log::info("DEM cohesive law: {} 条键启用 γ 控制双线性内聚律；"
            "δ_f/δ_p ∈ [{:.3f}, {:.3f}]，内聚区 ≥ {:.2f} 个键长；"
            "曲线下面积 ≡ fracture_energy（γ·A）",
            cohesive_bonds, min_ratio, max_ratio, min_zone_bonds);
    }

    // ==================================================================
    // 破碎形态统计：未断键连通分量 = 完整碎块；已断键连通 = 裂纹带
    // ==================================================================
    // 只看"断键率"会被骗：94% 断键既可能是 3 条主裂纹（含大量被切开的面），
    // 也可能是全盘粉化。碎块尺寸分布才是"几大块 + 飞溅细粒"的定量口径。
    DEMSolver::FragmentStats DEMSolver::ComputeFragmentStats() const
    {
        FragmentStats stats;
        const auto& particles = m_dem_data->GetParticles();
        const auto& bonds = m_dem_data->GetBonds();
        const index_type n = particles.size();
        const index_type e = bonds.size();
        if (n == 0)
        {
            return stats;
        }

        // 并查集（路径压缩 + 按大小合并）
        auto make_uf = [&](std::vector<index_type>& parent)
        {
            parent.resize(static_cast<size_t>(n));
            for (index_type i = 0; i < n; ++i) parent[static_cast<size_t>(i)] = i;
        };
        auto find = [&](std::vector<index_type>& parent, index_type a)
        {
            while (parent[static_cast<size_t>(a)] != a)
            {
                parent[static_cast<size_t>(a)] = parent[static_cast<size_t>(parent[static_cast<size_t>(a)])];
                a = parent[static_cast<size_t>(a)];
            }
            return a;
        };
        auto unite = [&](std::vector<index_type>& parent, index_type a, index_type b)
        {
            const index_type ra = find(parent, a);
            const index_type rb = find(parent, b);
            if (ra != rb) parent[static_cast<size_t>(ra)] = rb;
        };

        // ---- 完整碎块（未断键连通分量）----
        std::vector<index_type> pu;
        make_uf(pu);
        for (const auto& b : bonds)
        {
            if (b.active && b.idx_a < n && b.idx_b < n) unite(pu, b.idx_a, b.idx_b);
        }
        std::vector<index_type> count(static_cast<size_t>(n), 0);
        std::vector<double> mass(static_cast<size_t>(n), 0.0);
        double total_mass = 0.0;
        for (index_type i = 0; i < n; ++i)
        {
            const index_type r = find(pu, i);
            ++count[static_cast<size_t>(r)];
            mass[static_cast<size_t>(r)] += particles[static_cast<size_t>(i)].mass;
            total_mass += particles[static_cast<size_t>(i)].mass;
        }
        for (index_type i = 0; i < n; ++i)
        {
            if (count[static_cast<size_t>(i)] == 0) continue;
            ++stats.fragments;
            if (count[static_cast<size_t>(i)] == 1) ++stats.isolated_particles;
            const double frac = (total_mass > 0.0)
                ? mass[static_cast<size_t>(i)] / total_mass : 0.0;
            if (count[static_cast<size_t>(i)] > stats.largest_fragment)
            {
                stats.largest_fragment = count[static_cast<size_t>(i)];
                stats.largest_fragment_mass_fraction = frac;
            }
            if (count[static_cast<size_t>(i)] >= 20)
            {
                ++stats.fragments_ge_20;
                stats.fragments_ge_20_mass_fraction += frac;
            }
        }

        // ---- 裂纹带（已断键连通分量）----
        std::vector<index_type> pb;
        make_uf(pb);
        std::vector<char> touched(static_cast<size_t>(n), 0);
        for (const auto& b : bonds)
        {
            if (b.active || b.idx_a >= n || b.idx_b >= n) continue;
            unite(pb, b.idx_a, b.idx_b);
            touched[static_cast<size_t>(b.idx_a)] = 1;
            touched[static_cast<size_t>(b.idx_b)] = 1;
        }
        std::vector<index_type> bc(static_cast<size_t>(n), 0);
        for (index_type i = 0; i < n; ++i)
        {
            if (!touched[static_cast<size_t>(i)]) continue;
            ++bc[static_cast<size_t>(find(pb, i))];
        }
        for (index_type i = 0; i < n; ++i)
        {
            stats.largest_crack_band = std::max(stats.largest_crack_band,
                bc[static_cast<size_t>(i)]);
        }

        // ---- 径向断键率（12 个等宽环带）----
        // ⚠ 归一化必须用**固定的参考半径**（试件初始半径），不能用"当前最大半径"：
        //   一旦碎片飞出，max_radius 会随时间暴涨，环带会跟着漂移，
        //   而飞出去的高速碎块（键全断）会全部落进最外两环、把剖面刷成 100%，
        //   读起来像是"外缘被碎成粉末"，与事实正好相反。同时**只统计参考半径以内**的键，
        //   让这条曲线始终是"原始圆盘上的损伤分布"。
        const Eigen::Vector3d center = GetDEMParam()->GetGasCavityCenter();
        for (const auto& p : particles)
        {
            stats.max_radius = std::max(stats.max_radius,
                std::hypot(p.pos.x() - center.x(), p.pos.y() - center.y()));
        }
        const double ref_radius = (GetDEMParam()->GetSampleRadius() > 0.0)
            ? GetDEMParam()->GetSampleRadius() : stats.max_radius;
        if (ref_radius > 0.0)
        {
            double band_total[12] = { 0.0 };
            double band_broken[12] = { 0.0 };
            for (index_type k = 0; k < e; ++k)
            {
                const auto& b = bonds[static_cast<size_t>(k)];
                if (b.idx_a >= n || b.idx_b >= n) continue;
                // 用**初始**半径（拉格朗日坐标）：空腔膨胀/碎块飞出后，
                // 按当前位置统计会让环带漂移，读出来的图像完全反过来。
                double r;
                if (b.initial_radius >= 0.0)
                {
                    r = b.initial_radius;
                }
                else
                {
                    const double x = 0.5 * (particles[static_cast<size_t>(b.idx_a)].pos.x()
                        + particles[static_cast<size_t>(b.idx_b)].pos.x()) - center.x();
                    const double y = 0.5 * (particles[static_cast<size_t>(b.idx_a)].pos.y()
                        + particles[static_cast<size_t>(b.idx_b)].pos.y()) - center.y();
                    r = std::hypot(x, y);
                }
                if (r > ref_radius) continue;
                int ib = static_cast<int>(r / ref_radius * 12.0);
                if (ib < 0) ib = 0;
                if (ib > 11) ib = 11;
                band_total[ib] += 1.0;
                if (!b.active) band_broken[ib] += 1.0;
            }
            for (int i = 0; i < 12; ++i)
            {
                stats.radial_broken_fraction[i] = (band_total[i] > 0.0)
                    ? band_broken[i] / band_total[i] : 0.0;
            }
        }
        return stats;
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
            // 切向弹簧已独立断裂时不再累积切向位移：否则 δ_t 会重新增长，
            // 让后面所有"切向储能"的口径（判据与记账）都虚高。
            if (bond.shear_broken) bond.delta_t.setZero();
            bond.extension = length - bond.rest_length;

            // 弹性能。内聚模式下受力刚度已按 (1−d) 折减，储能必须同步折减，
            // 否则账本里的"可逆储能"与卸载路径不一致（这是账本闭合的必要条件）。
            bond.elastic_energy = 0.5 * solid_area_scale * (bond.cohesive ? (1.0 - bond.damage) : 1.0)
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
            // --- 压剪通道（LSM breakmod = 2）---
            // 压缩侧（rn ≤ r0 ⇒ extension < 0）且切向有明显错动、切向弹簧储能达到 UtIII
            // ⇒ 剪切型破坏。它与"压缩压溃"（breakmod = 3，纯法向压缩过深）是两条独立通道：
            // 前者由剪切滑移驱动（孔洞/裂隙壁面的错动正是这种），后者由法向压密驱动。
            // 只有 α > 0（显式启用切向通道）时才存在；默认关闭 ⇒ 逐位回退。
            if (bond.cohesive && bond.extension < 0.0 && bond.fracture_energy_shear > 0.0
                && !bond.shear_broken && bond.tangential_stiffness > 0.0
                && bond.delta_t.squaredNorm() > 0.0)
            {
                const double u_t_eff = bond.fracture_energy_shear * solid_area_scale
                    * strength_scale * strength_scale;
                const double Et_now = 0.5 * solid_area_scale * bond.tangential_stiffness
                    * bond.delta_t.squaredNorm();
                if (u_t_eff > 0.0 && Et_now >= u_t_eff)
                {
                    bond.break_mode = 2;
                    // 法向压缩储能与切向弹簧储能分别记账，避免与 dissipated_shear_energy 重复
                    const double E_n_elastic = 0.5 * solid_area_scale
                        * (1.0 - bond.damage) * bond.normal_stiffness
                        * bond.extension * bond.extension;
                    bond.dissipated_fracture_energy += E_n_elastic;
                    bond.dissipated_shear_energy += Et_now;
                    bond.elastic_energy = 0.0;
                    bond.damage = 1.0;
                    bond.active = false;
                    bond.force_a.setZero();
                    continue;
                }
            }
            if (compression_limit > 0.0 && compressive_strain > compression_limit)
            {
                bond.break_mode = 3;
                bond.dissipated_fracture_energy += bond.elastic_energy;
                // 断开的键不再储能：清零后 Σelastic_energy 只统计完好键，
                // 不会与 dissipated_fracture_energy 重复计数（否则能量账会算错）。
                bond.elastic_energy = 0.0;
                bond.damage = 1.0;
                bond.active = false;
                bond.force_a.setZero();
                continue;
            }

            if (bond.cohesive && bond.extension >= 0.0)
            {
                // ============ γ 控制的双线性内聚律 ============
                // 牵引-分离曲线：δ ≤ δ_p 线性升至 F_p = k_n·δ_p；δ_p<δ<δ_f 线性软化到 0；
                // δ ≥ δ_f 失效。构造保证曲线下**总面积 = fracture_energy·s²**（见
                // DEMBond::cohesive_failure_separation 的推导），所以"造新表面要花多少能量"
                // 是一个真正的材料输入，而不是由网格/阈值派生的量。
                //
                // 损伤由 δ 的解析式给出（形如 LSM 的二阶段软化，但面积严格守恒）：
                //   d(δ) = 1 − (δ_p/δ)·(δ_f−δ)/(δ_f−δ_p)
                // 于是 F = (1−d)·k_n·δ 恰好落在从 (δ_p,F_p) 到 (δ_f,0) 的直线上。
                // 可逆弹性能取卸载回原点的 ½(1−d)k_nδ²，耗散能取包络面积减可逆部分，
                // 二者之和逐位等于键力对外做的功 —— 这是能量账本能闭合的前提。
                const double dp = bond.cohesive_peak_separation;
                const double e = bond.extension;
                const double k_eff = solid_area_scale * bond.normal_stiffness;
                const double f_peak = k_eff * dp;
                // 断裂能的逐键缩放（强度按 s 缩放时应力应变同乘 s ⇒ 能量密度 ∝ s²）
                const double energy_scale = solid_area_scale * strength_scale * strength_scale;

                // ---- 切向独立断裂通道（LSM breakmod = 4）----
                // 切向弹簧储能 Et 达到 UtIII 时，切向先于法向失效：
                //   ① 切向不再承载（bond.shear_broken ⇒ 下面 force_t 取 0）；
                //   ② 按 LSM 的 GBratio 折减规则压低**法向能量预算**
                //      （Un_eff ← Un_eff − Ut_eff），于是"造切向新表面"与
                //      "造法向新表面"共享同一份材料断裂能。
                // ⚠ 折减的是能量阈值，**不是 δ_f** —— LSM 里 rnmax 是材料常数、
                // 全程不改动。若按面积重解 δ_f，会在"法向已经吸收了部分包络功之后
                // 才切向断裂"时把那段已吸收的功抹掉（实测丢 1.3e-5 相对能量）。
                // 只在 bond.fracture_energy_shear > 0（α > 0）时进入，默认关闭。
                const bool shear_independent = bond.fracture_energy_shear > 0.0;
                if (shear_independent && !bond.shear_broken
                    && bond.tangential_stiffness > 0.0)
                {
                    const double u_t_eff = bond.fracture_energy_shear * energy_scale;
                    const double Et = 0.5 * solid_area_scale * bond.tangential_stiffness
                        * bond.delta_t.squaredNorm();
                    if (u_t_eff > 0.0 && Et >= u_t_eff)
                    {
                        bond.shear_broken = true;
                        bond.break_mode = 4;
                        // 切向弹簧的储能全部不可逆地转为热（刚度归零 ⇒ 不再可逆）
                        bond.dissipated_shear_energy += Et;
                        bond.delta_t.setZero();
                        // LSM 的 GBratio 折减：只降**能量阈值**，δ_f（材料常数）不动：
                        //   GBratio_new = (UnIII − UtIII)·GBratio_old / UnIII
                        // 于是"造切向新表面"与"造法向新表面"共享同一份材料断裂能。
                        bond.cohesive_budget_scale *=
                            (bond.fracture_energy - bond.fracture_energy_shear)
                            / bond.fracture_energy;
                    }
                }

                const double df = bond.cohesive_failure_separation;
                // 损伤由**历史最大**拉伸分离度驱动 ⇒ 卸载不恢复损伤；卸载刚度取 (1−d)k_n，
                // 即线性卸载回原点（标准损伤力学口径）。
                const double e_drive = std::max(e,
                    bond.maximum_tensile_strain * bond.rest_length);

                double d_local;
                if (e_drive <= dp)
                {
                    d_local = 0.0;
                }
                else if (e_drive < df)
                {
                    d_local = 1.0 - (dp / e_drive) * ((df - e_drive) / (df - dp));
                }
                else
                {
                    d_local = 1.0;
                }
                if (d_local < 0.0) d_local = 0.0;
                if (d_local > 1.0) d_local = 1.0;
                bond.damage = d_local;

                // 包络（牵引-分离曲线）在 e_drive 处的累积面积 = 该键吸收过的总功
                double area_env;
                if (e_drive <= dp)
                {
                    area_env = 0.5 * k_eff * e_drive * e_drive;
                }
                else
                {
                    const double f_at = (df > dp)
                        ? f_peak * (df - e_drive) / (df - dp) : 0.0;
                    area_env = 0.5 * f_peak * dp + 0.5 * (f_peak + f_at) * (e_drive - dp);
                    if (area_env < 0.0) area_env = 0.5 * f_peak * dp;
                }

                // 可逆部分用**当前**分离度，耗散部分 = 吸收的总功 − 包络上那一点的可逆能。
                // 于是 (可逆 + 耗散) = 净吸收能量，卸载时会把能量还回去，账本自动闭合。
                const double reversible_normal = 0.5 * k_eff * (1.0 - d_local) * e * e;
                const double reversible_at_drive = 0.5 * k_eff * (1.0 - d_local)
                    * e_drive * e_drive;
                // 切向退化因子：默认与法向损伤绑定（历史口径，逐位不变）；
                // 启用切向独立通道后，切向弹簧有自己的强度极限、不随法向损伤折减
                // （与 LSM 的 En / Et 两条独立通道一致），切向断裂后直接归零。
                const double tangential_factor = shear_independent
                    ? (bond.shear_broken ? 0.0 : 1.0) : (1.0 - d_local);
                bond.elastic_energy = reversible_normal
                    + 0.5 * solid_area_scale * tangential_factor
                        * bond.tangential_stiffness * bond.delta_t.squaredNorm();
                double dissipated = area_env - reversible_at_drive;
                if (dissipated < 0.0) dissipated = 0.0;
                bond.dissipated_fracture_energy = dissipated;

                // 失效判据。默认（未启用切向通道）逐位等价于历史写法 `e_drive >= df`：
                // area_env 在 e_drive ≥ δ_f 处恰等于 Un_eff（曲线下总面积）。
                bool bond_failed = (e_drive >= df);
                // 启用切向通道时的两条额外判据（LSM）
                double envelope_total = 0.0, budget = 0.0;
                if (shear_independent)
                {
                    envelope_total = bond.fracture_energy * energy_scale;
                    budget = envelope_total * bond.cohesive_budget_scale;
                    const double Et_now = bond.shear_broken ? 0.0
                        : 0.5 * solid_area_scale * bond.tangential_stiffness
                            * bond.delta_t.squaredNorm();
                    // 拉伸混合通道（LSM breakmod = 1）：En + Et ≥ Un_eff。
                    // 剪切储能参与把键拉断 —— 这是"法向还没到软化终点、但剪切已经
                    // 帮着把键带坏"的情形，也让裂纹走向对剪切场敏感。
                    if (!bond_failed && budget > 0.0 && area_env + Et_now >= budget)
                    {
                        bond_failed = true;
                    }
                }
                if (bond_failed)
                {
                    // 完全失效。**默认路径保持历史表达式逐位不变**（闭式 Un_eff，
                    // 避免最后一步数值积分漂移）；启用切向通道时按"实际吸收"记账：
                    //   · 键断开后，当前储存在键里的可逆能全部转为不可逆耗散；
                    //   · 法向部分 = 已走完的包络面积，但不超过（已折减的）能量预算；
                    //   · 切向部分 = 切向弹簧储能（切向已断则此前已记账，此处为 0）。
                    // 于是 Σ(法向 + 切向) 恰好等于该键从系统吸收的总能量，账本逐位闭合。
                    bond.break_mode = 1;
                    bond.damage = 1.0;
                    if (shear_independent)
                    {
                        const double Et_now = bond.shear_broken ? 0.0
                            : 0.5 * solid_area_scale * bond.tangential_stiffness
                                * bond.delta_t.squaredNorm();
                        bond.dissipated_shear_energy += Et_now;
                        // e_drive ≥ δ_f 时包络已经走完 ⇒ 吸收量取整条曲线的面积
                        const double absorbed = (e_drive >= df) ? envelope_total : area_env;
                        bond.dissipated_fracture_energy =
                            (absorbed < budget) ? absorbed : budget;
                    }
                    else
                    {
                        bond.dissipated_fracture_energy = bond.fracture_energy
                            * solid_area_scale * strength_scale * strength_scale;
                    }
                    bond.elastic_energy = 0.0;
                    bond.active = false;
                    bond.force_a.setZero();
                    continue;
                }
            }
            else if (bond.fracture_energy > 0.0 && bond.extension >= 0.0
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
            if (!bond.cohesive && bond.fracture_energy <= 0.0)
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
            // 切向承载因子。默认与历史逐位一致：(1−d)·sas 作用在刚度上、阻尼不受损伤影响。
            // 启用切向独立通道后，切向不随法向损伤折减（独立强度通道），切向断裂后完全退出。
            const bool shear_carry_independent = bond.fracture_energy_shear > 0.0;
            const double tangential_stiffness_factor = shear_carry_independent
                ? (bond.shear_broken ? 0.0 : solid_area_scale)
                : (1.0 - bond.damage) * solid_area_scale;
            const double tangential_damping_factor =
                (shear_carry_independent && bond.shear_broken) ? 0.0 : 1.0;
            const Eigen::Vector3d force_t =
                -tangential_stiffness_factor * bond.tangential_stiffness * bond.delta_t
                - tangential_damping_factor * bond.tangential_damping * velocity_t;
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
        // 本步接触弹簧储能之和，同样每步清零后由 CalcContactForce 累加
        m_contact_elastic = 0.0;
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

        // 这是**速度改写**而不是力的做功：把它移出的动能单独记账，
        // 否则剧烈压缩段的能量账本会凭空缺口（实测可达十几个百分点）。
        const double ke_before = 0.5 * pa.mass * pa.vel.squaredNorm()
            + 0.5 * pb.mass * pb.vel.squaredNorm();

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
        const double ke_after = 0.5 * pa.mass * pa.vel.squaredNorm()
            + 0.5 * pb.mass * pb.vel.squaredNorm();
        m_rebound_loss += ke_before - ke_after; // >0 表示回弹移出了动能
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

    // ==================================================================
    // 碎后接触：滚动阻力矩
    // ==================================================================
    // 为什么需要它：DEM 里的"粒子"是光滑球，一旦碎块之间只剩球-球接触，
    // 它们可以自由滚动 ⇒ 碎块堆会像干砂一样摊平，而不是像真实的不规则
    // 碎块那样互相"角锁"、堆成有休止角的堆。滚动阻力矩就是角锁的最低阶
    // 等效描述（Ai, Mistry, Ibrahim 等 2011 的简化模型）：
    //     M_r = −μ_r · R* · |F_n| · ω̂_rel
    // 力矩与相对角速度反平行 ⇒ 耗散功率 |M_r|·|ω_rel| 恒 ≥ 0。
    // μ_r 有**长度量纲**（不是无量纲摩擦系数），量级 ~ 粒子半径的 0.01–0.1 倍。
    void DEMSolver::ApplyRollingResistance(const DEMParticle& pa, const DEMParticle& pb,
                                           DEMContact& contact, double dt) const
    {
        contact.torque_r.setZero();
        const double mu_r = GetDEMParam()->GetRollingFriction();
        if (!(mu_r > 0.0) || !(dt > 0.0)) return;

        const double Fn_mag = contact.force_n.norm();
        if (!(Fn_mag > 0.0)) return;

        // 相对角速度（A 相对 B）。注意接触的切向速度里用的是
        // (r_a·ω_a + r_b·ω_b)×n，方向约定与这里的 ω_a − ω_b 不同；
        // 滚动阻力只看两端角速度的**差**，与接触点位置无关。
        const Eigen::Vector3d omega_rel = pa.omega - pb.omega;
        const double w = omega_rel.norm();
        if (!(w > 1.0e-30)) return;

        // 等效半径与等效转动惯量（与法向等效质量同一套"不可动端按无限大"口径）
        const double denom = pa.radius + pb.radius;
        if (!(denom > 0.0)) return;
        const double R_star = (pa.radius * pb.radius) / denom;

        const double I_a = pa.inertia * pa.mass * pa.radius * pa.radius;
        const double I_b = pb.inertia * pb.mass * pb.radius * pb.radius;
        double I_eff = 0.0;
        if (!pa.IsDynamic())      I_eff = I_b;
        else if (!pb.IsDynamic()) I_eff = I_a;
        else                      I_eff = (I_a * I_b) / (I_a + I_b);

        double M = mu_r * R_star * Fn_mag;
        // 过冲截断：一步内最多把相对角速度恰好减到 0。没有这一步，力矩会在
        // 静止接触上把 ω_rel 推过零并反向，产生高频抖动（永远收敛不到"不滚"）。
        if (I_eff > 0.0)
        {
            const double M_limit = I_eff * w / dt;
            if (M > M_limit) M = M_limit;
        }
        contact.torque_r = -M * (omega_rel / w);
        contact.dissipation += M * w * dt;
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
            // 滚动阻力矩：必须在切向之后（它依赖 contact.force_n，且与切向力无关）
            ApplyRollingResistance(pa, pb, c, dt);
            SaveContactHistory(c);

            // 机械耗散（法向阻尼 + 切向摩擦 + 切向阻尼 + 滚动阻力）对半分给两端 → 后续计入温度
            const double dissipated = 0.5 * c.dissipation;
            m_tmp_dissipation[c.idx_a] += dissipated;
            m_tmp_dissipation[c.idx_b] += dissipated;

            // 接触弹簧的可逆储能（½k_nδ² + ½k_t|δt|²）：剧烈压缩下可达焦耳量级，
            // 必须在能量账本里单列，否则压碎段会"凭空少掉"十几个百分点。
            m_contact_elastic += c.elastic_energy;

            Eigen::Vector3d F = c.force_n + c.force_t;

            // 作用-反作用
            pa.force += F;
            pa.torque += c.contact_point.cross(c.force_t) - pa.pos.cross(c.force_t);
            pb.force -= F;
            pb.torque -= c.contact_point.cross(c.force_t) - pb.pos.cross(c.force_t);
            // 滚动阻力是**力偶**：两端力矩大小相等、方向相反（不产生净力矩）
            pa.torque += c.torque_r;
            pb.torque -= c.torque_r;
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
                ApplyRollingResistance(pa, pb_wall, c, dt);
                SaveContactHistory(c);

                // 墙面不可动（无限质量），该对的耗散能量全部计入粒子
                m_tmp_dissipation[pi] += c.dissipation;

                pa.force += c.force_n + c.force_t;
                pa.torque += c.contact_point.cross(c.force_t) - pa.pos.cross(c.force_t);
                pa.torque += c.torque_r;
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
        buf.append("id,particle_a_id,particle_b_id,rest_length,extension,fx,fy,fz,active,heat_flow_a,damage,maximum_tensile_strain,elastic_energy,fracture_energy,dissipated_fracture_energy,source,maximum_compressive_strain,strength_scale,shear_broken,dissipated_shear_energy,break_mode,tangential_displacement,cohesive_budget_scale\n");
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
            AppendNumber(buf, b.strength_scale);                        buf.push_back(',');
            AppendNumber(buf, b.shear_broken ? 1 : 0);                  buf.push_back(',');
            AppendNumber(buf, b.dissipated_shear_energy);               buf.push_back(',');
            AppendNumber(buf, b.break_mode);                            buf.push_back(',');
            AppendNumber(buf, b.delta_t.norm());                        buf.push_back(',');
            AppendNumber(buf, b.cohesive_budget_scale);                 buf.push_back('\n');
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
