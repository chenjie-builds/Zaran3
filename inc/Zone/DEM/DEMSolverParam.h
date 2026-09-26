/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file DEMSolverParam.h
 * \brief Solver parameters for DEM simulation, loaded from GlobalData.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "SolverPara.h"
#include "GlobalData.h"
#include "BasicType.h"
#include <cstdint>
namespace zaran
{
    /// @brief DEM 求解参数，从 GlobalData（zaran.toml）读取
    class DEMSolverParam : public SolverParam
    {
    public:
        DEMSolverParam() = default;
        ~DEMSolverParam() override = default;

        void Init() override;

        // --- getters ---
        double              GetTimeStep()     const { return m_dt; }
        const Eigen::Vector3d& GetGravity()  const { return m_gravity; }
        const std::string&  GetContactModel() const { return m_contact_model; }
        const std::string&  GetParticleFile() const { return m_particle_file; }
        int                 GetOutputIter()   const { return m_output_iter; }
        int                 GetMaxIter()      const { return m_max_iter; }
        double              GetEndTime()      const { return m_end_time; }

        double GetYoungModulus()     const { return m_young_modulus; }
        double GetPoissonRatio()     const { return m_poisson_ratio; }
        double GetFrictionCoeff()    const { return m_friction_coeff; }
        double GetRestitutionCoeff() const { return m_restitution_coeff; }
        double GetDensity()          const { return m_density; }
        bool GetReactionEnabled() const { return m_reaction_enabled; }
        double GetIgnitionTemperature() const { return m_ignition_temperature; }
        double GetMixedPhaseThreshold() const { return m_mixed_phase_threshold; }
        double GetBurnoutThreshold() const { return m_burnout_threshold; }
        double GetMaxTemperature() const { return m_max_temperature; }
        double GetJwlA() const { return m_jwl_a; }
        double GetJwlB() const { return m_jwl_b; }
        double GetJwlR1() const { return m_jwl_r1; }
        double GetJwlR2() const { return m_jwl_r2; }
        double GetJwlOmega() const { return m_jwl_omega; }
        double GetBurnA() const { return m_burn_a; }
        double GetBurnB() const { return m_burn_b; }
        double GetBurnPressureThreshold() const { return m_burn_pressure_threshold; }
        double GetBurnTemperatureThreshold() const { return m_burn_temperature_threshold; }
        double GetGasTemperatureReference() const { return m_gas_temperature_reference; }
        double GetPressureCap() const { return m_pressure_cap; }
        double GetBondBreakStrain() const { return m_bond_break_strain; }
        double GetBondPeakStrain() const { return m_bond_peak_strain; }
        double GetReactionWeakening() const { return m_reaction_weakening; }
        const Eigen::Vector3d& GetIgnitionCenter() const { return m_ignition_center; }
        double GetIgnitionRadius() const { return m_ignition_radius; }
        double GetIgnitionCenterTemperature() const { return m_ignition_center_temperature; }
        double GetLatticeThickness() const { return m_lattice_thickness; }
        double GetInternalHeatLimit() const { return m_internal_heat_limit; }
        bool GetVoronoiEnabled() const { return m_voronoi_enabled; }
        int GetVoronoiUpdateInterval() const { return m_voronoi_update_interval; }
        double GetVoronoiMaxVolumeRatio() const { return m_voronoi_max_volume_ratio; }
        double GetBoxXMin() const { return m_box_x_min; }
        double GetBoxXMax() const { return m_box_x_max; }
        double GetBoxYMin() const { return m_box_y_min; }
        double GetBoxYMax() const { return m_box_y_max; }
        double GetSurfaceEnergy() const { return m_surface_energy; }
        double GetGrainBoundaryEnergy() const { return m_grain_boundary_energy; }
        double GetHighPressureRnn() const { return m_high_pressure_rnn; }
        double GetHighPressureKnn() const { return m_high_pressure_knn; }
        double GetHighPressureExponent() const { return m_high_pressure_exponent; }
        double GetNormalViscosity() const { return m_normal_viscosity; }
        double GetTangentialViscosity() const { return m_tangential_viscosity; }

        // --- 初始弹簧连接网络（可选，用于脆性材料断裂/破坏；与逐步接触检测相互独立）---
        /// @brief 是否在仿真开始时按几何邻近关系自动建立弹簧连接网络
        bool GetSpringNetworkEnabled() const { return m_spring_network_enabled; }
        /// @brief 建链判据的相对间隙容差：dist ≤ (r_a+r_b)·(1+gap)
        double GetSpringNetworkGap() const { return m_spring_network_gap; }
        /// @brief 弹簧法向刚度 (N/m)；≤0 表示按接触模型的等效刚度自动推导
        double GetSpringNetworkStiffness() const { return m_spring_network_stiffness; }
        /// @brief 弹簧切向刚度 (N/m)；≤0 表示取法向刚度的一半
        double GetSpringNetworkTangentialStiffness() const { return m_spring_network_tangential_stiffness; }
        /// @brief 弹簧断裂应变阈值；≤0 表示回退到 dem.bond_break_strain
        double GetSpringNetworkFractureStrain() const { return m_spring_network_fracture_strain; }

        // --- 键合强度：压缩失效 + 逐键异质性 ---
        /// @brief 抗压/抗拉强度比 σ_c/σ_t（脆性材料典型 8–15）。线性弹性下
        ///        压缩断裂应变 = ratio × 抗拉断裂应变。0 表示不做压缩失效判定。
        double GetBondCompressionStrengthRatio() const { return m_bond_compression_strength_ratio; }
        /// @brief 直接指定压缩断裂应变；>0 时优先于 ratio（两者都按逐键强度系数缩放）
        double GetBondBreakStrainCompression() const { return m_bond_break_strain_compression; }
        /// @brief Weibull 形状参数（模量）m；0 表示关闭异质性（所有键共用同一阈值）。
        ///        抽样归一化到均值为 1，故不改变整体强度标定，只引入离散度。
        double GetBondWeibullModulus() const { return m_bond_weibull_modulus; }
        /// @brief Weibull 抽样种子（按 (seed, id_a, id_b) 哈希，保证可复现且与键序无关）
        std::uint64_t GetBondWeibullSeed() const { return m_bond_weibull_seed; }

        // --- γ 控制的双线性内聚律（断裂能-强度自洽）---
        /// @brief 是否启用内聚模式。启用条件是：本开关为真 **且** 该键 `fracture_energy > 0`
        ///        （即给了 dem.surface_energy）**且** `dem.bond_peak_strain` 有限（<1e29）。
        ///        三者缺一即回退到原有的"储能达断裂能即脆断"路径（保证历史算例逐位不变）。
        bool GetBondCohesiveEnabled() const { return m_bond_cohesive_enabled; }
        /// @brief 内聚模式的严格性：δ_f < δ_p 时是报错还是仅告警并退化为脆断。
        ///        默认严格（报错），因为静默退化正是最难发现的一类错误。
        bool GetBondCohesiveStrict() const { return m_bond_cohesive_strict; }
        /// @brief 软化段至少要有多少个键长（δ_f − δ_p ≥ ratio·L0）才认为网格分辨得开裂纹；
        ///        不足时告警（不报错）。0 表示不检查。
        double GetBondCohesiveMinSofteningBonds() const { return m_bond_cohesive_min_bonds; }

        // --- 多通道失效（切向独立断裂 / 压剪），对照 LSM breakmod 2/4 ---
        /// @brief 切向独立断裂能比 α = UtIII/UnIII。α > 0 时启用切向通道：
        ///        切向弹簧储能达 UtIII 即切向先断，并把 UtIII 从法向断裂能里扣除；
        ///        压缩侧的同一阈值给出"压剪破坏"通道。α ≤ 0（默认）表示关闭，
        ///        所有既有一致性逐位不变。
        double GetBondShearEnergyRatio() const { return m_bond_shear_energy_ratio; }

        // --- 碎后接触：滚动阻力矩（决定"碎块堆是否像光滑球堆一样摊平"）---
        /// @brief 滚动阻力系数 μ_r（长度量纲）。滚动阻力矩 M_r = −μ_r·R*·|F_n|·ω̂_rel，
        ///        按 I_eff·|ω_rel|/dt 截断以防过冲反向。≤ 0（默认）表示关闭。
        double GetRollingFriction() const { return m_rolling_friction; }

        // --- 接触切向阻尼（法向有按恢复系数标定的 c_n，切向此前完全没有）---
        /// @brief 切向阻尼与法向阻尼之比 λ：c_t = λ·c_n。λ ≤ 0（默认）表示关闭。
        ///        按同恢复系数标定时 c_t/c_n = sqrt(k_t/k_n)（线性模型 ≈ 0.707）。
        double GetTangentialDampingScale() const { return m_tangential_damping_scale; }

        // --- 刚性边界（把 prescribed-motion 粒子当作刚体平面）---
        /// @brief 是否将 kinematic（规定运动）粒子视为刚性边界平面
        bool GetRigidBoundaryEnabled() const { return m_rigid_boundary_enabled; }

        // --- 接触重叠限制：不允许两个粒子重叠过近（参考 LSM rn_limit / rn_rebound）---
        /// @brief 过深压缩时放大法向排斥力的阈值，以 (r_a+r_b) 的倍数给出；
        ///        当 dist < ratio·(r_a+r_b) 时把法向排斥力乘以 (ratio·(r_a+r_b)/dist)^20。
        ///        ≤0 表示关闭该机制。
        double GetContactStiffenRatio() const { return m_contact_stiffen_ratio; }
        /// @brief 深压缩刚性回弹阈值，以 (r_a+r_b) 的倍数给出；
        ///        当 dist < ratio·(r_a+r_b) 且两端仍相向运动时，
        ///        基于动量守恒按（准）弹性碰撞直接改写两端的法向速度。
        ///        ≤0 表示关闭该机制。
        double GetContactReboundRatio() const { return m_contact_rebound_ratio; }

        // --- 机械耗散生热 ---
        /// @brief 是否把机械耗散（接触法向阻尼、接触切向摩擦、键合阻尼）计入温度。
        /// 关闭后温度只由键合导热与外部体热源驱动（与旧行为一致）。
        bool GetMechanicalHeatingEnabled() const { return m_mechanical_heating; }

        // --- 高压气腔加载（孔洞内充压气体把脆性材料撑碎）---
        /// @brief 是否启用气腔加载
        bool GetGasCavityEnabled() const { return m_gas_cavity_enabled; }
        /// @brief 气腔中心 (m)
        const Eigen::Vector3d& GetGasCavityCenter() const { return m_gas_cavity_center; }
        /// @brief 初始气腔半径 a0 (m)
        double GetGasCavityRadius() const { return m_gas_cavity_radius; }
        /// @brief 初始气体压力 p0 (Pa)
        double GetGasPressureInitial() const { return m_gas_pressure_initial; }
        /// @brief 气体绝热（多变）指数 gamma；1.4 = 双原子理想气体
        double GetGasPolytropicIndex() const { return m_gas_polytropic_index; }
        /// @brief 升压时间 (s)；0 = 瞬时加载
        double GetGasRampTime() const { return m_gas_ramp_time; }
        /// @brief 气腔边界环识别容差（× 格距）
        double GetGasCavityShellTolerance() const { return m_gas_cavity_shell_tolerance; }
        /// @brief 压力上限 (Pa)；≤0 = 不限
        double GetGasCavityPressureCap() const { return m_gas_cavity_pressure_cap; }
        /// @brief 泄压判据：气腔面积 / (π·样品外半径²)；≤0 = 永不泄压
        double GetGasCavityVentAreaRatio() const { return m_gas_cavity_vent_area_ratio; }
        /// @brief 样品外半径 (m)，供泄压判据使用
        double GetSampleRadius() const { return m_sample_radius; }

    private:
        double         m_dt             = 1.0e-6;
        Eigen::Vector3d m_gravity       = Eigen::Vector3d(0.0, -9.81, 0.0);
        std::string    m_contact_model  = "HertzMindlin";
        std::string    m_particle_file  = "particles.csv";
        int            m_output_iter    = 100;
        int            m_max_iter       = 100000;
        double         m_end_time       = 1.0;

        double m_young_modulus     = 1.0e8;
        double m_poisson_ratio     = 0.3;
        double m_friction_coeff    = 0.4;
        double m_restitution_coeff = 0.9;
        double m_density             = 2500.0;
        bool   m_reaction_enabled     = false;
        double m_ignition_temperature = 500.0;
        double m_mixed_phase_threshold = 1.0e-3;
        double m_burnout_threshold    = 0.99;
        double m_max_temperature      = 1300.0;
        // HMX 基 PBX 的 A/B/R1/R2/omega 默认值来自 LSM gas_phase_mod.f90 注释。
        double m_jwl_a = 8.545e11;
        double m_jwl_b = 0.2493e11;
        double m_jwl_r1 = 4.60;
        double m_jwl_r2 = 1.35;
        double m_jwl_omega = 0.25;
        double m_burn_a = 0.639;
        double m_burn_b = 1.19;
        double m_burn_pressure_threshold = 9.0e4;
        double m_burn_temperature_threshold = 1000.0;
        double m_gas_temperature_reference = 3000.0;
        double m_pressure_cap = 4.0e10;
        double m_bond_break_strain = 1.0e30;
        double m_bond_peak_strain = 1.0e30;
        double m_reaction_weakening = 0.0;
        Eigen::Vector3d m_ignition_center = Eigen::Vector3d::Zero();
        double m_ignition_radius = 0.0;
        double m_ignition_center_temperature = 0.0;
        double m_lattice_thickness = 1.0;
        double m_internal_heat_limit = 0.25;
        bool m_voronoi_enabled = false;
        int m_voronoi_update_interval = 20;
        double m_voronoi_max_volume_ratio = 4.0;
        double m_box_x_min = -1.0;
        double m_box_x_max = 1.0;
        double m_box_y_min = -1.0;
        double m_box_y_max = 1.0;
        double m_surface_energy = 0.0;
        double m_grain_boundary_energy = 0.0;
        double m_high_pressure_rnn = 0.0;
        double m_high_pressure_knn = 0.0;
        double m_high_pressure_exponent = 0.0;
        double m_normal_viscosity = 0.0;
        double m_tangential_viscosity = 0.0;

        // 初始弹簧连接网络（默认关闭，不影响既有算例）
        bool   m_spring_network_enabled = false;
        double m_spring_network_gap = 0.02;
        double m_spring_network_stiffness = 0.0;
        double m_spring_network_tangential_stiffness = 0.0;
        double m_spring_network_fracture_strain = 0.02;

        // 刚性边界：把 prescribed-motion（kinematic）粒子当作刚体平面。
        // 关闭时保持原行为（球-球接触，法向随中心连线翻转）。
        bool m_rigid_boundary_enabled = false;

        // 键合强度：压缩失效（默认按 σ_c/σ_t = 8 由抗拉阈值派生）与 Weibull 异质性（默认关闭）。
        double        m_bond_compression_strength_ratio = 8.0;
        double        m_bond_break_strain_compression = 0.0;
        double        m_bond_weibull_modulus = 0.0;
        std::uint64_t m_bond_weibull_seed = 20260917ULL;
        bool          m_bond_cohesive_enabled = true;
        bool          m_bond_cohesive_strict = true;
        double        m_bond_cohesive_min_bonds = 1.0;

        // 多通道失效 / 碎后接触（全部默认关闭 ⇒ 历史算例逐位不变）
        double        m_bond_shear_energy_ratio = 0.0;   ///< α = UtIII/UnIII；≤0 = 关闭
        double        m_rolling_friction = 0.0;          ///< μ_r（长度量纲）；≤0 = 关闭
        double        m_tangential_damping_scale = 0.0;  ///< c_t = λ·c_n；≤0 = 关闭

        // 接触重叠限制。(r_a+r_b) 对应 LSM 的平衡间距 r0，默认取 LSM 的 0.6 作为统一阈值。
        //   ① 排斥力放大：dist < ratio·(r_a+r_b) 时 F_n *= (ratio·(r_a+r_b)/dist)^20
        //      （对应 LSM rn_limit = 0.6·r0）。从阈值处平滑起步、被回弹在同一阈值截住，
        //      实际倍数远小于上限；实测巴西盘算例开启后 max|v| 与关闭时一致，无额外动能注入。
        //   ② 刚性回弹：dist < ratio·(r_a+r_b) 且相向运动 → 按弹性碰撞改写法向速度
        //      （对应 LSM rn_rebound，LSM 取 0.5·r0）。这是硬约束：把最近接近距离
        //      截在约 ratio·(r_a+r_b)，接触法向不会翻转。
        // 阈值为何取 0.6 而非 LSM 的 0.5：当**加载边界由离散球排构成**时（如巴西盘的
        // 位移控制压板，球心间距 = 2·R_boundary），试件粒子可从相邻边界球的尖角缝隙穿过，
        // 除非下限超过 R_boundary/(R_boundary + R_particle)（本工程算例 = 0.533）。
        // 实测 ratio=0.5 时 10 万步后仍逃逸 30/30 个粒子，ratio=0.6 则为 0/0。
        // 若想严格复刻 LSM，可显式设 contact_stiffen_ratio=0.6 / contact_rebound_ratio=0.5。
        // 两项均以 0 表示关闭；同时启用时要求 回弹阈值 ≤ 放大阈值。
        double m_contact_stiffen_ratio  = 0.6;
        double m_contact_rebound_ratio  = 0.6;

        // 机械耗散（接触法向阻尼 / 切向摩擦 / 键合阻尼）是否计入温度，默认开启。
        bool m_mechanical_heating = true;

        // 高压气腔加载（默认关闭）。
        //   物理：空腔内的气体以压力 p 顶在空腔壁上，壁面由**空腔边界粒子环**离散，
        //         每个边界粒子承受 F_i = p·ℓ_i·t·n̂_i（ℓ_i 为其角向分片对应的弦长、
        //         t 为平面外厚度）。气腔体积随边界位移按散度定理增长，气体按绝热
        //         （多变）规律膨胀做功、压力下降。
        //   为什么必须单独实现：DEM 里既有的气相压力（CalcGasPressureForce）只在
        //         reaction_enabled 打开、且气体由 Arrhenius/燃烧反应**产生**时才生效；
        //         本算例的气体是**初始就存在**于空腔中的高压气体，且加载面是自由表面
        //         （空腔壁），而不是两个固相粒子之间的内部连接面。
        bool   m_gas_cavity_enabled = false;
        Eigen::Vector3d m_gas_cavity_center = Eigen::Vector3d::Zero();
        double m_gas_cavity_radius = 0.0;
        double m_gas_pressure_initial = 0.0;
        double m_gas_polytropic_index = 1.4;
        double m_gas_ramp_time = 0.0;
        double m_gas_cavity_shell_tolerance = 0.5;
        double m_gas_cavity_pressure_cap = 0.0;   ///< ≤0 = 不限
        double m_gas_cavity_vent_area_ratio = 0.0;///< ≤0 = 永不泄压
        double m_sample_radius = 0.0;
    };
} // namespace zaran
