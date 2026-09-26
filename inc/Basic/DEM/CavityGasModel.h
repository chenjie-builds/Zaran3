/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file CavityGasModel.h
 * \brief 高压气腔（孔洞内充压气体）驱动脆性材料破碎的加载模型。
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "BasicType.h"
#include "CommonPara.h"
#include "DEMParticle.h"
#include <vector>

namespace zaran
{
    /// @brief 气腔模型的可调参数（全部来自 zaran.toml 的 `[dem]` 段）
    struct CavityGasOptions
    {
        bool   enabled = false;                 ///< 总开关
        Eigen::Vector3d center = Eigen::Vector3d::Zero(); ///< 空腔中心 (m)
        double cavity_radius = 0.0;             ///< 初始空腔半径 a0 (m)
        double pressure_initial = 0.0;          ///< 初始气体压力 p0 (Pa)
        double polytropic_index = 1.4;          ///< 绝热（多变）指数 gamma
        double ramp_time = 0.0;                 ///< 升压时间 (s)；0 = 瞬时加载
        double thickness = 1.0;                 ///< 平面外厚度 (m)

        /// @brief 边界环识别容差：d_i ≤ a0 + r_p + shell_tolerance·L0
        double shell_tolerance = 0.5;
        /// @brief 边界环识别上限：d_i ≤ max_ring_radius_ratio·a0（防止容差过大时抓到远处粒子）
        double max_ring_radius_ratio = 2.0;

        /// @brief 样品外半径 (m)；仅用于泄压判据
        double disk_radius = 0.0;
        /// @brief 气腔面积超过 vent_area_ratio·πR² 时判定为"已经破到自由面"，压力归零；
        ///        ≤0 表示不做泄压判定（气体永久封闭）。
        double vent_area_ratio = 0.0;

        /// @brief 压力上限 (Pa)；≤0 表示不限制
        double pressure_cap = 0.0;
    };

    /// @brief 高压气腔加载模型
    ///
    /// ## 物理模型
    /// 空腔内的气体对空腔壁施加法向压力 p，壁面由**空腔边界粒子环**离散：
    /// 第 i 个边界粒子占据角向扇区 Δθ_i、内切半径 ρ_i = d_i − r_p，故其暴露面段
    /// 长度为 ℓ_i = ρ_i·Δθ_i，外法向 n̂_i 由腔心指向该粒子。气体作用力
    ///
    ///     F_i = p · ℓ_i · t · n̂_i        （t = 平面外厚度）
    ///
    /// ## 体积与做功的严格自洽
    /// 同一个边界离散同时给出空腔面积的变化率（二维散度定理）：
    ///
    ///     ΔA = Σ_i ℓ_i (Δu_i · n̂_i)，   ΔV = ΔA·t
    ///
    /// 于是 Σ_i F_i·Δu_i = p·t·Σ_i ℓ_i(Δu_i·n̂_i) = p·ΔV —— **气体做功与体积变化
    /// 逐位一致**，与内能更新 U ← U − W 合起来给出严格闭合的能量账本。
    ///
    /// ## 状态方程
    /// 气体状态用内能 U 表征（U = pV/(γ−1)），每步按绝热做功更新：
    /// U ← U − pΔV，再由 p = (γ−1)U/V 反解压力。连续极限下即 pV^γ = const。
    /// `ramp_time > 0` 时升压段内压力按 p0·t/t_ramp 强制给定，注入的能量单独计入
    /// `external_work`，能量账本仍然闭合。
    class CavityGasModel
    {
    public:
        CavityGasModel() = default;
        ~CavityGasModel() = default;

        /// @brief 识别边界环、建立角向分片、由几何给出初始空腔面积/体积与气体内能
        /// @param lattice_spacing 键合格距（用于边界环识别的容差）；≤0 时退化为 2·粒子半径
        void Initialize(const dynamic_array<DEMParticle>& particles,
                        const CavityGasOptions& options,
                        double lattice_spacing);

        bool IsEnabled()  const { return m_options.enabled; }
        bool IsVented()   const { return m_vented; }

        /// @brief 步首（积分之前）：按当前位置刷新环几何并施加气体压力
        void ApplyPressureForce(dynamic_array<DEMParticle>& particles, double time);

        /// @brief 步末（积分之后）：用**与受力完全相同**的边界离散更新体积与内能
        void UpdateState(const dynamic_array<DEMParticle>& particles, double dt, double time);

        // ---------------- 诊断 ----------------
        double GetPressure()            const { return m_pressure; }
        double GetVolume()              const { return m_volume; }
        double GetArea()                const { return m_volume / m_thickness; }
        double GetAreaGeometric()       const { return m_area_geom; }
        double GetEffectiveRadius()     const { return std::sqrt(std::max(0.0, GetArea()) / PI); }
        double GetInternalEnergy()      const { return m_internal_energy; }
        double GetInitialInternalEnergy() const { return m_internal_energy_0; }
        double GetInitialVolume()       const { return m_volume_0; }
        double GetGasWork()             const { return m_gas_work; }
        double GetExternalWork()        const { return m_external_work; }
        double GetRingMaxGap()          const { return m_ring_max_gap; }
        double GetAppliedForceSum()     const { return m_force_sum; }
        double GetVentTime()            const { return m_vent_time; }
        index_type GetRingCount()       const { return static_cast<index_type>(m_ring.size()); }
        const std::vector<index_type>& GetRingIndices() const { return m_ring; }

    private:
        /// @brief 按当前位置刷新 θ_i / Δθ_i / ρ_i / ℓ_i / n̂_i，并给出几何面积
        void UpdateRingGeometry(const dynamic_array<DEMParticle>& particles);
        /// @brief 本步使用的压力：升压段按时间给定，否则由内能反解
        double CurrentPressure(double time) const;

        CavityGasOptions m_options;
        double m_thickness = 1.0;

        std::vector<index_type> m_ring;   ///< 边界环粒子索引（t=0 确定，之后恒定）
        std::vector<double> m_theta;      ///< 当前角位置（按环序单调展开）
        std::vector<double> m_dtheta;     ///< 角向分片宽度
        std::vector<double> m_rho;        ///< 内切半径 d_i − r_p
        std::vector<double> m_chord;      ///< 面段长度 ℓ_i
        std::vector<double> m_normal;     ///< 单位外法向（3×ring）

        double m_volume = 0.0;
        double m_volume_0 = 0.0;
        double m_internal_energy = 0.0;
        double m_internal_energy_0 = 0.0;
        double m_pressure = 0.0;
        double m_gas_work = 0.0;
        double m_external_work = 0.0;
        double m_area_geom = 0.0;
        double m_ring_max_gap = 0.0;
        double m_force_sum = 0.0;

        bool   m_vented = false;
        double m_vent_time = -1.0;
        bool   m_warned_ring = false;
    };
} // namespace zaran
