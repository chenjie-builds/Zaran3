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

        // --- 刚性边界（把 prescribed-motion 粒子当作刚体平面）---
        /// @brief 是否将 kinematic（规定运动）粒子视为刚性边界平面
        bool GetRigidBoundaryEnabled() const { return m_rigid_boundary_enabled; }

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
    };
} // namespace zaran
