/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file DEMParticle.h
 * \brief DEMParticle class, representing a single DEM sphere particle.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "BasicType.h"
namespace zaran
{
    /// @brief 离散元球形粒子
    struct DEMParticle
    {
        index_type id    = 0;       ///< 全局粒子编号
        int        group = 0;       ///< 粒子分组（用于区分材料/边界等）

        // 几何量
        double radius = 1.0;        ///< 粒子半径 (m)

        // 力学量
        double mass    = 1.0;       ///< 质量 (kg)
        double inertia = 0.4;       ///< 转动惯量系数 I = inertia * m * r^2（球体 = 2/5）

        // 状态量
        Eigen::Vector3d pos   = Eigen::Vector3d::Zero();  ///< 质心位置 (m)
        Eigen::Vector3d vel   = Eigen::Vector3d::Zero();  ///< 平动速度 (m/s)
        Eigen::Vector3d rotation = Eigen::Vector3d::Zero(); ///< 累积转角 (rad)
        Eigen::Vector3d omega = Eigen::Vector3d::Zero();  ///< 角速度 (rad/s)

        // 力/力矩（每步清零后累积）
        Eigen::Vector3d force  = Eigen::Vector3d::Zero(); ///< 合力 (N)
        Eigen::Vector3d torque = Eigen::Vector3d::Zero(); ///< 合力矩 (N·m)

        // 材料参数
        double young_modulus     = 1.0e8;  ///< 杨氏模量 (Pa)
        double poisson_ratio     = 0.3;    ///< 泊松比
        double friction_coeff    = 0.4;    ///< 静/动摩擦系数
        double restitution_coeff = 0.9;    ///< 法向恢复系数

        // 热化学状态。字段采用有名称的成员，替代 LSM 中 heat/EC/EMC/decomposition
        // 等依赖数字下标的并行数组。
        bool   energetic = false;          ///< 是否为含能材料格点
        double temperature = 300.0;        ///< 固相温度 (K)，对应 heat(3,:)
        double reaction_progress = 0.0;    ///< 反应/分解质量分数 alpha，[0,1]
        double specific_heat = 1000.0;     ///< 定压/等效比热 (J/(kg K))
        double thermal_conductivity = 0.0; ///< 热导率 (W/(m K))
        double reaction_heat = 0.0;        ///< 单位质量完全反应热 (J/kg)
        double arrhenius_prefactor = 0.0;  ///< Arrhenius 指前因子 Z (1/s)
        double activation_temperature = 0.0; ///< Ea/R (K)，速率中使用 exp(-Ta/T)
        double heat_source = 0.0;          ///< 当前步其它体热源功率 (W)
        double reaction_rate = 0.0;        ///< 当前步 d(alpha)/dt (1/s)，用于输出
        int    phase = 0;                  ///< 0=固相，1=气固共存，2=反应产物
        double gas_temperature = 300.0;    ///< 共址反应产物温度 (K)
        double gas_pressure = 0.0;         ///< JWL 等效压力 (Pa)
        double reference_volume = 0.0;     ///< 初始控制体积 (m^3)
        double volume_ratio = 1.0;          ///< 当前/初始控制体积比
        double total_volume = 0.0;          ///< 当前局部总控制体积 (m^3)
        double solid_volume = 0.0;          ///< 气固共存点内剩余固相体积 (m^3)
        double gas_volume = 0.0;            ///< 气固共存点内产物气体体积 (m^3)
        double solid_core_radius = 0.0;      ///< 二维等效固相核半径 (m)
        double gas_radius = 0.0;             ///< 二维等效控制体半径 (m)
        double gas_internal_energy = 0.0;   ///< 产物气体内能 (J)
        double internal_heat_transfer = 0.0; ///< 本步从气相传给固相的热量 (J)
        double body_reaction_increment = 0.0; ///< 本步 Arrhenius 热解增量
        double core_burn_increment = 0.0;   ///< 本步气固共存固相核燃烧增量
        double neighbor_burn_increment = 0.0; ///< 本步相邻气相传播燃烧增量

        bool active = true; ///< 粒子是否参与计算（可用于标记固定边界粒子）
        bool kinematic = false; ///< 是否按给定速度运动而不受力加速度影响
        bool material_from_file = false; ///< 材料参数是否从文件显式指定（true 时 InitField 不覆盖）

        bool IsDynamic() const { return active && !kinematic; }
    };
} // namespace zaran
