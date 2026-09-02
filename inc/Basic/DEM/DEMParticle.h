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
        Eigen::Vector3d omega = Eigen::Vector3d::Zero();  ///< 角速度 (rad/s)

        // 力/力矩（每步清零后累积）
        Eigen::Vector3d force  = Eigen::Vector3d::Zero(); ///< 合力 (N)
        Eigen::Vector3d torque = Eigen::Vector3d::Zero(); ///< 合力矩 (N·m)

        // 材料参数
        double young_modulus     = 1.0e8;  ///< 杨氏模量 (Pa)
        double poisson_ratio     = 0.3;    ///< 泊松比
        double friction_coeff    = 0.4;    ///< 静/动摩擦系数
        double restitution_coeff = 0.9;    ///< 法向恢复系数

        bool active = true; ///< 粒子是否参与计算（可用于标记固定边界粒子）
        bool material_from_file = false; ///< 材料参数是否从文件显式指定（true 时 InitField 不覆盖）
    };
} // namespace zaran
