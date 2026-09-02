/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file DEMWall.h
 * \brief DEMWall class, representing a planar boundary wall for DEM simulations.
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
    /// @brief 简单平面墙（法向量 + 点确定平面）
    struct DEMWall
    {
        index_type id = 0;              ///< 墙面编号
        Eigen::Vector3d normal;         ///< 外法向（指向流体域）
        Eigen::Vector3d point;          ///< 平面上一点

        // 材料参数（墙面侧）
        double young_modulus     = 1.0e9;
        double poisson_ratio     = 0.3;
        double friction_coeff    = 0.4;
        double restitution_coeff = 0.9;

        DEMWall()
            : normal(Eigen::Vector3d(0.0, 1.0, 0.0)),
              point(Eigen::Vector3d::Zero())
        {}

        /// @brief 计算粒子质心到平面的有符号距离（正值=粒子在墙外侧）
        double SignedDist(const Eigen::Vector3d& pos) const
        {
            return normal.dot(pos - point);
        }
    };
} // namespace zaran
