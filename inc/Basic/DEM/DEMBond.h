/**
 * Zaran - A Totally Automatic CFD Software
 * \file DEMBond.h
 * \brief Permanent elastic connection between two inert DEM particles.
 */
#pragma once
#include "BasicType.h"

namespace zaran
{
    /// @brief 永久键合弹簧；端点使用 DEMFieldData 中的粒子数组索引。
    struct DEMBond
    {
        index_type id = 0;
        index_type idx_a = 0;
        index_type idx_b = 0;

        double rest_length = 0.0;          ///< 无应力长度 (m)
        double normal_stiffness = 0.0;     ///< 法向刚度 (N/m)
        double tangential_stiffness = 0.0; ///< 切向刚度 (N/m)
        double normal_damping = 0.0;       ///< 法向阻尼 (N s/m)
        double tangential_damping = 0.0;   ///< 切向阻尼 (N s/m)

        Eigen::Vector3d delta_t = Eigen::Vector3d::Zero(); ///< 累积切向位移 (m)
        double extension = 0.0;                            ///< 当前轴向伸长 (m)
        Eigen::Vector3d force_a = Eigen::Vector3d::Zero(); ///< 当前施加于 A 的力 (N)
        bool active = true;
    };
} // namespace zaran
