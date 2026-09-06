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
        double conduction_area = 0.0;      ///< 端点间等效导热面积 (m^2)，0 表示不导热

        Eigen::Vector3d delta_t = Eigen::Vector3d::Zero(); ///< 累积切向位移 (m)
        double extension = 0.0;                            ///< 当前轴向伸长 (m)
        Eigen::Vector3d force_a = Eigen::Vector3d::Zero(); ///< 当前施加于 A 的力 (N)
        double heat_flow_a = 0.0;          ///< 当前从 B 流向 A 的热功率 (W)
        double damage = 0.0;               ///< 不可逆双线性拉伸损伤 [0,1]
        double maximum_tensile_strain = 0.0; ///< 历史最大拉伸应变
        double elastic_energy = 0.0;        ///< 当前键储存的弹性能 (J)
        double fracture_energy = 0.0;       ///< 本键断裂能阈值 (J)，0 表示未启用能量准则
        double dissipated_fracture_energy = 0.0; ///< 断裂时耗散的能量 (J)
        bool active = true;
    };
} // namespace zaran
