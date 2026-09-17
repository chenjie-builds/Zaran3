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
        double maximum_compressive_strain = 0.0; ///< 历史最大压缩应变（取正值）
        double elastic_energy = 0.0;        ///< 当前键储存的弹性能 (J)
        double fracture_energy = 0.0;       ///< 本键断裂能阈值 (J)，0 表示未启用能量准则
        double dissipated_fracture_energy = 0.0; ///< 断裂时耗散的能量 (J)
        bool active = true;

        /// @brief 逐键强度折减系数（Weibull 异质性）。所有强度阈值（拉伸断裂应变、
        /// 软化峰值应变、压缩断裂应变）按 s 缩放，断裂能按 s² 缩放；
        /// s = 1 表示均质。由 dem.bond_weibull_modulus 抽样得到，抽样均值归一为 1，
        /// 因此不改变整体强度标定，只引入离散度。关闭异质性时恒为 1。
        double strength_scale = 1.0;

        /// @brief 连接来源标记：0 = 外部给定（输入 bonds.csv），
        /// 1 = 仿真初始时刻（t=0）按几何邻近关系自动生成的弹簧连接网络。
        /// 仅用于区分两种相互平行的机制与输出/后处理标记，不参与受力与断裂计算。
        int source = 0;
    };
} // namespace zaran
