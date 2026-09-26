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

        // ---------------- γ 控制的双线性内聚律（γ 与峰值应变同时给定才启用）----------------
        // 牵引-分离曲线：δ ≤ δ_p 线性上升（斜率 kn），δ_p < δ < δ_f 线性软化到 0，
        // δ ≥ δ_f 键失效。构造成**曲线下面积恒等于 fracture_energy**：
        //   δ_f = 2·fracture_energy·s² / (kn·δ_p)   （solid_area_scale 已约掉）
        // 于是断裂时耗散的能量恰好等于材料断裂能，不再是一个由网格决定的派生量。
        // 该模式下 damage 由 δ 的解析式给出（不是峰值-断裂应变的线性插值），
        // 因此 "键对外做的功 = 弹性能 + 耗散能" 逐位成立。见 LSM
        // covalent_ionic_spring_mod.f90（nn=1、β=1 时即为本式）与
        // read_material_properties.f90 的 rnmax 反解。
        bool   cohesive = false;                    ///< 是否走 γ 控制内聚律
        double cohesive_peak_separation = 0.0;      ///< 峰值分离度 δ_p (m)，含 strength_scale
        double cohesive_failure_separation = 0.0;   ///< 失效分离度 δ_f (m)，由断裂能反解

        // ---------------- 多通道失效（对照 LSM covalent_ionic_spring_mod.f90:489-527）----------------
        // LSM 每个键有四条互不相同的断裂通道（breakmod 1/2/3/4）：
        //   1 拉伸混合（En + Et ≥ UnIII）、2 压剪（rn ≤ r0 且 Et ≥ UtIII）、
        //   3 压缩压溃（rn < rnmax·0.6）、4 切向弹簧独立断裂（Et ≥ UtIII ⇒ tangent=.false.，
        //   并把剩余法向阈值按 GBratio 折减）。
        // 本工程原有的实现只有"拉伸（法向能量/应变）"与"压缩压溃"两条；下面补齐 2 与 4。
        //
        // 切向断裂能的守恒写法：LSM 折减规则
        //     GBratio_new = (UnIII − UtIII)·GBratio_old / UnIII
        // 代入 Un_eff = UnIII·GBratio 得
        //     Un_eff_new = Un_eff_old − Ut_eff
        // 即"切向先断时，把造切向新表面花掉的那份断裂能从法向阈值里扣掉"。
        // 于是全过程耗散（切向 + 法向）仍然等于材料分配给该键的总断裂能，账本自动闭合。

        /// @brief 切向独立断裂能阈值 UtIII (J)，来自 dem.bond_shear_energy_ratio × UnIII。
        /// **恒为 0 表示关闭切向通道**（此时下面整条多通道路径逐位回退到历史行为）。
        double fracture_energy_shear = 0.0;
        /// @brief 法向**断裂能预算**的折减系数（即 LSM 的 GBratio），初值 1。
        /// 切向独立断裂时按 LSM 规则乘上 (1 − Ut0/Un0)：
        ///     GBratio_new = (UnIII − UtIII)·GBratio_old / UnIII
        /// ⚠ **它折减的是能量阈值，不是 δ_f**。LSM 里 `rnmax` 是材料常数、全程不改动，
        /// 判据用的是 `UnIII·GBratio`。第一版把 δ_f 按面积重解，结果在"法向通道已经
        /// 吸收了部分包络功之后才切向断裂"时，会把那段已吸收的功从账上抹掉
        /// （实测丢 0.0103 J / 相对 1.3e-5）—— 改成只折减阈值后账本逐帧闭合。
        double cohesive_budget_scale = 1.0;
        /// @brief 切向弹簧是否已独立断裂（LSM breakmod=4）。断裂后切向不再承载，
        /// 且 cohesive_budget_scale 已按 LSM 规则折减过一次。
        bool   shear_broken = false;
        /// @brief 切向独立断裂时耗散掉的能量 (J)，与 dissipated_fracture_energy（法向部分）
        /// 相加才是该键的断裂耗散总量。
        double dissipated_shear_energy = 0.0;
        /// @brief 失效模式标记（LSM 的 break_mod，仅诊断用，不参与受力）：
        /// 0 = 未失效，1 = 拉伸混合，2 = 压剪，3 = 压缩压溃，4 = 切向独立断裂。
        int    break_mode = 0;

        /// @brief 键中点到空腔中心的**初始**半径 (m)，建键时赋值一次。
        /// 径向损伤剖面必须用它、而不是当前位置：空腔膨胀/碎片飞出后，
        /// 按当前位置统计会让环带跟着漂移，读出来的"外缘全断"其实是飞出去的碎块。
        double initial_radius = -1.0;

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
