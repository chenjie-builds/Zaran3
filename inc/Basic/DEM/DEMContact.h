/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file DEMContact.h
 * \brief DEMContact struct, representing a contact pair between two DEM entities.
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
    /// @brief 接触类型
    enum class ContactType
    {
        ParticleParticle, ///< 粒子-粒子接触
        ParticleWall,     ///< 粒子-墙面接触
    };

    /// @brief 接触对信息，包含接触几何量与历史切向位移（用于 Mindlin 模型）
    struct DEMContact
    {
        ContactType type = ContactType::ParticleParticle;

        index_type idx_a = 0; ///< 粒子 A 索引
        index_type idx_b = 0; ///< 粒子 B / 墙面索引

        // 接触几何
        double     overlap_n = 0.0;                          ///< 法向重叠量 δ_n (m)
        Eigen::Vector3d normal = Eigen::Vector3d::Zero();    ///< 单位法向量 n (A→B)
        Eigen::Vector3d contact_point = Eigen::Vector3d::Zero(); ///< 接触点坐标 (m)

        // 历史量（切向弹簧位移，用于 Mindlin 模型的增量切向力）
        Eigen::Vector3d delta_t = Eigen::Vector3d::Zero();   ///< 累积切向位移 (m)

        // 当前步接触力
        Eigen::Vector3d force_n = Eigen::Vector3d::Zero();   ///< 法向接触力（施加于 A）
        Eigen::Vector3d force_t = Eigen::Vector3d::Zero();   ///< 切向接触力（施加于 A）

        /// @brief 本步滚动阻力矩 (N·m)，**力偶**形式（施加于 A，B 取反）。
        /// 模型：M_r = −μ_r·R*·|F_n|·ω̂_rel（Ai 等 2011 的简化滚动阻力模型），
        /// 表示碎块因形状不规则/角锁而不能自由滚动 —— 缺了它，"碎块堆"会像
        /// 光滑球堆一样摊平，这正是"看起来太松散"的一个来源。
        /// 未启用（dem.rolling_friction ≤ 0）时恒为 0，不施加任何力矩。
        Eigen::Vector3d torque_r = Eigen::Vector3d::Zero();

        /// @brief 本步该接触对不可逆耗散的能量 (J)：法向阻尼耗散 + 切向摩擦耗散。
        /// 由接触模型写入（CalcNormalForce 先重置并累加法向项，CalcTangentialForce
        /// 再累加切向项），随后由 DEMSolver 对半分给两端并计入温度的能量账。
        /// 恒为非负（法向阻尼项 c_n·v_n²·dt ≥ 0；切向项取库仑截断前后的切向弹簧
        /// 储能之差，亦 ≥ 0）。
        double dissipation = 0.0;

        /// @brief 本步该接触对**可逆**储存在接触弹簧里的弹性能 (J)：
        /// ½·k_n·δ_n² + ½·k_t·|δ_t|²（切向项取库仑截断之后的值）。
        /// 由接触模型写入，供 DEM 求解器统计"接触弹性能"这一能量收支项。
        /// ⚠ 重叠保护开启时（dem.contact_stiffen_ratio > 0）实际法向力被非物理放大，
        /// 此时该值只是**名义**弹簧储能，放大带来的额外做功不在其中。
        double elastic_energy = 0.0;
    };
} // namespace zaran
