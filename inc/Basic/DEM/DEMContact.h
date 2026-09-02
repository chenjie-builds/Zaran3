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
    };
} // namespace zaran
