/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file ContactModel.h
 * \brief Abstract base class for DEM contact models.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "DEMParticle.h"
#include "DEMContact.h"
namespace zaran
{
    /// @brief DEM 接触力模型抽象接口
    class ContactModel
    {
    public:
        virtual ~ContactModel() = default;

        /// @brief 计算法向接触力，结果写入 contact.force_n
        /// @param pa 粒子 A
        /// @param pb 粒子 B（若为墙面接触则用等效粒子）
        /// @param contact 接触对（含重叠量、法向量）
        /// @param dt 时间步长 (s)
        virtual void CalcNormalForce(const DEMParticle& pa, const DEMParticle& pb,
                                     DEMContact& contact, double dt) = 0;

        /// @brief 计算切向接触力，结果写入 contact.force_t，更新 contact.delta_t
        virtual void CalcTangentialForce(const DEMParticle& pa, const DEMParticle& pb,
                                         DEMContact& contact, double dt) = 0;
    };
} // namespace zaran
