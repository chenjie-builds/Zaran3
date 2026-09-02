/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file LinearSpringDashpot.h
 * \brief Linear Spring-Dashpot (LSD) contact model for DEM.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "ContactModel.h"
#include <cmath>
namespace zaran
{
    /// @brief 线弹簧阻尼（Linear Spring-Dashpot）接触力模型
    ///
    /// 法向力：F_n = k_n * δ_n + c_n * v_n_rel（n 方向）
    /// 切向力：F_t = -k_t * δ_t  （Coulomb 摩擦截断）
    class LinearSpringDashpot : public ContactModel
    {
    public:
        LinearSpringDashpot() = default;
        ~LinearSpringDashpot() override = default;

        void CalcNormalForce(const DEMParticle& pa, const DEMParticle& pb,
                             DEMContact& contact, double dt) override;

        void CalcTangentialForce(const DEMParticle& pa, const DEMParticle& pb,
                                 DEMContact& contact, double dt) override;

    private:
        /// @brief 根据材料参数计算等效法向刚度 k_n 与阻尼系数 c_n
        void CalcNormalStiffness(const DEMParticle& pa, const DEMParticle& pb,
                                 double& k_n, double& c_n) const;
    };
} // namespace zaran
