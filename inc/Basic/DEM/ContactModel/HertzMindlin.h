/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file HertzMindlin.h
 * \brief Hertz-Mindlin nonlinear contact model for DEM.
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
    /// @brief 赫兹-敏德林（Hertz-Mindlin）非线性接触力模型
    ///
    /// 法向力：F_n = (4/3) E* sqrt(R*) δ_n^{3/2}  + 阻尼项
    /// 切向力：增量切向弹簧 + Coulomb 摩擦截断
    class HertzMindlin : public ContactModel
    {
    public:
        HertzMindlin() = default;
        ~HertzMindlin() override = default;

        void CalcNormalForce(const DEMParticle& pa, const DEMParticle& pb,
                             DEMContact& contact, double dt) override;

        void CalcTangentialForce(const DEMParticle& pa, const DEMParticle& pb,
                                 DEMContact& contact, double dt) override;

    private:
        /// @brief 计算等效弹性模量 E*、等效半径 R*、等效剪切模量 G*
        void CalcEffectiveParams(const DEMParticle& pa, const DEMParticle& pb,
                                 double& E_star, double& R_star, double& G_star) const;
    };
} // namespace zaran
