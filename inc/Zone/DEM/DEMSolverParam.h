/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file DEMSolverParam.h
 * \brief Solver parameters for DEM simulation, loaded from GlobalData.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "SolverPara.h"
#include "GlobalData.h"
#include "BasicType.h"
namespace zaran
{
    /// @brief DEM 求解参数，从 GlobalData（zaran.toml）读取
    class DEMSolverParam : public SolverParam
    {
    public:
        DEMSolverParam() = default;
        ~DEMSolverParam() override = default;

        void Init() override;

        // --- getters ---
        double              GetTimeStep()     const { return m_dt; }
        const Eigen::Vector3d& GetGravity()  const { return m_gravity; }
        const std::string&  GetContactModel() const { return m_contact_model; }
        const std::string&  GetParticleFile() const { return m_particle_file; }
        int                 GetOutputIter()   const { return m_output_iter; }
        int                 GetMaxIter()      const { return m_max_iter; }
        double              GetEndTime()      const { return m_end_time; }

        double GetYoungModulus()     const { return m_young_modulus; }
        double GetPoissonRatio()     const { return m_poisson_ratio; }
        double GetFrictionCoeff()    const { return m_friction_coeff; }
        double GetRestitutionCoeff() const { return m_restitution_coeff; }

    private:
        double         m_dt             = 1.0e-6;
        Eigen::Vector3d m_gravity       = Eigen::Vector3d(0.0, -9.81, 0.0);
        std::string    m_contact_model  = "HertzMindlin";
        std::string    m_particle_file  = "particles.csv";
        int            m_output_iter    = 100;
        int            m_max_iter       = 100000;
        double         m_end_time       = 1.0;

        double m_young_modulus     = 1.0e8;
        double m_poisson_ratio     = 0.3;
        double m_friction_coeff    = 0.4;
        double m_restitution_coeff = 0.9;
    };
} // namespace zaran
