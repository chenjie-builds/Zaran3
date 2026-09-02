/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file DEMSolver.h
 * \brief DEMSolver class: performs one DEM time step.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "FieldSolver.h"
#include "DEMFieldData.h"
#include "DEMSolverParam.h"
#include "ContactModel.h"
namespace zaran
{
    /// @brief DEM 求解器，继承 FieldSolver 基类
    ///
    /// 每调用一次 Solve() 推进一个 DEM 时间步：
    ///   ZeroForce → ContactDetection → CalcContactForce → CalcWallForce
    ///   → CalcGravity → Integrate
    class DEMSolver : public FieldSolver
    {
    public:
        DEMSolver(index_type index,
                  const std::string& name,
                  shared_ptr<DEMSolverParam> para,
                  shared_ptr<DEMFieldData>   dem_data);
        ~DEMSolver() override = default;

        // --- FieldSolver 接口 ---
        void Init()        override;
        void Solve()       override;
        void Preprocess()  override;
        void Postprocess() override;
        void InitSolver()  override;

        void InitField()                             override;
        void BoundaryCondition()                     override {}
        void UpdateField()                           override {}
        void BackupField(std::string& back_folder)   override;

        // --- 专用接口 ---
        /// @brief 备份粒子状态到文件（CSV）
        void BackupField(const std::string& back_folder) const;

        DEMFieldData*    GetDEMData()  const { return m_dem_data.get(); }
        DEMSolverParam*  GetDEMParam() const;

    protected:
        // 主流程各子步
        void ZeroForce();
        void ContactDetection();
        void CalcContactForce();
        void CalcWallForce();
        void CalcGravity();
        void Integrate();

    private:
        shared_ptr<DEMFieldData>   m_dem_data;
        unique_ptr<ContactModel>   m_contact_model;
    };
} // namespace zaran
