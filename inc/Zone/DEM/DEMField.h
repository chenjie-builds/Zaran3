/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file DEMField.h
 * \brief DEMField class: Field specialization for DEM (no background grid).
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "Field.h"
#include "DEMFieldData.h"
#include "DEMSolver.h"
#include "DEMSolverParam.h"
#include "DEMDataManager.h"
namespace zaran
{
    /// @brief DEM 场，继承自 Field，不需要背景网格
    class DEMField : public Field
    {
    public:
        DEMField();
        ~DEMField() override = default;

        void Allocate() override;

        /// @brief 获取粒子数据容器
        shared_ptr<DEMFieldData>   GetDEMData()        const { return m_dem_data; }
        shared_ptr<DEMSolver>      GetDEMSolver()      const { return m_dem_solver; }
        shared_ptr<DEMSolverParam> GetDEMSolverParam() const { return m_dem_param; }
        shared_ptr<DEMDataManager> GetDEMDataManager() const { return m_dem_data_manager; }

    protected:
        void AllocateFieldData()   override;
        void AllocateSolver()      override;
        void AllocateSolverPara()  override;
        void AllocateDataManager() override;

    private:
        shared_ptr<DEMFieldData>   m_dem_data;
        shared_ptr<DEMSolver>      m_dem_solver;
        shared_ptr<DEMSolverParam> m_dem_param;
        shared_ptr<DEMDataManager> m_dem_data_manager;
    };
} // namespace zaran
