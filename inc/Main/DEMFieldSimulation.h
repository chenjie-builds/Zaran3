/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file DEMFieldSimulation.h
 * \brief DEMFieldSimulation: main loop controller for DEM simulation.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "BasicType.h"
#include "GlobalField.h"
#include "DEMField.h"
namespace zaran
{
    /// @brief DEM 主循环控制器
    ///
    /// 流程：Init → 每步 Solve → 定期输出 VTP + 备份粒子 CSV
    class DEMFieldSimulation
    {
    public:
        explicit DEMFieldSimulation(shared_ptr<FieldManager> field_manager);
        ~DEMFieldSimulation() = default;

        void SolveField();

    protected:
        void Initialize();
        bool ContinueSolve() const;
        void SolveOneStep();
        void SaveFieldData(int iter) const;

    private:
        shared_ptr<FieldManager> m_field_manager;
        shared_ptr<DEMField>     m_dem_field; ///< 指向第一个 DEM 场的快捷指针
    };
} // namespace zaran
