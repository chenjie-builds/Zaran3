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
        /// @brief 写 result/cavity_gas.csv（气腔状态 + 能量账本，每帧一行）
        void SaveCavityGasReport(int iter, const std::string& result_folder) const;
        /// @brief 写 result/fragment_stats.csv（碎块/裂纹带/径向断键率，起裂后每帧一行）
        void SaveFragmentStats(int iter, const std::string& result_folder) const;

    private:
        shared_ptr<FieldManager> m_field_manager;
        shared_ptr<DEMField>     m_dem_field; ///< 指向第一个 DEM 场的快捷指针
    };
} // namespace zaran
