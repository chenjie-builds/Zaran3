/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file NSFieldSimulation.h
 * \brief NSFieldSimulation class, used to control the NS field simulation.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "BasicType.h"
#include "Visual.h"
#include "GlobalField.h"
#include "FieldDataCommInfo.h"
#include <iostream>
namespace zaran
{
    /// @brief 控制器类，用于控制流场求解
    /// @details 控制流场求解的前处理，求解，后处理，Feild之间的交互
    ///
    /// 阶段 3（DEM-CFD 单向耦合）把这层做成了**可继承的时间推进骨架**：
    /// 主流程（初始化 → 循环 { 预处理 / 推进一步 / 后处理 } → 收尾输出）保持不变，
    /// 只在三处留了耦合钩子。子类 `DEMCFDSimulation` 只实现这三个钩子，不重复时间循环。
    class NSFieldSimulation
    {
    public:
        explicit NSFieldSimulation(const shared_ptr<FieldManager> &field_manager);
        ~NSFieldSimulation();

    public:
        // 流场求解（虚函数：耦合算例继承后不需要重写时间推进）
        virtual void SolveField();

    protected:
        // 前处理
        void PreSolve();
        // 计算一步
        void SolveOneStep();
        // 后处理
        void PostSolve();

        // ------------------------------------------------------------
        // 耦合钩子（默认空实现 ⇒ 单相/两相算例的既有行为完全不变）
        // ------------------------------------------------------------
        /// @brief 求解开始前调用一次，且在 `Initialize()`（⇒ `solver->Init()`）**之前**。
        /// @details 供耦合算例注入 ε 场：外部 ε 必须在 InitField() 之前写好。
        virtual void PrepareCoupling() {}
        /// @brief 每个时间步**之前**调用（在 PreSolve 之前）。
        /// @details 双向耦合将在这里把上一时刻的相间反作用力散列进源项、
        ///          并按子循环推进 DEM；单向耦合下无事可做。
        virtual void CouplingPreStep(int iter) { (void)iter; }
        /// @brief 每个时间步**之后**调用（在 PostSolve 之后，流场输出已完成）。
        /// @details 单向耦合在这里用当前流场评估颗粒受力并输出。
        virtual void CouplingPostStep(int iter) { (void)iter; }

        /// @brief 取场管理器（子类需要访问场/求解器/数据管理器）
        shared_ptr<FieldManager> GetFieldManager() const { return m_field_manager; }

    protected:
		void CalcTimeStep();
        // 初始化，包括求解器初始化和流场初始化
        void Initialize() const;
        // 计算残差
        void CalcResidual();
        // field数据通信
        void CommFieldData();
        // 是否停止计算
        bool ContinueSolve();

    protected:
        // 输出流场数据
        void SaveFieldData();
        // 备份残差
        void SaveResidual() const;
        // 输出流场为Tecplot格式
        void SaveDataTecplot() const;
        // 备份流场数据，用于续算
        void BackupFieldData(std::string &back_folder);
        // 备份残差文件
        void BackupResidual(std::string &back_folder);
        // 备份日志文件
        void BackupLog(std::string &back_folder);
        // 备份全局参数
        void BackupGlobalData(std::string &back_folder);
    private:
        shared_ptr<FieldManager> m_field_manager;
        shared_ptr<Visual> m_visual;
        double m_res_max{};
        double m_res_ave{};
        bool m_res_flag;
    };

}