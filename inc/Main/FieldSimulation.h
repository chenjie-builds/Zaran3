//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Controller.h														||
//*	@brief	Control solver to modify field_manager data in grid							||
//*	@author	Chen Jie.															||
//==============================================================================||
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
    class NSFieldSimulation
    {
    public:
        NSFieldSimulation(shared_ptr<FieldManager> global_Field);
        ~NSFieldSimulation();

    public:
        // 流场求解
        void SolveField();

    protected:
        // 前处理
        void PreSolve();
        // 计算一步
        void SolveOneStep();
        // 后处理
        void PostSolve();

    protected:
		void CalcTimeStep();
        // 初始化，包括求解器初始化和流场初始化
        void Initialize();
        // 计算残差
        void CalcResidual();
        // field数据通信
        void CommFieldData();
        // 是否停止计算
        bool ContinueSolve();

    protected:
        // 输出流场数据
        void SaveFieldData();
        // 备份边界节点数据
        void SaveWallNode();
        // 备份残差
        void SaveResidual() const;
        // 输出流场为VTK格式
        void SaveDataVTK(std::ostream &os);
        // 输出流场为Tecplot格式
        void SaveDataTecplot();
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
        double m_res_max;
        double m_res_ave;
        bool m_res_flag;
    };

}