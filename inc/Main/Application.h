//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Simulation.h														||
//*	@brief	模拟类, 根据任务的不同执行不同任务										||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "GlobalData.h"
#include "FieldSimulation.h"
#include "GridFactory.h"
#include "GridList.h"
namespace zaran
{

    enum class TaskType
    {
		/// @brief 用于求解场问题
        SOLVE_FIELD = 0,
		/// @brief 网格转换，将网格转换为程序所需的格式
        CONVERT_GRID = 1,
		/// @brief 读取模型
        READ_MODEL = 2
    };

    //应用类，用于控制整个程序的功能
    class Application
    {
    public:
        Application() :
            m_control_file("zaran.ini"),
            m_min_ctrol_file_version("V3.3.7"),
            m_grid_builder(nullptr),
            m_task(TaskType::SOLVE_FIELD) // Ensure m_task is initialized
        {
            GlobalData::Init();
        }
        ~Application() {};
        void Run();
    protected:
        //读取控制文件参数
        void ReadGlobalData();
        void InitTask();
        void ShowInfo();
        void SolveField();
        void ConvertGrid();
        void ReadModel();
    private:
        TaskType m_task;
        //参数文件名
        string m_control_file;
        //网格生成工厂
        Ptr<GridBuilder> m_grid_builder;
        //程序支持的最小控制文件版本号, 如果读取的版本号低于此版本, 将不进行计算
        //作用: 防止程序更新后，需要更多的控制文件信息，老版本控制文件未添加导致计算失败
        string m_min_ctrol_file_version;
    };

}