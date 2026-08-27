/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file Application.h
 * \brief Application class, used to control the simulation.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "GlobalData.h"
#include "FieldSimulation.h"
#include "GridGenerator.h"
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
        Application(std::string work_dir) :
			m_work_dir(work_dir),
			m_control_file(work_dir + "/zaran.ini"),
            m_min_ctrol_file_version("V3.3.7"),
            m_grid_builder(nullptr),
            m_task(TaskType::SOLVE_FIELD) // Ensure m_task is initialized
        {
            GlobalData::Init();
			GlobalData::Update("work_dir", work_dir);
        }
        ~Application() {};
        void Run();
    protected:
        //读取控制文件参数
        void ReadGlobalData();
        void InitTask();
        void ShowInfo();
        bool IsVersionTooOld(const std::string& version) const;
        void SolveField() const;
        void ConvertGrid();
        void ReadModel();
    private:
		string m_work_dir;
        TaskType m_task;
        //参数文件名
        string m_control_file;
        //网格生成工厂
        shared_ptr<GridGenerator> m_grid_builder;
        //程序支持的最小控制文件版本号, 如果读取的版本号低于此版本, 将不进行计算
        //作用: 防止程序更新后，需要更多的控制文件信息，老版本控制文件未添加导致计算失败
        string m_min_ctrol_file_version;
    };

}