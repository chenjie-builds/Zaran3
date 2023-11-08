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
#include "controller.h"
#include "GridListFactory.h"
#include "SolverFactory.h"
#include "GridList.h"
namespace zaran
{

	enum class SimulationTask
	{
		SOLVE_FIELD = 0,
		CONVERT_GRID = 1,
		READ_MODEL = 2
	};

	//仿真类，用于控制计算流程
	class Simulation
	{
	public:
		Simulation() :
			globalDataFileName_("zaran.ini"),
			programVersion_("V3.3.1_alpha"),
			lastModifyDate_("2023/11/08"),
			minSupportContrlFileVersion_("V3.3.1_alpha")
		{
			GlobalData::Init();
		}
		~Simulation() {};
		void Start();
	protected:
		//读取控制文件参数
		void ReadGlobalData();
		void InitSimulationTask();
		void ShowInfo();
	private:
		SimulationTask task_;
		//参数文件名
		string globalDataFileName_;
		//网格生成工厂
		Ptr<GridListFactory> gridListFatory_;
		//求解器生成工厂
		Ptr<SolverFactory> solverFactory_;
		//程序版本号, 内部生成
		string programVersion_;
		//程序上次更改日期, 内部生成
		string lastModifyDate_;
		//程序支持的最小控制文件版本号, 如果读取的版本号低于此版本, 将不进行计算
		//作用: 防止程序更新后，需要更多的控制文件信息，老版本控制文件未添加导致计算失败
		string minSupportContrlFileVersion_;
	};

}