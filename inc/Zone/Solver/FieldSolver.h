//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	FieldSolver.h														||
//*	@brief	流场求解器															||
//*	@author	Chen Jie.															||
//==============================================================================||

#pragma once
#include"Solver.h"
#include"FieldDataManager.h"
namespace zaran
{
	enum class FieldSolverType
	{
		NS_Struct,
		NS_FNFDM,
		NS_ZaRan,
	};
	/// @brief 场求解器基类
	class FieldSolver :public Solver
	{
	public:
		FieldSolver(Id index, string name, std::shared_ptr<SolverParam> para, std::shared_ptr < GridBase> grid, std::shared_ptr < DataManager>data_manager);
		virtual ~FieldSolver();
	public:
		// 初始化流场数据
		virtual void InitField() = 0;
		// 边界处理
		virtual void BoundaryCondition() = 0;
		// 更新场数据, 根据迭代中间结果更新场数据
		virtual void UpdateField() = 0;
		// 备份场数据
		virtual void BackupField(std::string& back_folder) = 0;
		std::shared_ptr<DataManager>& GetDataManager() { return m_data_manager; }
	private:
		// 场数据
		std::shared_ptr<DataManager> m_data_manager;
	};
}
