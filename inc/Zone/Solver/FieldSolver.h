/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file FieldSolver.h
 * \brief FieldSolver class, used to solve field data.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
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
		DEM,
	};
	/// @brief 场求解器基类
	class FieldSolver :public Solver
	{
	public:
		FieldSolver(index_type index, string name, shared_ptr<SolverParam> para, shared_ptr < GridBase> grid, shared_ptr < DataManager>data_manager);
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
		shared_ptr<DataManager>& GetDataManager() { return m_data_manager; }
	private:
		// 场数据
		shared_ptr<DataManager> m_data_manager;
	};
}
