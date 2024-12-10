//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Field.h																||
//*	@brief	Control solver to modify field data in grid							||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "BasicType.h"
#include "Solver.h"
#include "GridBase.h"
#include "SolverPara.h"
#include "FieldDataManager.h"
#include <iostream>
#include "FieldSolver.h"
namespace zaran
{
	/// @brief 场的类型
	/// @details 场的类型包括非结构化网格场，结构化网格场，用于告诉外部场的类型，如可视化模块等
	enum class FieldType
	{
		NS_FlexibleNode,
		NS_Structured,
		NS_Zaran,
		Unset,
	};
	// 场类，包含网格，求解器，求解器参数，场数据
	class Field
	{
	public:
		Field(shared_ptr<GridBase> grid = nullptr, FieldType fieldType = FieldType::Unset);
		virtual ~Field();

	public:
		virtual void Allocate();
		void SetIdx(index_type idx) { m_idx = idx; }
		void SetGrid(shared_ptr<GridBase> grid) { m_grid = grid; }
		index_type GetIdx() const { return m_idx; }
		shared_ptr<GridBase> GetGrid() const { return m_grid; }
		shared_ptr<FieldSolver> GetSolver() const { return m_solver; }
		shared_ptr<FieldData> GetData() const { return m_data; }
		shared_ptr<SolverParam> GetSolverPara() const { return m_solver_para; }
		shared_ptr<DataManager> GetDataManager() const { return m_data_manager; }
		FieldType GetFieldType() const { return m_fieldType; }

	protected:
		virtual void AllocateFieldData();
		virtual void AllocateSolver() = 0;
		virtual void AllocateSolverPara() = 0;
		virtual void AllocateDataManager() = 0;

	protected:
		/// @brief 本场在全场数组中的索引
		index_type m_idx;
		/// @brief 场的网格，可能是结构化网格，也可能是非结构化网格
		shared_ptr<GridBase> m_grid;
		/// @brief 场的求解器，用于求解场
		shared_ptr<FieldSolver> m_solver;
		/// @brief 场的数据，包含场的值，梯度等
		shared_ptr<FieldData> m_data;
		/// @brief 场的求解器参数，用于设置求解器的参数
		shared_ptr<SolverParam> m_solver_para;
		/// @brief 场的数据管理器，用于管理场的数据
		shared_ptr<DataManager> m_data_manager;
		/// @brief 场的类型
		FieldType m_fieldType;
	};
}