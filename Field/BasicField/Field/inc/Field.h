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
#include"BasicType.h"
#include "Solver.h"
#include "GridBase.h"
#include "SolverPara.h"
#include "FieldDataManager.h"
#include <iostream>
#include"FieldSolver.h"
namespace zaran
{
	/// @brief 场的类型
	/// @details 场的类型包括非结构化网格场，结构化网格场，用于告诉外部场的类型，如可视化模块等
	enum class FieldType
	{
		NS_FlexibleNode,
		NS_Structured,
		Unset,
	};
	//场类，包含网格，求解器，求解器参数，场数据
	class Field
	{
	public:
		Field(GridBase* grid, FieldType fieldType = FieldType::Unset);
		virtual~Field();
	public:
		void SetIdx(int idx) { m_idx = idx; }
		int GetIdx() { return m_idx; }
		virtual GridBase* GetGrid() { return m_grid; }
		virtual FieldSolver* GetSolver() { return m_solver; }
		virtual FieldData* GetFieldData() { return m_fieldData; }
		virtual SolverPara* GetSolverPara() { return m_solver_para; }
		virtual DataManager* GetDataManager() { return m_data_manager; }
		FieldType GetFieldType() { return m_fieldType; }
	protected:
		virtual void Allocate();
	protected:
		/// @brief 本场在全场数组中的索引
		int m_idx;
		/// @brief 场的网格，可能是结构化网格，也可能是非结构化网格
		GridBase* m_grid;
		/// @brief 场的求解器，用于求解场
		FieldSolver* m_solver;
		/// @brief 场的数据，包含场的值，梯度等
		FieldData* m_fieldData;
		/// @brief 场的求解器参数，用于设置求解器的参数
		SolverPara* m_solver_para;
		/// @brief 场的数据管理器，用于管理场的数据
		DataManager* m_data_manager;
		/// @brief 场的类型
		FieldType m_fieldType;
	};
}