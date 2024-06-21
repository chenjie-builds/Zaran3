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
		Field(GridBase* grid,FieldType fieldType);
		virtual~Field();
	public:
		void SetIdx(int idx) { m_idx = idx; }
		int GetIdx() { return m_idx; }
		virtual GridBase* GetGrid() { return m_grid; }
		virtual FieldSolver* GetSolver() { return m_solver; }
		virtual FieldData* GetFieldData() { return m_fieldData; }
		virtual SolverPara* GetSolverPara() { return m_solver_para; }
		virtual DataManager* GetDataManager() { return m_dataManager; }
		FieldType GetFieldType() { return m_fieldType; }
		protected:
		virtual void Allocate();
	protected:
		/// @brief 本场在场数组中的索引
		int m_idx;
		GridBase* m_grid;
		FieldSolver* m_solver;
		FieldData* m_fieldData;
		SolverPara* m_solver_para;
		DataManager* m_dataManager;
		FieldType m_fieldType;
	};
}