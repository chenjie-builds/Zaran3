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
#include "Grid.h"
#include "Visual.h"
#include "SolverPara.h"
#include <iostream>
namespace zaran
{
	//场类，包含网格，求解器，求解器参数，场数据
	class Field
	{
	public:
		Field();
		~Field();
	public:
		void SetSolverPara(SolverPara* para) { m_solverPara = para; }
		void SetSolver(FieldSolver* solver) 
		{
			m_solver = solver;
			m_solver->SetFieldData(m_fieldData);
		}
		void SetGrid(Grid* grid) { m_grid = grid; }
		void SetFieldData(FieldData* fieldData) { m_fieldData = fieldData; }
	public:
		Grid* GetGrid() { return m_grid; }
		FieldSolver* GetSolver() { return m_solver; }
		FieldData* GetFieldData() { return m_fieldData; }
		SolverPara* GetSolverPara() { return m_solverPara; }
	private:
		Grid* m_grid;
		FieldSolver* m_solver;
		FieldData* m_fieldData;
		SolverPara* m_solverPara;
	};
}