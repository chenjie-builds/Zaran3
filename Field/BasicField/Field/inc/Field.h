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
	class Field
	{
	public:
		Field();
		~Field();
	public:
		void SetSolverPara(Ptr<SolverPara>& para) { m_solverPara = para; }
		void SetSolver(Ptr<FieldSolver>& solver) {
			m_solver = solver;
			m_solver->SetFieldData(m_fieldData);
		}
		void SetGrid(Ptr<Grid>& grid) { m_grid = grid; }
		void SetFieldData(Ptr<FieldData>& fieldData) { m_fieldData = fieldData; }
	public:
		Ptr<Grid>& GetGrid() { return m_grid; }
		Ptr<FieldSolver>& GetSolver() { return m_solver; }
		Ptr<FieldData>& GetFieldData() { return m_fieldData; }
		Ptr<SolverPara>& GetSolverPara() { return m_solverPara; }
	private:
		Ptr<Grid> m_grid;
		Ptr<FieldSolver> m_solver;
		Ptr<FieldData> m_fieldData;
		Ptr<SolverPara> m_solverPara;
	};
}