//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	SolverFactory.h														||
//*	@brief	求解器工厂类, 生成求解器												||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "SolverVec.h"
#include "solver.h"
#include "flowsolver.h"
#include "NSSolver.h"
#include "Solver_NS_3D.h"
#include"Solver_NS_2D.h"
#include "Solver_NS_3D_Struct.h"
#include "Solver_NS_2D_Struct.h"
#include"Solver_NS_3D_Zaran.h"
namespace zaran
{
	class SolverFactory
	{
	public:
		void Create(Grid* grid, Solver*& solver, FieldSolverType& solverType);
	};
}