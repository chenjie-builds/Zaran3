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
#include "NSSolverStruct.h"
#include "NSSolverUnstruct.h"
namespace zaran
{
	class SolverFactory
	{
	public:
		void Create(Ptr<GridList>& gridList, Ptr<SolverVec>& solverVecPtr);
	};
}