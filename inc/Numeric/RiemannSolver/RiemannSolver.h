//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	.h															||
//*	@brief						||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include"BasicType.h"
#include "RiemannSolverPara.h"
namespace zaran
{
	class RiemannSolver
	{
	public:
		virtual void Solver(Ptr<RiemannSolverPara>& para) = 0;
	};
}