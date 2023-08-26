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
#include "RiemannSolverPara.h"
#include <memory>
namespace zaran
{
	class RiemannSolver
	{
	public:
		virtual void Solver(std::shared_ptr<RiemannSolverPara>& para) = 0;
	};
}