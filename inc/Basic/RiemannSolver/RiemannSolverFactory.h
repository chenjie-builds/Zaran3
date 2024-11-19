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
#include "RiemannSolver.h"
namespace zaran
{

	class RiemannSolverBuilder
	{
	public:
		 RiemannSolver* Create(RiemannSolverType type);
	};
}