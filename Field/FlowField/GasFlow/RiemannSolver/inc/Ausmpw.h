//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Ausmpw.h															||
//*	@brief	AUSM+ ÀèÂüÇó½âÆ÷														||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "RiemannSolver.h"
namespace zaran
{
	class Ausmpw :public RiemannSolver
	{
	public:
		void Solver(RiemannSolverPara& para)override;
	};
}