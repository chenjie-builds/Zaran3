//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	HLLC.h																||
//*	@brief	HLLC ÀèÂüÇó½âÆ÷														||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "RiemannSolver.h"
namespace zaran
{
	class HLLC :public RiemannSolver
	{
	public:
		void Solver(Ptr<RiemannSolverPara>& para)override;
	};
}