//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	StegerWarming.h														||
//*	@brief	StegerWarming ÀèÂüÇó½âÆ÷				||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "RiemannSolver.h"
namespace zaran
{
	class StegerWarming :public RiemannSolver
	{
	public:
		void Solver(Ptr<RiemannSolverPara>& para)override;
	};
}