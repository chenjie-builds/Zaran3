//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	VanLeer.h															||
//*	@brief	van Leer 黎曼求解器													||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "RiemannSolver.h"
namespace zaran
{
	class Vanleer :public RiemannSolver
	{
	public:
		void Solver(RiemannSolverPara& para)override;
	};
}