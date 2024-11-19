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
	//数值通量计算方法
	enum class RiemannSolverType
	{
		VanLeer,
		HLLC,
		Roe,
		StegerWarming,
		Ausmpw
	};
	class RiemannSolver
	{
	public:
		virtual void Solver(RiemannSolverPara& para) = 0;
		virtual ~RiemannSolver() {}
	};
}