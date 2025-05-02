/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file RiemannSolver.h
 * \brief Riemann solver base class.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */

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