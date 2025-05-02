/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file HLLC.h
 * \brief HLLC Riemann solver class.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "RiemannSolver.h"
namespace zaran
{
	class HLLC :public RiemannSolver
	{
	public:
		void Solver(RiemannSolverPara& para)override;
	};
}