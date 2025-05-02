/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file RiemannSolverFactory.h
 * \brief Riemann solver Factory class, used to create Riemann solver.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */

#pragma once
#include "RiemannSolver.h"
#include <memory>
namespace zaran
{

	class RiemannSolverBuilder
	{
	public:
		 RiemannSolver* Create(RiemannSolverType type);
		 std::unique_ptr<RiemannSolver> CreateUnique(RiemannSolverType type);
	};
}