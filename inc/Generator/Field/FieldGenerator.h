/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file FieldGenerator.h
 * \brief Field generator class, used to generate field manager.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */

#pragma once
#include "GlobalField.h"
namespace zaran
{
	class FieldGenerator
	{
	public:
		FieldGenerator(GridType grid_type, FieldSolverType solver_type, Dimension dim)
			: m_grid_type(grid_type), m_solver_type(solver_type), m_dim(dim) {
		};
		virtual shared_ptr<FieldManager> Create();

	private:
		void CreateGrid();
		void CreateField();
		void CreateSolver();
		shared_ptr<FieldManager> CreateFieldZaran();

	private:
		GridType m_grid_type;
		FieldSolverType m_solver_type;
		Dimension m_dim;
	};
} // namespace zaran