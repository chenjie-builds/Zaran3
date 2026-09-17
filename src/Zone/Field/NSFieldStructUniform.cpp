/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file NSFieldStructUniform.cpp
 * \brief NS field on a uniform Cartesian structured grid.
 */

#include "NSFieldStructUniform.h"

namespace zaran
{
	NSFieldStructUniform::NSFieldStructUniform(shared_ptr<GridStruct> grid)
		: NSFieldStruct(grid)
	{
	}

	NSFieldStructUniform::~NSFieldStructUniform()
	{
	}

	shared_ptr<EulerSolverStructUniform> NSFieldStructUniform::GetSolver()
	{
		return std::static_pointer_cast<EulerSolverStructUniform>(Field::GetSolver());
	}

	shared_ptr<FlowSolverParamUniform> NSFieldStructUniform::GetSolverPara()
	{
		return std::static_pointer_cast<FlowSolverParamUniform>(Field::GetSolverPara());
	}

	void NSFieldStructUniform::AllocateSolver()
	{
		m_solver = make_shared<EulerSolverStructUniform>(
			GetIdx(), "Euler_Struct_Uniform", GetSolverPara(), GetGrid(), GetDataManager());
	}

	void NSFieldStructUniform::AllocateSolverPara()
	{
		m_solver_para = make_shared<FlowSolverParamUniform>();
		// GetSolverPara() 此处调用的是本类的隐藏版本，返回 FlowSolverParamUniform
		std::static_pointer_cast<FlowSolverParamUniform>(m_solver_para)->Init();
	}
}
