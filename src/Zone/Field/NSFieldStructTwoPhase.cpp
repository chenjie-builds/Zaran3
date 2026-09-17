/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file NSFieldStructTwoPhase.cpp
 * \brief Two-phase NS field on a uniform Cartesian structured grid.
 */

#include "NSFieldStructTwoPhase.h"
#include "Log.h"
#include "ZaranError.h"

namespace zaran
{
	NSFieldStructTwoPhase::NSFieldStructTwoPhase(shared_ptr<GridStruct> grid)
		: NSFieldStruct(grid)
	{
	}

	NSFieldStructTwoPhase::~NSFieldStructTwoPhase()
	{
	}

	shared_ptr<EulerTwoPhaseStructUniform> NSFieldStructTwoPhase::GetTwoPhaseSolver()
	{
		return std::dynamic_pointer_cast<EulerTwoPhaseStructUniform>(Field::GetSolver());
	}

	shared_ptr<FlowSolverParamTwoPhase> NSFieldStructTwoPhase::GetTwoPhaseSolverPara()
	{
		return std::dynamic_pointer_cast<FlowSolverParamTwoPhase>(Field::GetSolverPara());
	}

	shared_ptr<DataManagerNSTwoPhase> NSFieldStructTwoPhase::GetTwoPhaseDataManager()
	{
		return std::dynamic_pointer_cast<DataManagerNSTwoPhase>(Field::GetDataManager());
	}

	void NSFieldStructTwoPhase::AllocateSolverPara()
	{
		m_solver_para = make_shared<FlowSolverParamTwoPhase>();
		std::static_pointer_cast<FlowSolverParamTwoPhase>(m_solver_para)->Init();
	}

	void NSFieldStructTwoPhase::AllocateDataManager()
	{
		auto grid_struct = std::dynamic_pointer_cast<GridStruct>(GetGrid());
		if (!grid_struct)
		{
			throw ZaranError("NSFieldStructTwoPhase requires a structured grid");
		}
		auto data_manager = make_shared<DataManagerNSTwoPhase>(
			GetData(), static_cast<int>(grid_struct->GetNi()),
			static_cast<int>(grid_struct->GetNj()), static_cast<int>(grid_struct->GetNk()));
		data_manager->CreateData();
		data_manager->RegisterData();
		m_data_manager = data_manager;
	}

	void NSFieldStructTwoPhase::AllocateSolver()
	{
		auto para = std::dynamic_pointer_cast<FlowSolverParamTwoPhase>(m_solver_para);
		auto data_manager = std::dynamic_pointer_cast<DataManagerNS>(m_data_manager);
		if (!para || !data_manager)
		{
			throw ZaranError("NSFieldStructTwoPhase: inconsistent solver parameter / data manager");
		}
		m_solver = make_shared<EulerTwoPhaseStructUniform>(
			GetIdx(), "Euler_TwoPhase_Struct_Uniform", para, GetGrid(), data_manager);
	}
}
