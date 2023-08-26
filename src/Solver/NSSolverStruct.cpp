#include "NSSolverStruct.h"
#include "GlobalData.h"
#include "log.h"
namespace zaran
{

	void NSSolverStruct::InitField()
	{

	}

	void NSSolverStruct::InitSolver()
	{

	}


	void NSSolverStruct::ComputeTimeStep()
	{

	}

	void NSSolverStruct::TimeAdvance()
	{

	}

	void NSSolverStruct::BoundaryCondition()
	{

	}

	void NSSolverStruct::ComputePrimtiveGradient()
	{

	}

	NSSolverStruct::StructGridPtr NSSolverStruct::GetGrid()
	{
		return std::static_pointer_cast<StructGrid> (Solver::GetGrid());
	}


	void NSSolverStruct::ComputeTimeStepLocal()
	{

	}

}
