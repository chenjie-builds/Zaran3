#include "NSSolverStruct.h"
#include "GlobalData.h"
#include "log.h"
using namespace zaran;
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

Ptr<StructGrid> NSSolverStruct::NSSolverStruct::GetGrid()
{
	return std::static_pointer_cast<StructGrid> (Solver::GetGrid());
}


void NSSolverStruct::ComputeTimeStepLocal()
{

}

