#include"flowsolver.h"
#include "GlobalData.h"
using namespace zaran;

FlowSolver::FlowSolver(int index, string name, FlowSolverPara* para, GridBase* grid, FieldData* fieldData) :FieldSolver(index, name, para, grid, fieldData)
{
}
FlowSolver::~FlowSolver()
{
}
void FlowSolver::InitField()
{

}

void FlowSolver::InitSolver()
{

}

void FlowSolver::Solve()
{

}

void FlowSolver::Postprocess()
{

}

double FlowSolver::ComputeCFL()
{
	return GlobalData::GetDouble("cfl");
}

FlowSolverPara* FlowSolver::GetPara()
{
	return dynamic_cast<FlowSolverPara*>(para_);
}

