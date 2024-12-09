#include "FlowSolver.h"
#include "GlobalData.h"
using namespace zaran;

FlowFieldSolver::FlowFieldSolver(Id index, string name, std::shared_ptr<FlowSolverParam> para, std::shared_ptr < GridBase> grid, std::shared_ptr < DataManager>data_manager) : FieldSolver(index, name, para, grid, data_manager)
{
}
FlowFieldSolver::~FlowFieldSolver()
{
}
void FlowFieldSolver::InitField()
{
}

void FlowFieldSolver::InitSolver()
{
}

void FlowFieldSolver::Solve()
{
}

void FlowFieldSolver::Postprocess()
{
}

double FlowFieldSolver::ComputeCFL()
{
	return GlobalData::GetDouble("cfl");
}

std::shared_ptr<FlowSolverParam> FlowFieldSolver::GetPara()
{
	return std::static_pointer_cast<FlowSolverParam>(FieldSolver::GetPara());
}
