#include "FlowSolver.h"
#include "GlobalData.h"
using namespace zaran;

FlowFieldSolver::FlowFieldSolver(int index, string name, FlowSolverPara *para, GridBase *grid, DataManager *data_manager) : FieldSolver(index, name, para, grid, data_manager)
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

FlowSolverPara *FlowFieldSolver::GetPara()
{
	return dynamic_cast<FlowSolverPara *>(m_para);
}
