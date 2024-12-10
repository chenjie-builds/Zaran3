#include "FlowSolver.h"
#include "GlobalData.h"
using namespace zaran;

FlowFieldSolver::FlowFieldSolver(index_type index, string name, shared_ptr<FlowSolverParam> para, shared_ptr < GridBase> grid, shared_ptr < DataManager>data_manager) : FieldSolver(index, name, para, grid, data_manager)
{
	m_para = para;
	
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

FlowSolverParam* FlowFieldSolver::GetPara()
{
	return m_para.get();
}