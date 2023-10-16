#include"flowsolver.h"
#include "GlobalData.h"
using namespace zaran;
void FlowSolver::InitField()
{

}

void FlowSolver::InitSolver()
{
	para_->Init();
	GridPtr& grid = GetGrid();
	int nodeNum = grid->GetTotalNodeNum();
	CreateFieldData();
	RegisterFieldData();
}

void FlowSolver::Solve()
{

}

void FlowSolver::Post()
{

}

double FlowSolver::ComputeCFL()
{
	return GlobalData::GetDouble("cfl");
}

void FlowSolver::ZeroResidual()
{
	for (auto& res : m_Residual)
	{
		for (auto& resData : *res)
		{
			resData = 0.0;
		}
	}
}

FlowSolverParaPtr FlowSolver::GetPara()
{
	return std::static_pointer_cast<FlowSolverPara>(para_);
}

