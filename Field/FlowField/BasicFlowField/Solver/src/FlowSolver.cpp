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
	for (int i = 0;i < m_Residual.size();i++)
	{
#pragma omp parallel for
		for (int j = 0;j < m_Residual[i]->size();j++)
		{
			(*m_Residual[i])[j] = 0.0;
		}
	}

}

FlowSolverParaPtr FlowSolver::GetPara()
{
	return std::static_pointer_cast<FlowSolverPara>(para_);
}

