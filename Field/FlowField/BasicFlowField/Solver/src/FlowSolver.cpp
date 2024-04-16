#include"flowsolver.h"
#include "GlobalData.h"
using namespace zaran;
zaran::FlowSolver::FlowSolver()
{
	para_ = new FlowSolverPara;
}
zaran::FlowSolver::~FlowSolver()
{
	delete[] para_;
}
void FlowSolver::InitData()
{

}

void FlowSolver::InitSolver()
{
	para_->Init();
	CreateData();
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
	GridPtr& grid = GetGrid();
	int n_equation = GetNumberOfEquations();
	int n_node = grid->GetTotalNodeNum();

	for (int i = 0;i < n_equation;i++)
	{
#pragma omp parallel for
		for (int j = 0;j < n_node;j++)
		{
			m_residual[i][j] = 0.0;
		}
	}

}

FlowSolverPara* FlowSolver::GetPara()
{
	return dynamic_cast<FlowSolverPara*>(para_);
}

