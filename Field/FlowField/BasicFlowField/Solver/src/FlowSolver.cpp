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
void FlowSolver::InitField()
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

void FlowSolver::Postprocess()
{

}

double FlowSolver::ComputeCFL()
{
	return GlobalData::GetDouble("cfl");
}

void FlowSolver::ZeroResidual()
{
	GridPtr& grid = GetGrid();
	int n_equ = GetNumberOfEquations();
	int n_node = grid->GetTotalNodeNum();

	for (int iEqu = 0;iEqu < n_equ;iEqu++)
	{
#pragma omp parallel for
		for (int iNode = 0;iNode < n_node;iNode++)
		{
			GetResidual(iNode, iEqu) = 0.0;
		}
	}

}

FlowSolverPara* FlowSolver::GetPara()
{
	return dynamic_cast<FlowSolverPara*>(para_);
}

