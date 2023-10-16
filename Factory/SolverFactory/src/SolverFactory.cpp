#include "SolverFactory.h"
#include "GlobalData.h"
#include "GridList.h"

using namespace zaran;
void SolverFactory::Create(Ptr<GridList>& gridList, Ptr<SolverVec>& solverVecPtr)
{
	if (solverVecPtr.get() == nullptr)
	{
		solverVecPtr = std::make_shared<SolverVec>();
	}
	for (int iSolver = 0; iSolver < gridList->GetGridNumber(); ++iSolver)
	{
		Ptr<Solver> newSolver = std::make_shared<Solver_NS_2D>();
		newSolver->SetGridList(gridList);
		newSolver->SetGridIndex(iSolver);
		newSolver->SetName("test");
		FlowSolverParaPtr para = std::make_shared<FlowSolverPara>();
		para->Init();
		newSolver->SetPara(para);
		solverVecPtr->AddSolver(newSolver);
	}
}
