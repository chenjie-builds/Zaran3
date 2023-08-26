#include "SolverFactory.h"
#include "GlobalData.h"
#include "GridList.h"

void zaran::SolverFactory::Create(std::shared_ptr<GridList>& gridList, std::shared_ptr<SolverVec>& solverVecPtr)
{
	if (solverVecPtr.get() == nullptr)
	{
		solverVecPtr = std::make_shared<SolverVec>();
	}
	for (int iSolver = 0; iSolver < gridList->GetGridNumber(); ++iSolver)
	{
		std::shared_ptr<Solver> newSolver = std::make_shared<NSSolver>();
		newSolver->SetGridList(gridList);
		newSolver->SetGridIndex(iSolver);
		newSolver->SetName("test");
		FlowSolverParaPtr para = std::make_shared<FlowSolverPara>();
		para->Init();
		newSolver->SetPara(para);
		solverVecPtr->AddSolver(newSolver);
	}
}
