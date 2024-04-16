#include "SolverFactory.h"
#include "GlobalData.h"
#include "GridList.h"

using namespace zaran;
void SolverFactory::Create(Ptr<GridList>& gridList, Ptr<SolverVec>& solverVecPtr, FieldSolverType& solverType)
{
	if (solverVecPtr.get() == nullptr)
	{
		solverVecPtr = std::make_shared<SolverVec>();
	}
	for (int iSolver = 0; iSolver < gridList->GetGridNumber(); ++iSolver)
	{
		Ptr<Solver> newSolver = std::make_shared<Solver_NS_2D>();
		if (solverType == FieldSolverType::NS_2D)
			newSolver = std::make_shared<Solver_NS_2D>();
		else if (solverType == FieldSolverType::NS_3D)
			newSolver = std::make_shared<Solver_NS_3D>();
		else if (solverType == FieldSolverType::NS_2D_Struct)
			newSolver = std::make_shared<Solver_NS_2D_Struct>();
		else if (solverType == FieldSolverType::NS_3D_Struct)
			newSolver = std::make_shared<Solver_NS_3D_Struct>();
		else if (solverType == FieldSolverType::NS_ZaRan_3D)
			newSolver = std::make_shared<Solver_NS_3D_Zaran>();
		else
		{
			Log::warn("Unsupported Solver Type! Please Check!");
			system("pause");
		}
		newSolver->SetGridList(gridList);
		newSolver->SetGridIndex(iSolver);
		newSolver->SetName("test");
		FlowSolverParaPtr para = std::make_shared<FlowSolverPara>();
		para->Init();
		newSolver->SetPara(para);
		solverVecPtr->AddSolver(newSolver);
	}
}
