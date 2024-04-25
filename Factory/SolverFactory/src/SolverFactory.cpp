#include "SolverFactory.h"
#include "GlobalData.h"
#include "GridList.h"

using namespace zaran;
void SolverFactory::Create(Grid* grid, Solver*& solver, FieldSolverType& solverType)
{
	if (solver != nullptr)
	{
		delete[] solver;
		solver = nullptr;
	}
	if (solverType == FieldSolverType::NS_2D)
		solver = new Solver_NS_2D();
	else if (solverType == FieldSolverType::NS_3D)
		solver = new Solver_NS_3D();
	else if (solverType == FieldSolverType::NS_2D_Struct)
		solver = new Solver_NS_2D_Struct();
	else if (solverType == FieldSolverType::NS_3D_Struct)
		solver = new Solver_NS_3D_Struct();
	else if (solverType == FieldSolverType::NS_ZaRan_3D)
		solver = new Solver_NS_3D_Zaran();
	else
	{
		Log::warn("Unsupported Solver Type! Please Check!");
		system("pause");
	}
	solver->SetGrid(grid);
	solver->SetName("test");
}

