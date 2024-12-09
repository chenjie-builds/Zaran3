#include "FieldSolver.h"
#include "Log.h"
namespace zaran
{
	FieldSolver::FieldSolver(Id index, string name, std::shared_ptr<SolverParam> para, std::shared_ptr < GridBase> grid, std::shared_ptr < DataManager>data_manager) :
		Solver(index, name, para, grid), m_data_manager(data_manager)
	{
	}

	FieldSolver::~FieldSolver()
	{
	}

}