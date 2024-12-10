#include "FieldSolver.h"
#include "Log.h"
namespace zaran
{
	FieldSolver::FieldSolver(index_type index, string name, shared_ptr<SolverParam> para, shared_ptr < GridBase> grid, shared_ptr < DataManager>data_manager) :
		Solver(index, name, para, grid), m_data_manager(data_manager)
	{
	}

	FieldSolver::~FieldSolver()
	{
	}

}