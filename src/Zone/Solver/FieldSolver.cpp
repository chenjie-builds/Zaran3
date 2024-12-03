#include "FieldSolver.h"
#include "Log.h"
namespace zaran
{
    FieldSolver::FieldSolver(int index, string name, SolverParam* para, GridBase* grid,DataManager* data_manager) :
        Solver(index, name, para, grid), m_data_manager(data_manager)
    {
    }

    FieldSolver::~FieldSolver()
    {
    }

}