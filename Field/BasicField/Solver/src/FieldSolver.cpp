#include "FieldSolver.h"
#include "Log.h"
namespace zaran
{
    FieldSolver::FieldSolver(int index, string name, SolverPara* para, GridBase* grid, FieldData* fieldData) :
        Solver(index, name, para, grid)
    {
        m_field_data = fieldData;
        m_equ_num = 0;
    }

    FieldSolver::~FieldSolver()
    {
    }

    void FieldSolver::SetFieldData(FieldData* fieldData)
    {
        m_field_data = fieldData;
    }
}