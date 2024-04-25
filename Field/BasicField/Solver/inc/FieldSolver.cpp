#include "FieldSolver.h"
#include "Log.h"
zaran::FieldSolver::FieldSolver()
{
    m_field_data=nullptr;

}

zaran::FieldSolver::~FieldSolver()
{
    Log::info("FieldSolverÎö¹¹º¯Êý");
}

void zaran::FieldSolver::SetFieldData(FieldData* fieldData)
{
    m_field_data=fieldData;
}

