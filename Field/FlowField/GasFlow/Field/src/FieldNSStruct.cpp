#include "FieldNSStruct.h"
namespace zaran
{
    FieldNS_Struct::FieldNS_Struct(GridBase* grid):Field(grid,FieldType::NS_Structured)
    {
        Allocate();
    }
    FieldNS_Struct::~FieldNS_Struct()
    {
    }
    GridStruct* FieldNS_Struct::GetGrid()
    {
        return static_cast<GridStruct*>(Field::GetGrid());
    }
    NSSolverStruct* FieldNS_Struct::GetSolver()
    {
        return static_cast<NSSolverStruct*>(Field::GetSolver());
    }
    FlowSolverPara* FieldNS_Struct::GetSolverPara()
    {
        return static_cast<FlowSolverPara*>(Field::GetSolverPara());
    }
    DataManagerNS* FieldNS_Struct::GetDataManager()
    {
        return static_cast<DataManagerNS*>(Field::GetDataManager());
    }
    void FieldNS_Struct::Allocate()
    {
        Field::Allocate();
        m_solver_para = new FlowSolverPara();
        m_solver_para->Init();
        m_dataManager = new DataManagerNS(GetFieldData(), GetGrid()->GetTotalNodeNum());
        m_solver = new NSSolverStruct(1, "NS_Struct", GetSolverPara(), GetGrid(), GetFieldData(), GetDataManager());
    }
} // namespace zaran