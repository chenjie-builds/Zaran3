#include "NSField.h"
namespace zaran
{
    FieldNS::FieldNS(GridBase* grid, FieldType fieldType) :Field(grid, fieldType)
    {
        m_res_info = nullptr;
        Allocate();
    }
    FieldNS::~FieldNS()
    {
        if(m_res_info != nullptr)
        delete m_res_info;
    }
    FlowSolverPara* FieldNS::GetSolverPara()
    {
        return static_cast<FlowSolverPara*>(Field::GetSolverPara());
    }
    DataManagerNS* FieldNS::GetDataManager()
    {
        return static_cast<DataManagerNS*>(Field::GetDataManager());
    }
    void FieldNS::Allocate()
    {
        Field::Allocate();
        m_solver_para = new FlowSolverPara();
        m_solver_para->Init();
    }
}