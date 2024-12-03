#include "NSField.h"
namespace zaran
{
    FieldNS::FieldNS(GridBase* grid, FieldType fieldType) :Field(grid, fieldType)
    {
        m_res_info = nullptr;
    }
    FieldNS::~FieldNS()
    {
        if (m_res_info != nullptr)
            delete m_res_info;
    }
    FlowSolverParam* FieldNS::GetSolverPara()
    {
        return static_cast<FlowSolverParam*>(Field::GetSolverPara());
    }
    DataManagerNS* FieldNS::GetDataManager()
    {
        return static_cast<DataManagerNS*>(Field::GetDataManager());
    }
    void FieldNS::AllocateSolverPara()
    {
        if (m_solver_para != nullptr)
        {
            delete m_solver_para;
        }
        m_solver_para = new FlowSolverParam();
        GetSolverPara()->Init();
    }
    void FieldNS::Allocate()
    {
        Field::Allocate();
        AllocateResInfo();
    }

    void FieldNS::AllocateResInfo()
    {
        if (m_res_info != nullptr)
        {
            delete m_res_info;
        }
        m_res_info = new ResInfo(GetSolverPara()->GetEqNum());
    }

}