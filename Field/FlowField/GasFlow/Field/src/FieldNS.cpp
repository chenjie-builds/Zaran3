#include "FieldNS.h"

namespace zaran
{
    FieldNS_FNFDM::FieldNS_FNFDM(GridBase* grid)
        : Field(grid, FieldType::NS_FlexibleNode)
    {
        Allocate();
    }

    FieldNS_FNFDM::~FieldNS_FNFDM()
    {
        delete m_res_analyze;
    }

    GridFN* FieldNS_FNFDM::GetGrid()
    {
        return static_cast<GridFN*>(Field::GetGrid());
    }

    NSSolverFNFDM* FieldNS_FNFDM::GetSolver()
    {
        return static_cast<NSSolverFNFDM*>(Field::GetSolver());
    }
    FlowSolverPara* FieldNS_FNFDM::GetSolverPara()
    {
        return static_cast<FlowSolverPara*>(Field::GetSolverPara());
    }
    DataManagerNS* FieldNS_FNFDM::GetDataManager()
    {
        return static_cast<DataManagerNS*>(Field::GetDataManager());
    }
    void FieldNS_FNFDM::Allocate()
    {
        Field::Allocate();
        m_solver_para = new FlowSolverPara();
        m_solver_para->Init();
        m_dataManager = new DataManagerNS(GetFieldData(), GetGrid()->GetTotalNodeNum());
        m_solver = new NSSolverFNFDM(1, "NS_FNFDM", GetSolverPara(), GetGrid(), GetFieldData(), GetDataManager());
        double** res = new double* [GetSolverPara()->GetEquNum()];
        for (int i = 0; i < GetSolverPara()->GetEquNum(); ++i)
        {
            res[i] = GetDataManager()->GetResidual(i);
        }
        m_res_analyze = new ResAnalyzerFN(GetGrid(), res, GetSolverPara()->GetEquNum());
    }
}