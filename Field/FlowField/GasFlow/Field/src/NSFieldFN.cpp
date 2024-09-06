#include "NSFieldFN.h"

namespace zaran
{
    NSFieldFNFDM::NSFieldFNFDM(GridBase* grid)
        : FieldNS(grid, FieldType::NS_FlexibleNode)
    {

    }

    NSFieldFNFDM::~NSFieldFNFDM()
    {

    }

    GridFN* NSFieldFNFDM::GetGrid()
    {
        return static_cast<GridFN*>(Field::GetGrid());
    }


    NSSolverFNFDM* NSFieldFNFDM::GetSolver()
    {
        return static_cast<NSSolverFNFDM*>(Field::GetSolver());
    }
    void NSFieldFNFDM::CalcResidual()
    {
        auto grid = GetGrid();
        auto node = grid->GetNode();
        auto para = GetSolverPara();
        int equ_num = para->GetEqNum();
        int inner_node_num = grid->GetInnerNodeNum();
        auto inner_node = grid->GetInnerNode();
        double norm_inf = -LARGE_NUMBER;
        double norm_L2 = 0;
        int norm_inf_node = 0;
        double norm_inf_coord[3];
        for (int iEqu = 0;iEqu < equ_num;iEqu++)
        {
            norm_inf = -LARGE_NUMBER;
            norm_L2 = 0;
            norm_inf_node = -1;
            auto res = GetDataManager()->GetResidual(iEqu);
#pragma omp parallel for reduction(max:norm_inf) reduction(+:norm_L2)
            for (int idx = 0; idx < inner_node_num; idx++)
            {
                int iNode = inner_node[idx];
                if (abs(res[iNode]) > norm_inf)
                {
                    norm_inf = abs(res[iNode]);
                    norm_inf_node = iNode;
                    for (int iDim = 0;iDim < grid->GetDim();iDim++)
                    {
                        norm_inf_coord[iDim] = node->GetCoord(iNode)[iDim];
                    }
                }
                norm_L2 += res[iNode] * res[iNode];
            }
            norm_L2 = sqrt(norm_L2 / inner_node_num);
            m_res_info->SetInfNorm(iEqu, norm_inf);
            m_res_info->SetL2Norm(iEqu, norm_L2);
            m_res_info->SetInfNormCoord(iEqu, norm_inf_coord);
        }
    }

    void NSFieldFNFDM::AllocateSolver()
    {
        if (m_solver != nullptr)
        {
            delete m_solver;
        }
        m_solver = new NSSolverFNFDM(1, "NS_FNFDM", GetSolverPara(), GetGrid(), GetDataManager());
    }
    void NSFieldFNFDM::AllocateDataManager()
    {
        if (m_data_manager != nullptr)
        {
            delete m_data_manager;
        }
        m_data_manager = new DataManagerNS(GetFieldData(), GetGrid()->GetTotalNodeNum());
        m_data_manager->CreateData();
        m_data_manager->RegisterData();
    }
}