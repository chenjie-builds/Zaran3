#include "NSFieldFN.h"

namespace zaran
{
    NSFieldFNFDM::NSFieldFNFDM(GridBase* grid)
        : FieldNS(grid, FieldType::NS_FlexibleNode)
    {
        Allocate();
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
        int equ_num = para->GetEquNum();
        int inner_node_num = grid->GetInnerNodeNum();
        auto inner_node = grid->GetInnerNode();
        double max_res = -LARGE_NUMBER;
        double ave_res = 0;
        int max_res_node = 0;
        double max_res_coord[3];
        for (int iEqu = 0;iEqu < equ_num;iEqu++)
        {
            max_res = -LARGE_NUMBER;
            ave_res = 0;
            max_res_node = -1;
            auto res = GetDataManager()->GetResidual(iEqu);
#pragma omp parallel for reduction(max:max_res) reduction(+:ave_res)
            for (int idx = 0; idx < inner_node_num; idx++)
            {
                int iNode = inner_node[idx];
                if (abs(res[iNode]) > max_res)
                {
                    max_res = abs(res[iNode]);
                    max_res_node = iNode;
                    for (int iDim = 0;iDim < grid->GetDim();iDim++)
                    {
                        max_res_coord[iDim] = node->GetCoord(iNode)[iDim];
                    }
                }
                ave_res += res[iNode] * res[iNode];
            }
            ave_res = sqrt(ave_res / inner_node_num);
            m_res_info->SetMaxRes(iEqu, max_res);
            m_res_info->SetAveRes(iEqu, ave_res);
            m_res_info->SetMaxResCoord(iEqu, max_res_coord);
        }
    }
    void NSFieldFNFDM::Allocate()
    {
        FieldNS::Allocate();
        m_dataManager = new DataManagerNS(GetFieldData(), GetGrid()->GetTotalNodeNum());
        m_solver = new NSSolverFNFDM(1, "NS_FNFDM", GetSolverPara(), GetGrid(), GetFieldData(), GetDataManager());
        m_res_info = new ResInfo(GetSolverPara()->GetEquNum());
    }
}