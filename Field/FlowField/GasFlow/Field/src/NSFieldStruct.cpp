#include "NSFieldStruct.h"
namespace zaran
{
    NSFieldStruct::NSFieldStruct(GridBase* grid) :FieldNS(grid, FieldType::NS_Structured)
    {
        m_idx_proxy = nullptr;
        Allocate();
    }
    NSFieldStruct::~NSFieldStruct()
    {
        if (m_idx_proxy != nullptr)
            delete m_idx_proxy;
    }
    GridStruct* NSFieldStruct::GetGrid()
    {
        return static_cast<GridStruct*>(Field::GetGrid());
    }
    NSSolverStruct* NSFieldStruct::GetSolver()
    {
        return static_cast<NSSolverStruct*>(Field::GetSolver());
    }
    void NSFieldStruct::CalcResidual()
    {
        auto grid = GetGrid();
        auto node = grid->GetNode();
        auto para = GetSolverPara();
        int equ_num = para->GetEquNum();
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        double norm_inf = -LARGE_NUMBER;
        double norm_l2 = 0;
        int norm_inf_node = -1;
        double norm_inf_coord[3];
        for (int iEqu = 0; iEqu < equ_num; iEqu++)
        {
            norm_inf = -LARGE_NUMBER;
            norm_l2 = 0;
            auto res = GetDataManager()->GetResidual(iEqu);
            for (int k = ks; k < ke; k++)
            {
                for (int j = js; j < je; j++)
                {
                    for (int i = is; i < ie; i++)
                    {
                        int iNode = m_idx_proxy->GetIdx(i, j, k);
                        if (abs(res[iNode]) > norm_inf)
                        {
                            norm_inf = abs(res[iNode]);
                            norm_inf_node = iNode;
                            for (int iDim = 0; iDim < grid->GetDim(); iDim++)
                            {
                                norm_inf_coord[iDim] = node->GetCoord(i, j, k)[iDim];
                            }
                        }
                        norm_l2 += res[iNode] * res[iNode];
                    }
                }
            }
            norm_l2 = sqrt(norm_l2 / grid->GetTotalNodeNum());
            m_res_info->SetInfNorm(iEqu, norm_inf);
            m_res_info->SetL2Norm(iEqu, norm_l2);
            m_res_info->SetInfNormCoord(iEqu, norm_inf_coord);
            m_res_info->SetInfNormIdx(iEqu, norm_inf_node);
        }
    }
    void NSFieldStruct::Allocate()
    {
        FieldNS::Allocate();
        m_dataManager = new DataManagerNS(GetFieldData(), GetGrid()->GetTotalNodeNum());
        m_solver = new NSSolverStruct(1, "NS_Struct", GetSolverPara(), GetGrid(), GetFieldData(), GetDataManager());
        m_res_info = new ResInfo(GetSolverPara()->GetEquNum());
        m_idx_proxy = new StructIdxProxy(GetGrid());
    }
} // namespace zaran