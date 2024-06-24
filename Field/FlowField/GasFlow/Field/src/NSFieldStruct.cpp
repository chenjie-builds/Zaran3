#include "NSFieldStruct.h"
namespace zaran
{
    NSFieldStruct::NSFieldStruct(GridBase* grid) :FieldNS(grid, FieldType::NS_Structured)
    {
        Allocate();
    }
    NSFieldStruct::~NSFieldStruct()
    {
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
        for (int iEqu = 0; iEqu < equ_num; iEqu++)
        {
            double max_res = -LARGE_NUMBER;
            double ave_res = 0;
            int max_res_node = -1;
            double max_res_coord[3];
            auto res = GetDataManager()->GetResidual(iEqu);
            for (int k = ks; k < ke; k++)
            {
                for (int j = js; j < je; j++)
                {
                    for (int i = is; i < ie; i++)
                    {
                        int iNode = m_idx_proxy->GetIdx(i, j, k);
                        if (abs(res[iNode]) > max_res)
                        {
                            max_res = abs(res[iNode]);
                            max_res_node = iNode;
                            for (int iDim = 0; iDim < grid->GetDim(); iDim++)
                            {
                                max_res_coord[iDim] = node->GetCoord(i, j, k)[iDim];
                            }
                        }
                        ave_res += res[iNode] * res[iNode];
                    }
                }
            }
            ave_res = sqrt(ave_res / grid->GetTotalNodeNum());
            m_res_info->SetMaxRes(iEqu, max_res);
            m_res_info->SetAveRes(iEqu, ave_res);
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