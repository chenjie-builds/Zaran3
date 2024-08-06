#include "NSFieldStruct.h"
#include"DataManagerNSStruct.h"
namespace zaran
{
    NSFieldStruct::NSFieldStruct(GridBase* grid) :FieldNS(grid, FieldType::NS_Structured)
    {
        auto grid_struct = static_cast<GridStruct*>(grid);
        int ni = grid_struct->GetNi();
        int nj = grid_struct->GetNj();
        int nk = grid_struct->GetNk();
        m_data_manager = new DataManagerNSStruct(GetFieldData(), ni, nj, nk);
        m_data_manager->CreateData();
        m_data_manager->RegisterData();
        m_solver = new NSSolverStruct(1, "NS_Struct", GetSolverPara(), GetGrid(), GetDataManager());
        m_res_info = new ResInfo(GetSolverPara()->GetEquNum());
        m_idx_proxy = new StructIdxProxy(GetGrid());
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
    DataManagerNSStruct* NSFieldStruct::GetDataManager()
    {
        return static_cast<DataManagerNSStruct*>(Field::GetDataManager());
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
        if (m_data_manager != nullptr)
        {
            delete m_data_manager;
        }
        m_data_manager = new DataManagerNS(GetFieldData(), GetGrid()->GetTotalNodeNum());
        if (m_solver != nullptr)
        {
            delete m_solver;
        }
        m_solver = new NSSolverStruct(1, "NS_Struct", GetSolverPara(), GetGrid(), GetDataManager());
        if (m_res_info != nullptr)
        {
            delete m_res_info;
        }
        m_res_info = new ResInfo(GetSolverPara()->GetEquNum());
        if (m_idx_proxy != nullptr)
        {
            delete m_idx_proxy;
        }
        m_idx_proxy = new StructIdxProxy(GetGrid());
    }
} // namespace zaran