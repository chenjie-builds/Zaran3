#include "NSFieldStruct.h"
#include "DataManagerNSStruct.h"
#include "FlowSolverStructPara.h"
#include "NSSolverStructCMM.h"
#include "NSSolverStructDEER.h"
#include "Log.h"
namespace zaran {
NSFieldStruct::NSFieldStruct(GridBase* grid)
    : FieldNS(grid, FieldType::NS_Structured)
{
    m_idx_proxy = nullptr;
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
FlowSolverStructPara* NSFieldStruct::GetSolverPara()
{
    return static_cast<FlowSolverStructPara*>(Field::GetSolverPara());
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
    int equ_num = para->GetEqNum();
    int is, ie, js, je, ks, ke;
    grid->GetRange(is, ie, js, je, ks, ke);
    int total_node_num = (ie - is + 1) * (je - js + 1) * (ke - ks + 1);
    double norm_inf = -LARGE_NUMBER;
    double norm_l2 = 0;
    int norm_inf_node = -1;
    double norm_inf_coord[3];
    for (int iEqu = 0; iEqu < equ_num; iEqu++) {
        norm_inf = -LARGE_NUMBER;
        norm_l2 = 0;
        auto res = GetDataManager()->GetResidual(iEqu);
        #pragma omp parallel for reduction(max : norm_inf) reduction(+ : norm_l2)
        for (int k = ks; k <= ke; k++) {
            for (int j = js; j <= je; j++) {
                for (int i = is; i <= ie; i++) {
                    int idx = m_idx_proxy->GetIdx(i, j, k);
                    if (abs(res[idx]) > norm_inf) {
                        norm_inf = abs(res[idx]);
                        norm_inf_node = idx;
                        for (int iDim = 0; iDim < grid->GetDim(); iDim++) {
                            norm_inf_coord[iDim] = node->GetCoord(i, j, k)[iDim];
                        }
                    }
                    norm_l2 += res[idx] * res[idx];
                }
            }
        }
        norm_l2 = sqrt(norm_l2 / total_node_num);
        m_res_info->SetInfNorm(iEqu, norm_inf);
        m_res_info->SetL2Norm(iEqu, norm_l2);
        m_res_info->SetInfNormCoord(iEqu, norm_inf_coord);
        m_res_info->SetInfNormIdx(iEqu, norm_inf_node);
    }
}
void NSFieldStruct::Allocate()
{
    FieldNS::Allocate();
    AllocateIdxProxy();
}
void NSFieldStruct::AllocateResInfo()
{
    if (m_res_info != nullptr) {
        delete m_res_info;
    }
    m_res_info = new ResInfo(GetSolverPara()->GetEqNum());
}
void NSFieldStruct::AllocateSolver()
{
    if (m_solver != nullptr) {
        delete m_solver;
    }
    auto para = GetSolverPara();
    if (para->GetMidMetricsScheme() == MidMetricsScheme::CMM) {
        m_solver = new NSSolverStructCMM(1, "NS_Struct", GetSolverPara(), GetGrid(), GetDataManager());
    } else if (para->GetMidMetricsScheme() == MidMetricsScheme::DEER) {
        m_solver = new NSSolverStructDEER(1, "NS_Struct", GetSolverPara(), GetGrid(), GetDataManager());
    } else {
        Log::error("Invalid metrics scheme");
    }
}
void NSFieldStruct::AllocateDataManager()
{
    if (m_data_manager != nullptr) {
        delete m_data_manager;
    }
    int ni = GetGrid()->GetNi();
    int nj = GetGrid()->GetNj();
    int nk = GetGrid()->GetNk();
    m_data_manager = new DataManagerNSStruct(GetFieldData(), ni, nj, nk);
    m_data_manager->CreateData();
    m_data_manager->RegisterData();
}
void NSFieldStruct::AllocateIdxProxy()
{
    if (m_idx_proxy != nullptr) {
        delete m_idx_proxy;
    }
    m_idx_proxy = new StructIdxProxy(GetGrid());
}
void NSFieldStruct::AllocateSolverPara()
{
    if (m_solver_para != nullptr) {
        delete m_solver_para;
    }
    m_solver_para = new FlowSolverStructPara();
    GetSolverPara()->Init();
}
} // namespace zaran