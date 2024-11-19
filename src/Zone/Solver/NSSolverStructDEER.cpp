#include "NSSolverStructDEER.h"
namespace zaran
{
NSSolverStructDEER::NSSolverStructDEER(int index, string name, FlowSolverPara *para, GridStruct *grid,
                                       DataManagerNSStruct *data_manager)
    : NSSolverStruct(index, name, para, grid, data_manager)
{
}
NSSolverStructDEER::~NSSolverStructDEER()
{
}
void NSSolverStructDEER::CalcFluxDifference2ndOrder()
{
    auto para = GetPara();
    if (para->GetInterSchme() == InterpolationScheme::Grad)
    {
        FlueDifference2nd_MoveMidNode();
        return;
    }
    auto grid = GetGrid();
    auto data_manager = GetDataManager();
    auto node_metrics = GetNodeMetrics();
    auto idx_proxy = GetIdxProxy();
    auto gas = GetGas();
    auto equ_num = para->GetEqNum();
    int i_start, i_end, j_start, j_end, k_start, k_end;
    grid->GetRange(i_start, i_end, j_start, j_end, k_start, k_end);
    RiemannSolverPara riemann_para[6];

    // 当前节点的编号
    int idx;
    // 差分模板的值
    double value[5];
    double res_tmp[5];
#pragma omp parallel for private(idx, value, res_tmp, riemann_para)
    for (int k = k_start; k <= k_end; ++k)
    {
        for (int j = j_start; j <= j_end; ++j)
        {
            for (int i = i_start; i <= i_end; ++i)
            {
                for (int iRie = 0; iRie < 6; ++iRie)
                {
                    riemann_para[iRie].gamma_left = riemann_para[iRie].gamma_right = gas->GetGamma();
                }
                idx = idx_proxy->GetIdx(i, j, k);
                for (int iVal = 0; iVal < equ_num; ++iVal)
                {
                    res_tmp[iVal] = data_manager->GetResidual(iVal, idx);
                }
                // i direction
                auto xi = node_metrics->GetXi(idx);
                for (int iRie = 0; iRie < 2; ++iRie)
                {
                    riemann_para[iRie].norm(0) = xi[0];
                    riemann_para[iRie].norm(1) = xi[1];
                    riemann_para[iRie].norm(2) = xi[2];
                    riemann_para[iRie].nt = xi[3];
                }
                int idx_left = idx_proxy->GetIdx(i - 1, j, k);
                for (int iEqu = 0; iEqu < equ_num; ++iEqu)
                {
                    riemann_para[0].prim_left(iEqu) = data_manager->GetMidNodePrimLeft(iEqu, 0, idx);
                    riemann_para[0].prim_right(iEqu) = data_manager->GetMidNodePrimRight(iEqu, 0, idx);
                    riemann_para[1].prim_left(iEqu) = data_manager->GetMidNodePrimLeft(iEqu, 0, idx_left);
                    riemann_para[1].prim_right(iEqu) = data_manager->GetMidNodePrimRight(iEqu, 0, idx_left);
                }
                m_riemann_solver->Solver(riemann_para[0]);
                m_riemann_solver->Solver(riemann_para[1]);
                for (int iVal = 0; iVal < equ_num; ++iVal)
                {
                    res_tmp[iVal] -= (riemann_para[0].flux[iVal] - riemann_para[1].flux[iVal]);
                }
                // j direction
                auto eta = node_metrics->GetEta(idx);
                for (int iRie = 2; iRie < 4; ++iRie)
                {
                    riemann_para[iRie].norm(0) = eta[0];
                    riemann_para[iRie].norm(1) = eta[1];
                    riemann_para[iRie].norm(2) = eta[2];
                    riemann_para[iRie].nt = eta[3];
                }
                int idx_down = idx_proxy->GetIdx(i, j - 1, k);
                for (int iEqu = 0; iEqu < equ_num; ++iEqu)
                {
                    riemann_para[2].prim_left(iEqu) = data_manager->GetMidNodePrimLeft(iEqu, 1, idx);
                    riemann_para[2].prim_right(iEqu) = data_manager->GetMidNodePrimRight(iEqu, 1, idx);
                    riemann_para[3].prim_left(iEqu) = data_manager->GetMidNodePrimLeft(iEqu, 1, idx_down);
                    riemann_para[3].prim_right(iEqu) = data_manager->GetMidNodePrimRight(iEqu, 1, idx_down);
                }
                m_riemann_solver->Solver(riemann_para[2]);
                m_riemann_solver->Solver(riemann_para[3]);
                for (int iVal = 0; iVal < equ_num; ++iVal)
                {
                    res_tmp[iVal] -= (riemann_para[2].flux[iVal] - riemann_para[3].flux[iVal]);
                }
                // k direction
                if (grid->GetDim() == 3)
                {
                    auto zeta = node_metrics->GetZeta(idx);
                    for (int iRie = 4; iRie < 6; ++iRie)
                    {
                        riemann_para[iRie].norm(0) = zeta[0];
                        riemann_para[iRie].norm(1) = zeta[1];
                        riemann_para[iRie].norm(2) = zeta[2];
                        riemann_para[iRie].nt = zeta[3];
                    }
                    int idx_back = idx_proxy->GetIdx(i, j, k - 1);
                    for (int iEqu = 0; iEqu < equ_num; ++iEqu)
                    {
                        riemann_para[4].prim_left(iEqu) = data_manager->GetMidNodePrimLeft(iEqu, 2, idx);
                        riemann_para[4].prim_right(iEqu) = data_manager->GetMidNodePrimRight(iEqu, 2, idx);
                        riemann_para[5].prim_left(iEqu) = data_manager->GetMidNodePrimLeft(iEqu, 2, idx_back);
                        riemann_para[5].prim_right(iEqu) = data_manager->GetMidNodePrimRight(iEqu, 2, idx_back);
                    }
                    m_riemann_solver->Solver(riemann_para[4]);
                    m_riemann_solver->Solver(riemann_para[5]);
                    for (int iVal = 0; iVal < equ_num; ++iVal)
                    {
                        res_tmp[iVal] -= (riemann_para[4].flux[iVal] - riemann_para[5].flux[iVal]);
                    }
                }
                for (int iVar = 0; iVar < equ_num; ++iVar)
                {
                    data_manager->SetResidual(iVar, idx, res_tmp[iVar]);
                }
            }
        }
    }
}

void NSSolverStructDEER::FlueDifference2nd_MoveMidNode()
{
    auto grid = GetGrid();
    auto node = grid->GetNode();
    auto para = GetPara();
    auto data_manager = GetDataManager();
    auto node_metrics = GetNodeMetrics();
    auto idx_proxy = GetIdxProxy();
    auto gas = GetGas();
    auto equ_num = para->GetEqNum();
    int is, ie, js, je, ks, ke;
    grid->GetRange(is, ie, js, je, ks, ke);
    RiemannSolverPara riemann_para[6];
    for (int i = 0; i < 6; ++i)
    {
        riemann_para[i].gamma_left = riemann_para[i].gamma_right = gas->GetGamma();
    }
    // 当前节点的编号
    int idx;
    // 差分模板的编号
    int idx_temp[5];
    // 差分模板的值
    double value[5];
    double res_tmp[5];
    double mid_coord_left[3], mid_coord_right[3]; // i-1/2,i+1/2
    double mid_coord[3];                          // i-1/2,i+1/2的中点
    double move_vector[3];                        // mid_coord到i的向量用于移动mid_coord_left和mid_coord_right
#pragma omp parallel for private(idx, idx_temp, value, res_tmp, mid_coord_left, mid_coord_right, mid_coord,            \
                                     move_vector, riemann_para)
    for (int k = ks; k <= ke; ++k)
    {
        for (int j = js; j <= je; ++j)
        {
            for (int i = is; i <= ie; ++i)
            {
                for (int iVal = 0; iVal < equ_num; ++iVal)
                {
                    res_tmp[iVal] = data_manager->GetResidual(iVal, idx_proxy->GetIdx(i, j, k));
                }
                // i direction
                auto left_coord = node->GetCoord(i - 1, j, k);
                auto coord = node->GetCoord(i, j, k);
                auto right_coord = node->GetCoord(i + 1, j, k);
                for (int iDim = 0; iDim < 3; ++iDim)
                {
                    mid_coord_left[iDim] = (left_coord[iDim] + coord[iDim]) / 2;
                    mid_coord_right[iDim] = (right_coord[iDim] + coord[iDim]) / 2;
                    mid_coord[iDim] = (mid_coord_left[iDim] + mid_coord_right[iDim]) / 2;
                    move_vector[iDim] = coord[iDim] - mid_coord[iDim];
                    mid_coord_left[iDim] += move_vector[iDim];
                    mid_coord_right[iDim] += move_vector[iDim];
                }
                idx = idx_proxy->GetIdx(i, j, k);
                riemann_para[0].norm(0) = node_metrics->GetXi(idx)[0];
                riemann_para[0].norm(1) = node_metrics->GetXi(idx)[1];
                riemann_para[0].norm(2) = node_metrics->GetXi(idx)[2];
                riemann_para[0].nt = node_metrics->GetXi(idx)[3];
                riemann_para[1].norm(0) = node_metrics->GetXi(idx)[0];
                riemann_para[1].norm(1) = node_metrics->GetXi(idx)[1];
                riemann_para[1].norm(2) = node_metrics->GetXi(idx)[2];
                riemann_para[1].nt = node_metrics->GetXi(idx)[3];
                MidNodeGrad(idx_proxy->GetIdx(i, j, k), idx_proxy->GetIdx(i + 1, j, k), coord, mid_coord_right,
                            right_coord, riemann_para[0].prim_left.data(), riemann_para[0].prim_right.data());
                MidNodeGrad(idx_proxy->GetIdx(i - 1, j, k), idx_proxy->GetIdx(i, j, k), left_coord, mid_coord_left,
                            coord, riemann_para[1].prim_left.data(), riemann_para[1].prim_right.data());
                m_riemann_solver->Solver(riemann_para[0]);
                m_riemann_solver->Solver(riemann_para[1]);
                for (int iVal = 0; iVal < equ_num; ++iVal)
                {
                    res_tmp[iVal] -= (riemann_para[0].flux[iVal] - riemann_para[1].flux[iVal]);
                }
                // j direction
                left_coord = node->GetCoord(i, j - 1, k);
                coord = node->GetCoord(i, j, k);
                right_coord = node->GetCoord(i, j + 1, k);
                for (int iDim = 0; iDim < 3; ++iDim)
                {
                    mid_coord_left[iDim] = (left_coord[iDim] + coord[iDim]) / 2;
                    mid_coord_right[iDim] = (right_coord[iDim] + coord[iDim]) / 2;
                    mid_coord[iDim] = (mid_coord_left[iDim] + mid_coord_right[iDim]) / 2;
                    move_vector[iDim] = coord[iDim] - mid_coord[iDim];
                    mid_coord_left[iDim] += move_vector[iDim];
                    mid_coord_right[iDim] += move_vector[iDim];
                }
                riemann_para[2].norm(0) = node_metrics->GetEta(idx)[0];
                riemann_para[2].norm(1) = node_metrics->GetEta(idx)[1];
                riemann_para[2].norm(2) = node_metrics->GetEta(idx)[2];
                riemann_para[2].nt = node_metrics->GetEta(idx)[3];
                riemann_para[3].norm(0) = node_metrics->GetEta(idx)[0];
                riemann_para[3].norm(1) = node_metrics->GetEta(idx)[1];
                riemann_para[3].norm(2) = node_metrics->GetEta(idx)[2];
                riemann_para[3].nt = node_metrics->GetEta(idx)[3];
                MidNodeGrad(idx_proxy->GetIdx(i, j, k), idx_proxy->GetIdx(i, j + 1, k), coord, mid_coord_right,
                            right_coord, riemann_para[2].prim_left.data(), riemann_para[2].prim_right.data());
                MidNodeGrad(idx_proxy->GetIdx(i, j - 1, k), idx_proxy->GetIdx(i, j, k), left_coord, mid_coord_left,
                            coord, riemann_para[3].prim_left.data(), riemann_para[3].prim_right.data());
                m_riemann_solver->Solver(riemann_para[2]);
                m_riemann_solver->Solver(riemann_para[3]);
                for (int iVal = 0; iVal < equ_num; ++iVal)
                {
                    res_tmp[iVal] -= (riemann_para[2].flux[iVal] - riemann_para[3].flux[iVal]);
                }
                // k direction
                if (grid->GetDim() == 3)
                {
                    left_coord = node->GetCoord(i, j, k - 1);
                    coord = node->GetCoord(i, j, k);
                    right_coord = node->GetCoord(i, j, k + 1);
                    for (int iDim = 0; iDim < 3; ++iDim)
                    {
                        mid_coord_left[iDim] = (left_coord[iDim] + coord[iDim]) / 2;
                        mid_coord_right[iDim] = (right_coord[iDim] + coord[iDim]) / 2;
                        mid_coord[iDim] = (mid_coord_left[iDim] + mid_coord_right[iDim]) / 2;
                        move_vector[iDim] = coord[iDim] - mid_coord[iDim];
                        mid_coord_left[iDim] += move_vector[iDim];
                        mid_coord_right[iDim] += move_vector[iDim];
                    }
                    riemann_para[4].norm(0) = node_metrics->GetZeta(idx)[0];
                    riemann_para[4].norm(1) = node_metrics->GetZeta(idx)[1];
                    riemann_para[4].norm(2) = node_metrics->GetZeta(idx)[2];
                    riemann_para[4].nt = node_metrics->GetZeta(idx)[3];
                    riemann_para[5].norm(0) = node_metrics->GetZeta(idx)[0];
                    riemann_para[5].norm(1) = node_metrics->GetZeta(idx)[1];
                    riemann_para[5].norm(2) = node_metrics->GetZeta(idx)[2];
                    riemann_para[5].nt = node_metrics->GetZeta(idx)[3];
                    MidNodeGrad(idx_proxy->GetIdx(i, j, k), idx_proxy->GetIdx(i, j, k + 1), coord, mid_coord_right,
                                right_coord, riemann_para[4].prim_left.data(), riemann_para[4].prim_right.data());
                    MidNodeGrad(idx_proxy->GetIdx(i, j, k - 1), idx_proxy->GetIdx(i, j, k), left_coord, mid_coord_left,
                                coord, riemann_para[5].prim_left.data(), riemann_para[5].prim_right.data());
                    m_riemann_solver->Solver(riemann_para[4]);
                    m_riemann_solver->Solver(riemann_para[5]);
                    for (int iVal = 0; iVal < equ_num; ++iVal)
                    {
                        res_tmp[iVal] -= (riemann_para[4].flux[iVal] - riemann_para[5].flux[iVal]);
                    }
                }
                for (int iVar = 0; iVar < equ_num; ++iVar)
                {
                    data_manager->SetResidual(iVar, idx, res_tmp[iVar]);
                }
            }
        }
    }
}

void NSSolverStructDEER::CalcFluxDifference6thOrder()
{
    auto grid = GetGrid();
    auto para = GetPara();
    auto node = grid->GetNode();
    auto data_manager = GetDataManager();
    auto idx_proxy = GetIdxProxy();
    auto node_metrics = GetNodeMetrics();
    auto gas = GetGas();
    auto equ_num = para->GetEqNum();
    int is, ie, js, je, ks, ke;
    grid->GetRange(is, ie, js, je, ks, ke);
    RiemannSolverPara riemann_para[6];
    // 当前节点的编号
    int idx;
    // 差分模板的编号
    int idx_temp[5];
    // 差分模板的值
    double value[5];
    double res_tmp[5];
#pragma omp parallel for private(idx, idx_temp, value, res_tmp, riemann_para)
    for (int k = ks; k <= ke; ++k)
    {
        for (int j = js; j <= je; ++j)
        {
            for (int i = is; i <= ie; ++i)
            {
                for (int iRie = 0; iRie < 6; ++iRie)
                {
                    riemann_para[iRie].gamma_left = riemann_para[iRie].gamma_right = gas->GetGamma();
                }
                for (int iVal = 0; iVal < equ_num; ++iVal)
                {
                    res_tmp[iVal] = data_manager->GetResidual(iVal, idx_proxy->GetIdx(i, j, k));
                }
                // i direction
                idx = idx_proxy->GetIdx(i, j, k);
                for (int i = 0; i < 6; ++i)
                {
                    riemann_para[i].norm(0) = node_metrics->GetXi(idx)[0];
                    riemann_para[i].norm(1) = node_metrics->GetXi(idx)[1];
                    riemann_para[i].norm(2) = node_metrics->GetXi(idx)[2];
                    riemann_para[i].nt = node_metrics->GetXi(idx)[3];
                }
                for (int iVal = 0; iVal < equ_num; ++iVal)
                {
                    riemann_para[0].prim_left(iVal) =
                        data_manager->GetMidNodePrimLeft(iVal, 0, idx_proxy->GetIdx(i - 3, j, k));
                    riemann_para[0].prim_right(iVal) =
                        data_manager->GetMidNodePrimRight(iVal, 0, idx_proxy->GetIdx(i - 3, j, k));
                    riemann_para[1].prim_left(iVal) =
                        data_manager->GetMidNodePrimLeft(iVal, 0, idx_proxy->GetIdx(i - 2, j, k));
                    riemann_para[1].prim_right(iVal) =
                        data_manager->GetMidNodePrimRight(iVal, 0, idx_proxy->GetIdx(i - 2, j, k));
                    riemann_para[2].prim_left(iVal) =
                        data_manager->GetMidNodePrimLeft(iVal, 0, idx_proxy->GetIdx(i - 1, j, k));
                    riemann_para[2].prim_right(iVal) =
                        data_manager->GetMidNodePrimRight(iVal, 0, idx_proxy->GetIdx(i - 1, j, k));
                    riemann_para[3].prim_left(iVal) =
                        data_manager->GetMidNodePrimLeft(iVal, 0, idx_proxy->GetIdx(i, j, k));
                    riemann_para[3].prim_right(iVal) =
                        data_manager->GetMidNodePrimRight(iVal, 0, idx_proxy->GetIdx(i, j, k));
                    riemann_para[4].prim_left(iVal) =
                        data_manager->GetMidNodePrimLeft(iVal, 0, idx_proxy->GetIdx(i + 1, j, k));
                    riemann_para[4].prim_right(iVal) =
                        data_manager->GetMidNodePrimRight(iVal, 0, idx_proxy->GetIdx(i + 1, j, k));
                    riemann_para[5].prim_left(iVal) =
                        data_manager->GetMidNodePrimLeft(iVal, 0, idx_proxy->GetIdx(i + 2, j, k));
                    riemann_para[5].prim_right(iVal) =
                        data_manager->GetMidNodePrimRight(iVal, 0, idx_proxy->GetIdx(i + 2, j, k));
                }
                for (int i = 0; i < 6; ++i)
                {
                    m_riemann_solver->Solver(riemann_para[i]);
                }
                for (int iVar = 0; iVar < equ_num; ++iVar)
                {
                    res_tmp[iVar] -= 75.0 / 64.0 * (riemann_para[3].flux[iVar] - riemann_para[2].flux[iVar]) -
                                     25.0 / 384.0 * (riemann_para[4].flux[iVar] - riemann_para[1].flux[iVar]) +
                                     3.0 / 640.0 * (riemann_para[5].flux[iVar] - riemann_para[0].flux[iVar]);
                }

                // j direction
                for (int i = 0; i < 6; ++i)
                {
                    riemann_para[i].norm(0) = node_metrics->GetEta(idx)[0];
                    riemann_para[i].norm(1) = node_metrics->GetEta(idx)[1];
                    riemann_para[i].norm(2) = node_metrics->GetEta(idx)[2];
                    riemann_para[i].nt = node_metrics->GetEta(idx)[3];
                }
                for (int iVal = 0; iVal < equ_num; ++iVal)
                {
                    riemann_para[0].prim_left(iVal) =
                        data_manager->GetMidNodePrimLeft(iVal, 1, idx_proxy->GetIdx(i, j - 3, k));
                    riemann_para[0].prim_right(iVal) =
                        data_manager->GetMidNodePrimRight(iVal, 1, idx_proxy->GetIdx(i, j - 3, k));
                    riemann_para[1].prim_left(iVal) =
                        data_manager->GetMidNodePrimLeft(iVal, 1, idx_proxy->GetIdx(i, j - 2, k));
                    riemann_para[1].prim_right(iVal) =
                        data_manager->GetMidNodePrimRight(iVal, 1, idx_proxy->GetIdx(i, j - 2, k));
                    riemann_para[2].prim_left(iVal) =
                        data_manager->GetMidNodePrimLeft(iVal, 1, idx_proxy->GetIdx(i, j - 1, k));
                    riemann_para[2].prim_right(iVal) =
                        data_manager->GetMidNodePrimRight(iVal, 1, idx_proxy->GetIdx(i, j - 1, k));
                    riemann_para[3].prim_left(iVal) =
                        data_manager->GetMidNodePrimLeft(iVal, 1, idx_proxy->GetIdx(i, j, k));
                    riemann_para[3].prim_right(iVal) =
                        data_manager->GetMidNodePrimRight(iVal, 1, idx_proxy->GetIdx(i, j, k));
                    riemann_para[4].prim_left(iVal) =
                        data_manager->GetMidNodePrimLeft(iVal, 1, idx_proxy->GetIdx(i, j + 1, k));
                    riemann_para[4].prim_right(iVal) =
                        data_manager->GetMidNodePrimRight(iVal, 1, idx_proxy->GetIdx(i, j + 1, k));
                    riemann_para[5].prim_left(iVal) =
                        data_manager->GetMidNodePrimLeft(iVal, 1, idx_proxy->GetIdx(i, j + 2, k));
                    riemann_para[5].prim_right(iVal) =
                        data_manager->GetMidNodePrimRight(iVal, 1, idx_proxy->GetIdx(i, j + 2, k));
                }
                for (int i = 0; i < 6; ++i)
                {
                    m_riemann_solver->Solver(riemann_para[i]);
                }
                for (int iVar = 0; iVar < equ_num; ++iVar)
                {
                    res_tmp[iVar] -= 75.0 / 64.0 * (riemann_para[3].flux[iVar] - riemann_para[2].flux[iVar]) -
                                     25.0 / 384.0 * (riemann_para[4].flux[iVar] - riemann_para[1].flux[iVar]) +
                                     3.0 / 640.0 * (riemann_para[5].flux[iVar] - riemann_para[0].flux[iVar]);
                }
                // k direction
                if (grid->GetDim() == 3)
                {
                    for (int i = 0; i < 6; ++i)
                    {
                        riemann_para[i].norm(0) = node_metrics->GetZeta(idx)[0];
                        riemann_para[i].norm(1) = node_metrics->GetZeta(idx)[1];
                        riemann_para[i].norm(2) = node_metrics->GetZeta(idx)[2];
                        riemann_para[i].nt = node_metrics->GetZeta(idx)[3];
                    }
                    for (int iVal = 0; iVal < equ_num; ++iVal)
                    {
                        riemann_para[0].prim_left(iVal) =
                            data_manager->GetMidNodePrimLeft(iVal, 2, idx_proxy->GetIdx(i, j, k - 3));
                        riemann_para[0].prim_right(iVal) =
                            data_manager->GetMidNodePrimRight(iVal, 2, idx_proxy->GetIdx(i, j, k - 3));
                        riemann_para[1].prim_left(iVal) =
                            data_manager->GetMidNodePrimLeft(iVal, 2, idx_proxy->GetIdx(i, j, k - 2));
                        riemann_para[1].prim_right(iVal) =
                            data_manager->GetMidNodePrimRight(iVal, 2, idx_proxy->GetIdx(i, j, k - 2));
                        riemann_para[2].prim_left(iVal) =
                            data_manager->GetMidNodePrimLeft(iVal, 2, idx_proxy->GetIdx(i, j, k - 1));
                        riemann_para[2].prim_right(iVal) =
                            data_manager->GetMidNodePrimRight(iVal, 2, idx_proxy->GetIdx(i, j, k - 1));
                        riemann_para[3].prim_left(iVal) =
                            data_manager->GetMidNodePrimLeft(iVal, 2, idx_proxy->GetIdx(i, j, k));
                        riemann_para[3].prim_right(iVal) =
                            data_manager->GetMidNodePrimRight(iVal, 2, idx_proxy->GetIdx(i, j, k));
                        riemann_para[4].prim_left(iVal) =
                            data_manager->GetMidNodePrimLeft(iVal, 2, idx_proxy->GetIdx(i, j, k + 1));
                        riemann_para[4].prim_right(iVal) =
                            data_manager->GetMidNodePrimRight(iVal, 2, idx_proxy->GetIdx(i, j, k + 1));
                        riemann_para[5].prim_left(iVal) =
                            data_manager->GetMidNodePrimLeft(iVal, 2, idx_proxy->GetIdx(i, j, k + 2));
                        riemann_para[5].prim_right(iVal) =
                            data_manager->GetMidNodePrimRight(iVal, 2, idx_proxy->GetIdx(i, j, k + 2));
                    }
                    for (int i = 0; i < 6; ++i)
                    {
                        m_riemann_solver->Solver(riemann_para[i]);
                    }
                    for (int iVar = 0; iVar < equ_num; ++iVar)
                    {
                        res_tmp[iVar] -= 75.0 / 64.0 * (riemann_para[3].flux[iVar] - riemann_para[2].flux[iVar]) -
                                         25.0 / 384.0 * (riemann_para[4].flux[iVar] - riemann_para[1].flux[iVar]) +
                                         3.0 / 640.0 * (riemann_para[5].flux[iVar] - riemann_para[0].flux[iVar]);
                    }
                }
                for (int iVar = 0; iVar < equ_num; ++iVar)
                {
                    data_manager->SetResidual(iVar, idx, res_tmp[iVar]);
                }
            }
        }
    }
}
void NSSolverStructDEER::MidNodeGrad(int idx_left, int idx_right, const double *lef_coord, const double *mid_coord,
                                     const double *right_coord, double *value_left, double *value_right)
{
    int equ_num = GetPara()->GetEqNum();
    auto data_manager = GetDataManager();

    for (int iEqu = 0; iEqu < equ_num; ++iEqu)
    {
        value_left[iEqu] = data_manager->GetPrim(iEqu, idx_left);
        value_right[iEqu] = data_manager->GetPrim(iEqu, idx_right);

        for (int iDim = 0; iDim < 3; ++iDim)
        {
            double left_grad =
                data_manager->GetPrimGrad(iEqu, iDim, idx_left) * data_manager->GetLimiter(iEqu, idx_left);
            double right_grad =
                data_manager->GetPrimGrad(iEqu, iDim, idx_right) * data_manager->GetLimiter(iEqu, idx_right);

            value_left[iEqu] += (mid_coord[iDim] - lef_coord[iDim]) * left_grad;
            value_right[iEqu] += (mid_coord[iDim] - right_coord[iDim]) * right_grad;
        }
    }
}
} // namespace zaran