#include "NSSolverStructDEER.h"
#include "Log.h"
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
    void NSSolverStructDEER::FluxDifference2nd()
    {
        auto grid = GetGrid();
        auto para = GetPara();
        auto data_manager = GetDataManager();
        auto node_metrics = GetNodeMetrics();
        auto idx_proxy = GetIdxProxy();
        auto gas = GetGas();
        auto equ_num = para->GetEquNum();
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
                    idx = idx_proxy->GetIdx(i, j, k);
                    riemann_para[0].norm(0) = node_metrics->GetXi(idx)[0];
                    riemann_para[0].norm(1) = node_metrics->GetXi(idx)[1];
                    riemann_para[0].norm(2) = node_metrics->GetXi(idx)[2];
                    riemann_para[0].nt = node_metrics->GetXi(idx)[3];
                    riemann_para[1].norm(0) = node_metrics->GetXi(idx)[0];
                    riemann_para[1].norm(1) = node_metrics->GetXi(idx)[1];
                    riemann_para[1].norm(2) = node_metrics->GetXi(idx)[2];
                    riemann_para[1].nt = node_metrics->GetXi(idx)[3];
                    for (int iVal = 0; iVal < equ_num; ++iVal)
                    {
                        riemann_para[0].prim_left(iVal) =
                            data_manager->GetMidNodePrimLeft(iVal, 0, idx_proxy->GetIdx(i, j, k));
                        riemann_para[0].prim_right(iVal) =
                            data_manager->GetMidNodePrimRight(iVal, 0, idx_proxy->GetIdx(i, j, k));
                        riemann_para[1].prim_left(iVal) =
                            data_manager->GetMidNodePrimLeft(iVal, 0, idx_proxy->GetIdx(i - 1, j, k));
                        riemann_para[1].prim_right(iVal) =
                            data_manager->GetMidNodePrimRight(iVal, 0, idx_proxy->GetIdx(i - 1, j, k));
                    }
                    m_riemann_solver->Solver(riemann_para[0]);
                    m_riemann_solver->Solver(riemann_para[1]);
                    for (int iVal = 0; iVal < equ_num; ++iVal)
                    {
                        res_tmp[iVal] -= (riemann_para[0].flux[iVal] - riemann_para[1].flux[iVal]);
                    }
                    // j direction
                    riemann_para[2].norm(0) = node_metrics->GetEta(idx)[0];
                    riemann_para[2].norm(1) = node_metrics->GetEta(idx)[1];
                    riemann_para[2].norm(2) = node_metrics->GetEta(idx)[2];
                    riemann_para[2].nt = node_metrics->GetEta(idx)[3];
                    riemann_para[3].norm(0) = node_metrics->GetEta(idx)[0];
                    riemann_para[3].norm(1) = node_metrics->GetEta(idx)[1];
                    riemann_para[3].norm(2) = node_metrics->GetEta(idx)[2];
                    riemann_para[3].nt = node_metrics->GetEta(idx)[3];
                    for (int iVal = 0; iVal < equ_num; ++iVal)
                    {
                        riemann_para[2].prim_left(iVal) =
                            data_manager->GetMidNodePrimLeft(iVal, 1, idx_proxy->GetIdx(i, j, k));
                        riemann_para[2].prim_right(iVal) =
                            data_manager->GetMidNodePrimRight(iVal, 1, idx_proxy->GetIdx(i, j, k));
                        riemann_para[3].prim_left(iVal) =
                            data_manager->GetMidNodePrimLeft(iVal, 1, idx_proxy->GetIdx(i, j - 1, k));
                        riemann_para[3].prim_right(iVal) =
                            data_manager->GetMidNodePrimRight(iVal, 1, idx_proxy->GetIdx(i, j - 1, k));
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
                        riemann_para[4].norm(0) = node_metrics->GetZeta(idx)[0];
                        riemann_para[4].norm(1) = node_metrics->GetZeta(idx)[1];
                        riemann_para[4].norm(2) = node_metrics->GetZeta(idx)[2];
                        riemann_para[4].nt = node_metrics->GetZeta(idx)[3];
                        riemann_para[5].norm(0) = node_metrics->GetZeta(idx)[0];
                        riemann_para[5].norm(1) = node_metrics->GetZeta(idx)[1];
                        riemann_para[5].norm(2) = node_metrics->GetZeta(idx)[2];
                        riemann_para[5].nt = node_metrics->GetZeta(idx)[3];
                        idx_temp[0] = idx_proxy->GetIdx(i, j, k - 2);
                        idx_temp[1] = idx_proxy->GetIdx(i, j, k - 1);
                        idx_temp[2] = idx_proxy->GetIdx(i, j, k);
                        idx_temp[3] = idx_proxy->GetIdx(i, j, k + 1);
                        idx_temp[4] = idx_proxy->GetIdx(i, j, k + 2);
                        for (int iVal = 0; iVal < equ_num; ++iVal)
                        {
                            riemann_para[4].prim_left(iVal) =
                                data_manager->GetMidNodePrimLeft(iVal, 2, idx_proxy->GetIdx(i, j, k));
                            riemann_para[4].prim_right(iVal) =
                                data_manager->GetMidNodePrimRight(iVal, 2, idx_proxy->GetIdx(i, j, k));
                            riemann_para[5].prim_left(iVal) =
                                data_manager->GetMidNodePrimLeft(iVal, 2, idx_proxy->GetIdx(i, j, k - 1));
                            riemann_para[5].prim_right(iVal) =
                                data_manager->GetMidNodePrimRight(iVal, 2, idx_proxy->GetIdx(i, j, k - 1));
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

    void NSSolverStructDEER::FluxDifference6th()
    {
        auto grid = GetGrid();
        auto para = GetPara();
        auto node = grid->GetNode();
        auto data_manager = GetDataManager();
        auto idx_proxy = GetIdxProxy();
        auto node_metrics = GetNodeMetrics();
        auto gas = GetGas();
        auto equ_num = para->GetEquNum();
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
} // namespace zaran