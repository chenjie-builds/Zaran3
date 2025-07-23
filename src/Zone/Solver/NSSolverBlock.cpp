#include "NSSolverBlock.h"
#include <omp.h>
namespace zaran
{
NSSolverBlock::NSSolverBlock(index_type index, string name, shared_ptr<FlowSolverParamStruct> para,
                             shared_ptr<GridBlock> grid, shared_ptr<DataManagerNSStruct> data_manager)
    : NSSolverStructDEER(index, name, para, grid, data_manager)
{
}

NSSolverBlock::~NSSolverBlock()
{
}

GridBlock *NSSolverBlock::GetGrid()
{
    return static_cast<GridBlock *>(NSSolver::GetGrid());
}

void NSSolverBlock::FluxDifference2nd()
{
    auto grid = GetGrid();
    auto node = grid->GetNode();
    auto para = GetPara();
    auto data_manager = GetDataManager();
    auto node_metrics = GetNodeMetrics();
    IdProxyStruct &idx_proxy = GetIdxProxy();
    auto gas = GetGas();
    auto equ_num = para->GetEqNum();
    index_type is, ie, js, je, ks, ke;
    grid->GetRange(is, ie, js, je, ks, ke);
    RiemannSolverPara riemann_para[2];
    double direction[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    double res_tmp[5];
#ifdef USE_OMP
#pragma omp parallel for collapse(3) private(res_tmp, riemann_para)
#endif // USE_OMP
    for (index_type k = ks; k <= ke; ++k)
    {
        for (index_type j = js; j <= je; ++j)
        {
            for (index_type i = is; i <= ie; ++i)
            {
                if (grid->GetIBlank(i, j, k) != IBlank::Fluid)
                {
                    continue;
                }
                int idx = idx_proxy(i, j, k);
                for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
                {
                    res_tmp[idx_eq] = data_manager->GetResidual(idx_eq, idx);
                }
                for (dimension_type dim = 0; dim < grid->GetDim(); dim++)
                {
                    for (int iTemp = 0; iTemp < 2; iTemp++)
                    {
                        riemann_para[iTemp].gamma_left = riemann_para[iTemp].gamma_right = gas->GetGamma();
                        riemann_para[iTemp].norm(0) = node_metrics->GetMetrics(dim, idx)[0];
                        riemann_para[iTemp].norm(1) = node_metrics->GetMetrics(dim, idx)[1];
                        riemann_para[iTemp].norm(2) = node_metrics->GetMetrics(dim, idx)[2];
                        riemann_para[iTemp].nt = node_metrics->GetMetrics(dim, idx)[3];
                        int idx_temp =
                            idx_proxy(i + (iTemp - 1) * direction[dim][0], j + (iTemp - 1) * direction[dim][1],
                                      k + (iTemp - 1) * direction[dim][2]);
                        for (int iVal = 0; iVal < equ_num; ++iVal)
                        {
                            riemann_para[iTemp].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, dim, idx_temp);
                            riemann_para[iTemp].prim_right(iVal) =
                                data_manager->GetMidNodePrimRight(iVal, dim, idx_temp);
                        }
                        m_riemann_solver->Solver(riemann_para[iTemp]);
                    }
                    for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
                    {
                        res_tmp[idx_eq] -= riemann_para[1].flux[idx_eq] - riemann_para[0].flux[idx_eq];
                    }
                }
                data_manager->SetResidual(idx, res_tmp);
            }
        }
    }
}

//void NSSolverBlock::InterMidNodePrim_Grad(int idx_left, int idx_right, const double *lef_coord, const double *mid_coord,
//                                          const double *right_coord, double *value_left, double *value_right)
//{
//    double coord_vec[3] = {0, 0, 0};
//    coord_vec[0] = right_coord[0] - lef_coord[0];
//    coord_vec[1] = right_coord[1] - lef_coord[1];
//    coord_vec[2] = right_coord[2] - lef_coord[2];
//    int equ_num = GetPara()->GetEqNum();
//    auto data_manager = GetDataManager();
//    for (int iVal = 0; iVal < equ_num; ++iVal)
//    {
//        value_left[iVal] =
//            data_manager->GetPrim(iVal, idx_left) + 0.5 * data_manager->GetLimiter(iVal, idx_left) *
//                                                        (coord_vec[0] * data_manager->GetPrimGrad(iVal, 0, idx_left) +
//                                                         coord_vec[1] * data_manager->GetPrimGrad(iVal, 1, idx_left) +
//                                                         coord_vec[2] * data_manager->GetPrimGrad(iVal, 2, idx_left));
//        value_right[iVal] =
//            data_manager->GetPrim(iVal, idx_right) - 0.5 * data_manager->GetLimiter(iVal, idx_right) *
//                                                         (coord_vec[0] * data_manager->GetPrimGrad(iVal, 0, idx_right) +
//                                                          coord_vec[1] * data_manager->GetPrimGrad(iVal, 1, idx_right) +
//                                                          coord_vec[2] * data_manager->GetPrimGrad(iVal, 2, idx_right));
//    }
//    if (value_left[0] < 0 || value_left[4] < 0)
//    {
//        for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
//        {
//            value_left[idx_eq] = data_manager->GetPrim(idx_eq, idx_left);
//        }
//    }
//    if (value_right[0] < 0 || value_right[4] < 0)
//    {
//        for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
//        {
//            value_right[idx_eq] = data_manager->GetPrim(idx_eq, idx_right);
//        }
//    }
//}

//void NSSolverBlock::CalcConvectionRes_Grad()
//{
//    FluxDifference2nd();
//}

void NSSolverBlock::Preprocess()
{
    NSSolver::Preprocess();
}

} // namespace zaran