#include "NSSolverStructCMM.h"
#include "NSSolverStruct.h"
#include "NodeMetricsFN.h"
#include "Log.h"
namespace zaran
{
    NSSolverStructCMM::NSSolverStructCMM(int index, string name, FlowSolverPara *para, GridStruct *grid,
                                         DataManagerNSStruct *data_manager)
        : NSSolverStruct(index, name, para, grid, data_manager)
    {
    }
    NSSolverStructCMM::~NSSolverStructCMM()
    {
    }

    void NSSolverStructCMM::CalcCoordTransCoef()
    {
        NSSolverStruct::CalcCoordTransCoef();
        CalcMidNodeMetrics();
    }
    void NSSolverStructCMM::CalcMidNodeMetrics()
    {
        auto para = GetPara();
        auto flux_diff_scheme = para->GetDifferenceScheme();
        if (flux_diff_scheme == DifferenceScheme::SecondOrder)
        {
            CalcMidNodeMetrics2nd();
        }
        else if (flux_diff_scheme == DifferenceScheme::SixthOrder)
        {
            CalcMidNodeMetrics6th();
        }
        else
        {
            CalcMidNodeMetrics2nd();
        }
    }
    void NSSolverStructCMM::CalcMidNodeMetrics2nd()
    {
        auto grid = GetGrid();
        auto node = grid->GetNode();
        auto coef = GetNodeMetrics();
        Metrics *coef_mid[3];
        coef_mid[0] = GetMidMetricsI();
        coef_mid[1] = GetMidMetricsJ();
        coef_mid[2] = GetMidMetricsK();
        auto idx_proxy = GetIdxProxy();
        int ni = grid->GetNi();
        int nj = grid->GetNj();
        int nk = grid->GetNk();
        int idx_temp[3];
        // 用于i,j,k方向寻找模板的系数
        int idx_direct[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        // 内部半点：根据整点度量系数计算半点度量系数
        for (int k = 1; k < nk - 2; ++k)
        {
            for (int j = 1; j < nj - 2; ++j)
            {
                for (int i = 1; i < ni - 2; ++i)
                {
                    int idx = idx_proxy->GetIdx(i, j, k);
                    for (int iDim = 0; iDim < 3; iDim++)
                    {
                        auto coef_xi = coef_mid[iDim]->GetXi(idx);
                        auto coef_eta = coef_mid[iDim]->GetEta(idx);
                        auto coef_zeta = coef_mid[iDim]->GetZeta(idx);
                        idx_temp[0] = idx_proxy->GetIdx(i, j, k);
                        idx_temp[1] = idx_proxy->GetIdx(i + 1 * idx_direct[iDim][0], j + 1 * idx_direct[iDim][1],
                                                        k + 1 * idx_direct[iDim][2]);
                        for (int jDim = 0; jDim < 3; ++jDim)
                        {
                            coef_xi[jDim] =
                                0.5 * (coef->GetXi(idx_temp[0])[jDim] + coef->GetXi(idx_temp[1])[jDim]);
                            coef_eta[jDim] =
                                0.5 * (coef->GetEta(idx_temp[0])[jDim] + coef->GetEta(idx_temp[1])[jDim]);
                            coef_zeta[jDim] =
                                0.5 * (coef->GetZeta(idx_temp[0])[jDim] + coef->GetZeta(idx_temp[1])[jDim]);
                        }
                    }
                }
            }
        }
        // 边界半点：根据整点度量系数计算半点度量系数
        // i=1/2，i=(ni-1)-1/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                auto coef_xi_left = coef_mid[0]->GetXi(idx_proxy->GetIdx(0, j, k));           // xi[1/2]
                auto coef_eta_left = coef_mid[0]->GetEta(idx_proxy->GetIdx(0, j, k));         // eta[1/2]
                auto coef_zeta_left = coef_mid[0]->GetZeta(idx_proxy->GetIdx(0, j, k));       // zeta[1/2]
                auto coef_xi_right = coef_mid[0]->GetXi(idx_proxy->GetIdx(ni - 2, j, k));     // xi[(ni-1)-1/2]
                auto coef_eta_right = coef_mid[0]->GetEta(idx_proxy->GetIdx(ni - 2, j, k));   // eta[(ni-1)-1/2]
                auto coef_zeta_right = coef_mid[0]->GetZeta(idx_proxy->GetIdx(ni - 2, j, k)); // zeta[(ni-1)-1/2]
                for (int iDim = 0; iDim < 3; ++iDim)
                {
                    coef_xi_left[iDim] = 0.5 * (3.0 * coef->GetXi(idx_proxy->GetIdx(1, j, k))[iDim] -
                                                coef->GetXi(idx_proxy->GetIdx(2, j, k))[iDim]);
                    coef_eta_left[iDim] = 0.5 * (3.0 * coef->GetEta(idx_proxy->GetIdx(1, j, k))[iDim] -
                                                 coef->GetEta(idx_proxy->GetIdx(2, j, k))[iDim]);
                    coef_zeta_left[iDim] = 0.5 * (3.0 * coef->GetZeta(idx_proxy->GetIdx(1, j, k))[iDim] -
                                                  coef->GetZeta(idx_proxy->GetIdx(2, j, k))[iDim]);
                    coef_xi_right[iDim] = 0.5 * (3.0 * coef->GetXi(idx_proxy->GetIdx(ni - 2, j, k))[iDim] -
                                                 coef->GetXi(idx_proxy->GetIdx(ni - 3, j, k))[iDim]);
                    coef_eta_right[iDim] = 0.5 * (3.0 * coef->GetEta(idx_proxy->GetIdx(ni - 2, j, k))[iDim] -
                                                  coef->GetEta(idx_proxy->GetIdx(ni - 3, j, k))[iDim]);
                    coef_zeta_right[iDim] = 0.5 * (3.0 * coef->GetZeta(idx_proxy->GetIdx(ni - 2, j, k))[iDim] -
                                                   coef->GetZeta(idx_proxy->GetIdx(ni - 3, j, k))[iDim]);
                }
            }
        }
        // j=1/2,j=(nj-1)-1/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                auto coef_xi_left = coef_mid[1]->GetXi(idx_proxy->GetIdx(i, 0, k));           // xi[1/2]
                auto coef_eta_left = coef_mid[1]->GetEta(idx_proxy->GetIdx(i, 0, k));         // eta[1/2]
                auto coef_zeta_left = coef_mid[1]->GetZeta(idx_proxy->GetIdx(i, 0, k));       // zeta[1/2]
                auto coef_xi_right = coef_mid[1]->GetXi(idx_proxy->GetIdx(i, nj - 2, k));     // xi[(nj-1)-1/2]
                auto coef_eta_right = coef_mid[1]->GetEta(idx_proxy->GetIdx(i, nj - 2, k));   // eta[(nj-1)-1/2]
                auto coef_zeta_right = coef_mid[1]->GetZeta(idx_proxy->GetIdx(i, nj - 2, k)); // zeta[(nj-1)-1/2]
                for (int iDim = 0; iDim < 3; ++iDim)
                {
                    coef_xi_left[iDim] = 0.5 * (3.0 * coef->GetXi(idx_proxy->GetIdx(i, 1, k))[iDim] -
                                                coef->GetXi(idx_proxy->GetIdx(i, 2, k))[iDim]);
                    coef_eta_left[iDim] = 0.5 * (3.0 * coef->GetEta(idx_proxy->GetIdx(i, 1, k))[iDim] -
                                                 coef->GetEta(idx_proxy->GetIdx(i, 2, k))[iDim]);
                    coef_zeta_left[iDim] = 0.5 * (3.0 * coef->GetZeta(idx_proxy->GetIdx(i, 1, k))[iDim] -
                                                  coef->GetZeta(idx_proxy->GetIdx(i, 2, k))[iDim]);
                    coef_xi_right[iDim] = 0.5 * (3.0 * coef->GetXi(idx_proxy->GetIdx(i, nj - 2, k))[iDim] -
                                                 coef->GetXi(idx_proxy->GetIdx(i, nj - 3, k))[iDim]);
                    coef_eta_right[iDim] = 0.5 * (3.0 * coef->GetEta(idx_proxy->GetIdx(i, nj - 2, k))[iDim] -
                                                  coef->GetEta(idx_proxy->GetIdx(i, nj - 3, k))[iDim]);
                    coef_zeta_right[iDim] = 0.5 * (3.0 * coef->GetZeta(idx_proxy->GetIdx(i, nj - 2, k))[iDim] -
                                                   coef->GetZeta(idx_proxy->GetIdx(i, nj - 3, k))[iDim]);
                }
            }
        }
        // k=1/2,k=(nk-1)-1/2
        for (int j = 1; j < nj - 1; ++j)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                auto coef_xi_left = coef_mid[2]->GetXi(idx_proxy->GetIdx(i, j, 0));           // xi[1/2]
                auto coef_eta_left = coef_mid[2]->GetEta(idx_proxy->GetIdx(i, j, 0));         // eta[1/2]
                auto coef_zeta_left = coef_mid[2]->GetZeta(idx_proxy->GetIdx(i, j, 0));       // zeta[1/2]
                auto coef_xi_right = coef_mid[2]->GetXi(idx_proxy->GetIdx(i, j, nk - 2));     // xi[(nk-1)-1/2]
                auto coef_eta_right = coef_mid[2]->GetEta(idx_proxy->GetIdx(i, j, nk - 2));   // eta[(nk-1)-1/2]
                auto coef_zeta_right = coef_mid[2]->GetZeta(idx_proxy->GetIdx(i, j, nk - 2)); // zeta[(nk-1)-1/2]
                for (int iDim = 0; iDim < 3; ++iDim)
                {
                    coef_xi_left[iDim] = 0.5 * (3.0 * coef->GetXi(idx_proxy->GetIdx(i, j, 1))[iDim] -
                                                coef->GetXi(idx_proxy->GetIdx(i, j, 2))[iDim]);
                    coef_eta_left[iDim] = 0.5 * (3.0 * coef->GetEta(idx_proxy->GetIdx(i, j, 1))[iDim] -
                                                 coef->GetEta(idx_proxy->GetIdx(i, j, 2))[iDim]);
                    coef_zeta_left[iDim] = 0.5 * (3.0 * coef->GetZeta(idx_proxy->GetIdx(i, j, 1))[iDim] -
                                                  coef->GetZeta(idx_proxy->GetIdx(i, j, 2))[iDim]);
                    coef_xi_right[iDim] = 0.5 * (3.0 * coef->GetXi(idx_proxy->GetIdx(i, j, nk - 2))[iDim] -
                                                 coef->GetXi(idx_proxy->GetIdx(i, j, nk - 3))[iDim]);
                    coef_eta_right[iDim] = 0.5 * (3.0 * coef->GetEta(idx_proxy->GetIdx(i, j, nk - 2))[iDim] -
                                                  coef->GetEta(idx_proxy->GetIdx(i, j, nk - 3))[iDim]);
                    coef_zeta_right[iDim] = 0.5 * (3.0 * coef->GetZeta(idx_proxy->GetIdx(i, j, nk - 2))[iDim] -
                                                   coef->GetZeta(idx_proxy->GetIdx(i, j, nk - 3))[iDim]);
                }
            }
        }
    }
    void NSSolverStructCMM::CalcMidNodeMetrics6th()
    {
        auto grid = GetGrid();
        auto node = grid->GetNode();
        auto coef = GetNodeMetrics();
        Metrics *coef_mid[3] = {GetMidMetricsI(), GetMidMetricsJ(), GetMidMetricsK()};
        auto idx_proxy = GetIdxProxy();
        auto Idx = [&](int i, int j, int k)
        {
            return idx_proxy->GetIdx(i, j, k);
        };
        int ni = grid->GetNi();
        int nj = grid->GetNj();
        int nk = grid->GetNk();
        int idx_temp[6];
        // 内部半点：根据整点度量系数计算半点度量系数
        // 第三步：根据整点逆变换度量系数计算半点逆变换度量系数
        double inter_temp[6];
        // i+1/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 3; i < ni - 3; ++i)
                {
                    for (int iDim = 0; iDim < 3; iDim++)
                    {
                        for (int iTemp = 0; iTemp < 6; iTemp++)
                        {
                            inter_temp[iTemp] = coef->GetXi(Idx(i + iTemp - 2, j, k))[iDim];
                        }
                        coef_mid[0]->GetXi(Idx(i, j, k))[iDim] = MidNodeInter6th(inter_temp);
                        for (int iTemp = 0; iTemp < 6; iTemp++)
                        {
                            inter_temp[iTemp] = coef->GetEta(Idx(i + iTemp - 2, j, k))[iDim];
                        }
                        coef_mid[0]->GetEta(Idx(i, j, k))[iDim] = MidNodeInter6th(inter_temp);
                        for (int iTemp = 0; iTemp < 6; iTemp++)
                        {
                            inter_temp[iTemp] = coef->GetZeta(Idx(i + iTemp - 2, j, k))[iDim];
                        }
                        coef_mid[0]->GetZeta(Idx(i, j, k))[iDim] = MidNodeInter6th(inter_temp);
                    }
                }
            }
        }
        // j+1/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 3; j < nj - 3; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    for (int iDim = 0; iDim < 3; iDim++)
                    {
                        for (int iTemp = 0; iTemp < 6; iTemp++)
                        {
                            inter_temp[iTemp] = coef->GetXi(Idx(i, j + iTemp - 2, k))[iDim];
                        }
                        coef_mid[1]->GetXi(Idx(i, j, k))[iDim] = MidNodeInter6th(inter_temp);
                        for (int iTemp = 0; iTemp < 6; iTemp++)
                        {
                            inter_temp[iTemp] = coef->GetEta(Idx(i, j + iTemp - 2, k))[iDim];
                        }
                        coef_mid[1]->GetEta(Idx(i, j, k))[iDim] = MidNodeInter6th(inter_temp);
                        for (int iTemp = 0; iTemp < 6; iTemp++)
                        {
                            inter_temp[iTemp] = coef->GetZeta(Idx(i, j + iTemp - 2, k))[iDim];
                        }
                        coef_mid[1]->GetZeta(Idx(i, j, k))[iDim] = MidNodeInter6th(inter_temp);
                    }
                }
            }
        }
        // k+1/2
        for (int k = 3; k < nk - 3; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    for (int iDim = 0; iDim < 3; iDim++)
                    {
                        for (int iTemp = 0; iTemp < 6; iTemp++)
                        {
                            inter_temp[iTemp] = coef->GetXi(Idx(i, j, k + iTemp - 2))[iDim];
                        }
                        coef_mid[2]->GetXi(Idx(i, j, k))[iDim] = MidNodeInter6th(inter_temp);
                        for (int iTemp = 0; iTemp < 6; iTemp++)
                        {
                            inter_temp[iTemp] = coef->GetEta(Idx(i, j, k + iTemp - 2))[iDim];
                        }
                        coef_mid[2]->GetEta(Idx(i, j, k))[iDim] = MidNodeInter6th(inter_temp);
                        for (int iTemp = 0; iTemp < 6; iTemp++)
                        {
                            inter_temp[iTemp] = coef->GetZeta(Idx(i, j, k + iTemp - 2))[iDim];
                        }
                        coef_mid[2]->GetZeta(Idx(i, j, k))[iDim] = MidNodeInter6th(inter_temp);
                    }
                }
            }
        }
        // 第四步：计算边界半点的逆变换度量系数(二阶偏置插值)
        // i=1/2,3/2,5/2;ni-3/2,ni-5/2,ni-7/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int iDim = 0; iDim < 3; ++iDim)
                {
                    // i=1/2，3/2，5/2
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetXi(Idx(1 + iTemp, j, k))[iDim];
                    }
                    coef_mid[0]->GetXi(Idx(0, j, k))[iDim] = MidNodeInter4thRight1(inter_temp);
                    coef_mid[0]->GetXi(Idx(1, j, k))[iDim] = MidNodeInter4thRight2(inter_temp);
                    coef_mid[0]->GetXi(Idx(2, j, k))[iDim] = MidNodeInter4th(inter_temp);
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetEta(Idx(1 + iTemp, j, k))[iDim];
                    }
                    coef_mid[0]->GetEta(Idx(0, j, k))[iDim] = MidNodeInter4thRight1(inter_temp);
                    coef_mid[0]->GetEta(Idx(1, j, k))[iDim] = MidNodeInter4thRight2(inter_temp);
                    coef_mid[0]->GetEta(Idx(2, j, k))[iDim] = MidNodeInter4th(inter_temp);
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZeta(Idx(1 + iTemp, j, k))[iDim];
                    }
                    coef_mid[0]->GetZeta(Idx(0, j, k))[iDim] = MidNodeInter4thRight1(inter_temp);
                    coef_mid[0]->GetZeta(Idx(1, j, k))[iDim] = MidNodeInter4thRight2(inter_temp);
                    coef_mid[0]->GetZeta(Idx(2, j, k))[iDim] = MidNodeInter4th(inter_temp);
                    // i=ni-3/2,ni-5/2,ni-7/2
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetXi(Idx(ni - 5 + iTemp, j, k))[iDim];
                    }
                    coef_mid[0]->GetXi(Idx(ni - 2, j, k))[iDim] = MidNodeInter4thLeft1(inter_temp);
                    coef_mid[0]->GetXi(Idx(ni - 3, j, k))[iDim] = MidNodeInter4thLeft2(inter_temp);
                    coef_mid[0]->GetXi(Idx(ni - 4, j, k))[iDim] = MidNodeInter4th(inter_temp);
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetEta(Idx(ni - 5 + iTemp, j, k))[iDim];
                    }
                    coef_mid[0]->GetEta(Idx(ni - 2, j, k))[iDim] = MidNodeInter4thLeft1(inter_temp);
                    coef_mid[0]->GetEta(Idx(ni - 3, j, k))[iDim] = MidNodeInter4thLeft2(inter_temp);
                    coef_mid[0]->GetEta(Idx(ni - 4, j, k))[iDim] = MidNodeInter4th(inter_temp);
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZeta(Idx(ni - 5 + iTemp, j, k))[iDim];
                    }
                    coef_mid[0]->GetZeta(Idx(ni - 2, j, k))[iDim] = MidNodeInter4thLeft1(inter_temp);
                    coef_mid[0]->GetZeta(Idx(ni - 3, j, k))[iDim] = MidNodeInter4thLeft2(inter_temp);
                    coef_mid[0]->GetZeta(Idx(ni - 4, j, k))[iDim] = MidNodeInter4th(inter_temp);
                }
            }
        }
        // j=1/2,3/2,5/2;nj-3/2,nj-5/2,nj-7/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                for (int iDim = 0; iDim < 3; ++iDim)
                {
                    // j=1/2，3/2，5/2
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetXi(Idx(i, 1 + iTemp, k))[iDim];
                    }
                    coef_mid[1]->GetXi(Idx(i, 0, k))[iDim] = MidNodeInter4thRight1(inter_temp);
                    coef_mid[1]->GetXi(Idx(i, 1, k))[iDim] = MidNodeInter4thRight2(inter_temp);
                    coef_mid[1]->GetXi(Idx(i, 2, k))[iDim] = MidNodeInter4th(inter_temp);
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetEta(Idx(i, 1 + iTemp, k))[iDim];
                    }
                    coef_mid[1]->GetEta(Idx(i, 0, k))[iDim] = MidNodeInter4thRight1(inter_temp);
                    coef_mid[1]->GetEta(Idx(i, 1, k))[iDim] = MidNodeInter4thRight2(inter_temp);
                    coef_mid[1]->GetEta(Idx(i, 2, k))[iDim] = MidNodeInter4th(inter_temp);
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZeta(Idx(i, 1 + iTemp, k))[iDim];
                    }
                    coef_mid[1]->GetZeta(Idx(i, 0, k))[iDim] = MidNodeInter4thRight1(inter_temp);
                    coef_mid[1]->GetZeta(Idx(i, 1, k))[iDim] = MidNodeInter4thRight2(inter_temp);
                    coef_mid[1]->GetZeta(Idx(i, 2, k))[iDim] = MidNodeInter4th(inter_temp);
                    // j=nj-3/2,nj-5/2,nj-7/2
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetXi(Idx(i, nj - 2 - iTemp, k))[iDim];
                    }
                    coef_mid[1]->GetXi(Idx(i, nj - 2, k))[iDim] = MidNodeInter4thLeft1(inter_temp);
                    coef_mid[1]->GetXi(Idx(i, nj - 3, k))[iDim] = MidNodeInter4thLeft2(inter_temp);
                    coef_mid[1]->GetXi(Idx(i, nj - 4, k))[iDim] = MidNodeInter4th(inter_temp);
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetEta(Idx(i, nj - 5 + iTemp, k))[iDim];
                    }
                    coef_mid[1]->GetEta(Idx(i, nj - 2, k))[iDim] = MidNodeInter4thLeft1(inter_temp);
                    coef_mid[1]->GetEta(Idx(i, nj - 3, k))[iDim] = MidNodeInter4thLeft2(inter_temp);
                    coef_mid[1]->GetEta(Idx(i, nj - 4, k))[iDim] = MidNodeInter4th(inter_temp);
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZeta(Idx(i, nj - 5 + iTemp, k))[iDim];
                    }
                    coef_mid[1]->GetZeta(Idx(i, nj - 2, k))[iDim] = MidNodeInter4thLeft1(inter_temp);
                    coef_mid[1]->GetZeta(Idx(i, nj - 3, k))[iDim] = MidNodeInter4thLeft2(inter_temp);
                    coef_mid[1]->GetZeta(Idx(i, nj - 4, k))[iDim] = MidNodeInter4th(inter_temp);
                }
            }
        }
        // k=1/2,3/2,5/2;nk-3/2,nk-5/2,nk-7/2
        for (int j = 1; j < nj - 1; ++j)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                for (int iDim = 0; iDim < 3; ++iDim)
                {
                    // k=1/2，3/2，5/2
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetXi(Idx(i, j, 1 + iTemp))[iDim];
                    }
                    coef_mid[2]->GetXi(Idx(i, j, 0))[iDim] = MidNodeInter4thRight1(inter_temp);
                    coef_mid[2]->GetXi(Idx(i, j, 1))[iDim] = MidNodeInter4thRight2(inter_temp);
                    coef_mid[2]->GetXi(Idx(i, j, 2))[iDim] = MidNodeInter4th(inter_temp);
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetEta(Idx(i, j, 1 + iTemp))[iDim];
                    }
                    coef_mid[2]->GetEta(Idx(i, j, 0))[iDim] = MidNodeInter4thRight1(inter_temp);
                    coef_mid[2]->GetEta(Idx(i, j, 1))[iDim] = MidNodeInter4thRight2(inter_temp);
                    coef_mid[2]->GetEta(Idx(i, j, 2))[iDim] = MidNodeInter4th(inter_temp);
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZeta(Idx(i, j, 1 + iTemp))[iDim];
                    }
                    coef_mid[2]->GetZeta(Idx(i, j, 0))[iDim] = MidNodeInter4thRight1(inter_temp);
                    coef_mid[2]->GetZeta(Idx(i, j, 1))[iDim] = MidNodeInter4thRight2(inter_temp);
                    coef_mid[2]->GetZeta(Idx(i, j, 2))[iDim] = MidNodeInter4th(inter_temp);
                    // k=nk-3/2,nk-5/2,nk-7/2
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetXi(Idx(i, j, nk - 5 + iTemp))[iDim];
                    }
                    coef_mid[2]->GetXi(Idx(i, j, nk - 2))[iDim] = MidNodeInter4thLeft1(inter_temp);
                    coef_mid[2]->GetXi(Idx(i, j, nk - 3))[iDim] = MidNodeInter4thLeft2(inter_temp);
                    coef_mid[2]->GetXi(Idx(i, j, nk - 4))[iDim] = MidNodeInter4th(inter_temp);
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetEta(Idx(i, j, nk - 5 + iTemp))[iDim];
                    }
                    coef_mid[2]->GetEta(Idx(i, j, nk - 2))[iDim] = MidNodeInter4thLeft1(inter_temp);
                    coef_mid[2]->GetEta(Idx(i, j, nk - 3))[iDim] = MidNodeInter4thLeft2(inter_temp);
                    coef_mid[2]->GetEta(Idx(i, j, nk - 4))[iDim] = MidNodeInter4th(inter_temp);
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZeta(Idx(i, j, nk - 5 + iTemp))[iDim];
                    }
                    coef_mid[2]->GetZeta(Idx(i, j, nk - 2))[iDim] = MidNodeInter4thLeft1(inter_temp);
                    coef_mid[2]->GetZeta(Idx(i, j, nk - 3))[iDim] = MidNodeInter4thLeft2(inter_temp);
                    coef_mid[2]->GetZeta(Idx(i, j, nk - 4))[iDim] = MidNodeInter4th(inter_temp);
                }
            }
        }
    }

    void NSSolverStructCMM::CalcInviscidResidual1st()
    {
        CalcMidNode1st();
        CalcMidGhostNodePrimMUSCL();
        CalcInviscidFlux();
        FluxDifference2nd();
    }

    void NSSolverStructCMM::CalcInviscidResidualMUSCL()
    {
        CalcMidNodePrimMUSCL();
        CalcMidGhostNodePrimMUSCL();
        CalcInviscidFlux();
        FluxDifference2nd();
    }

    void NSSolverStructCMM::CalcInviscidResidualWCNS5()
    {
        CalcMidNodePrimWCNS5();
        CalcMidGhostNodePrimWCNS5();
        CalcInviscidFlux();
        FluxDifference6th();
    }

    void NSSolverStructCMM::CalcInviscidFlux()
    {
        auto grid = GetGrid();
        auto data_manager = GetDataManager();
        auto node_metrics = GetNodeMetrics();
        auto idx_proxy = GetIdxProxy();
        auto gas = GetGas();
        auto equ_num = GetPara()->GetEqNum();
        int ni, nj, nk;
        ni = idx_proxy->GetNi();
        nj = idx_proxy->GetNj();
        nk = idx_proxy->GetNk();
        RiemannSolverPara riemann_para;
        riemann_para.gamma_left = riemann_para.gamma_right = gas->GetGamma();
        for (int k = 0; k < nk - 1; ++k)
        {
            for (int j = 0; j < nj - 1; ++j)
            {
                for (int i = 0; i < ni - 1; ++i)
                {
                    int idx = idx_proxy->GetIdx(i, j, k);
                    double jacobi = node_metrics->GetJacobian(idx);
                    // i direction
                    riemann_para.norm[0] = GetMidMetricsI()->GetXi(idx)[0];
                    riemann_para.norm[1] = GetMidMetricsI()->GetXi(idx)[1];
                    riemann_para.norm[2] = GetMidMetricsI()->GetXi(idx)[2];
                    riemann_para.nt = GetMidMetricsI()->GetXi(idx)[3];
                    for (int iEqu = 0; iEqu < equ_num; ++iEqu)
                    {
                        riemann_para.prim_left[iEqu] = data_manager->GetMidNodePrimLeft(iEqu, 0, idx);
                        riemann_para.prim_right[iEqu] = data_manager->GetMidNodePrimRight(iEqu, 0, idx);
                    }
                    m_riemann_solver->Solver(riemann_para);
                    for (int iEqu = 0; iEqu < equ_num; ++iEqu)
                    {
                        data_manager->SetMidNodeFlux(iEqu, 0, idx, riemann_para.flux[iEqu]);
                    }
                    // j direction
                    riemann_para.norm[0] = GetMidMetricsJ()->GetEta(idx)[0];
                    riemann_para.norm[1] = GetMidMetricsJ()->GetEta(idx)[1];
                    riemann_para.norm[2] = GetMidMetricsJ()->GetEta(idx)[2];
                    riemann_para.nt = GetMidMetricsJ()->GetEta(idx)[3];
                    for (int iEqu = 0; iEqu < equ_num; ++iEqu)
                    {
                        riemann_para.prim_left[iEqu] = data_manager->GetMidNodePrimLeft(iEqu, 1, idx);
                        riemann_para.prim_right[iEqu] = data_manager->GetMidNodePrimRight(iEqu, 1, idx);
                    }
                    m_riemann_solver->Solver(riemann_para);
                    for (int iEqu = 0; iEqu < equ_num; ++iEqu)
                    {
                        data_manager->SetMidNodeFlux(iEqu, 1, idx, riemann_para.flux[iEqu]);
                    }
                    // k direction
                    if (grid->GetDim() == 3)
                    {
                        riemann_para.norm[0] = GetMidMetricsK()->GetZeta(idx)[0];
                        riemann_para.norm[1] = GetMidMetricsK()->GetZeta(idx)[1];
                        riemann_para.norm[2] = GetMidMetricsK()->GetZeta(idx)[2];
                        riemann_para.nt = GetMidMetricsK()->GetZeta(idx)[3];
                        for (int iEqu = 0; iEqu < equ_num; ++iEqu)
                        {
                            riemann_para.prim_left[iEqu] = data_manager->GetMidNodePrimLeft(iEqu, 2, idx);
                            riemann_para.prim_right[iEqu] = data_manager->GetMidNodePrimRight(iEqu, 2, idx);
                        }
                        m_riemann_solver->Solver(riemann_para);
                        for (int iEqu = 0; iEqu < equ_num; ++iEqu)
                        {
                            data_manager->SetMidNodeFlux(iEqu, 2, idx, riemann_para.flux[iEqu]);
                        }
                    }
                }
            }
        }
    }
    void NSSolverStructCMM::FluxDifference2nd()
    {
        auto grid = GetGrid();
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
                    idx_temp[0] = idx_proxy->GetIdx(i - 1, j, k);
                    idx_temp[1] = idx_proxy->GetIdx(i, j, k);
                    idx_temp[2] = idx_proxy->GetIdx(i + 1, j, k);
                    for (int iVal = 0; iVal < equ_num; ++iVal)
                    {
                        res_tmp[iVal] -= (data_manager->GetMidNodeFlux(iVal, 0, idx_temp[1]) -
                                          data_manager->GetMidNodeFlux(iVal, 0, idx_temp[0]));
                    }
                    // j direction
                    idx_temp[0] = idx_proxy->GetIdx(i, j - 1, k);
                    idx_temp[1] = idx_proxy->GetIdx(i, j, k);
                    idx_temp[2] = idx_proxy->GetIdx(i, j + 1, k);
                    for (int iVal = 0; iVal < equ_num; ++iVal)
                    {
                        res_tmp[iVal] -= (data_manager->GetMidNodeFlux(iVal, 1, idx_temp[1]) -
                                          data_manager->GetMidNodeFlux(iVal, 1, idx_temp[0]));
                    }
                    // k direction
                    if (grid->GetDim() == 3)
                    {
                        idx_temp[0] = idx_proxy->GetIdx(i, j, k - 1);
                        idx_temp[1] = idx_proxy->GetIdx(i, j, k);
                        idx_temp[2] = idx_proxy->GetIdx(i, j, k + 1);
                        for (int iVal = 0; iVal < equ_num; ++iVal)
                        {
                            res_tmp[iVal] -= (data_manager->GetMidNodeFlux(iVal, 2, idx_temp[1]) -
                                              data_manager->GetMidNodeFlux(iVal, 2, idx_temp[0]));
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
    void NSSolverStructCMM::FluxDifference6th()
    {
        auto grid = GetGrid();
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
        int idx_temp[6];
        // 差分模板的值
        double value[5];
        double res_tmp[5];
        double coef1 = 75.0 / 64.0;
        double coef2 = -25.0 / 384.0;
        double coef3 = 3.0 / 640.0;
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
                    idx_temp[0] = idx_proxy->GetIdx(i - 3, j, k);
                    idx_temp[1] = idx_proxy->GetIdx(i - 2, j, k);
                    idx_temp[2] = idx_proxy->GetIdx(i - 1, j, k);
                    idx_temp[3] = idx_proxy->GetIdx(i, j, k);
                    idx_temp[4] = idx_proxy->GetIdx(i + 1, j, k);
                    idx_temp[5] = idx_proxy->GetIdx(i + 2, j, k);
                    for (int iVal = 0; iVal < equ_num; ++iVal)
                    {
                        res_tmp[iVal] -= coef1 * (data_manager->GetMidNodeFlux(iVal, 0, idx_temp[3]) -
                                                  data_manager->GetMidNodeFlux(iVal, 0, idx_temp[2])) +
                                         coef2 * (data_manager->GetMidNodeFlux(iVal, 0, idx_temp[4]) -
                                                  data_manager->GetMidNodeFlux(iVal, 0, idx_temp[1])) +
                                         coef3 * (data_manager->GetMidNodeFlux(iVal, 0, idx_temp[5]) -
                                                  data_manager->GetMidNodeFlux(iVal, 0, idx_temp[0]));
                    }
                    // j direction
                    idx_temp[0] = idx_proxy->GetIdx(i, j - 3, k);
                    idx_temp[1] = idx_proxy->GetIdx(i, j - 2, k);
                    idx_temp[2] = idx_proxy->GetIdx(i, j - 1, k);
                    idx_temp[3] = idx_proxy->GetIdx(i, j, k);
                    idx_temp[4] = idx_proxy->GetIdx(i, j + 1, k);
                    idx_temp[5] = idx_proxy->GetIdx(i, j + 2, k);
                    for (int iVal = 0; iVal < equ_num; ++iVal)
                    {
                        res_tmp[iVal] -= coef1 * (data_manager->GetMidNodeFlux(iVal, 1, idx_temp[3]) -
                                                  data_manager->GetMidNodeFlux(iVal, 1, idx_temp[2])) +
                                         coef2 * (data_manager->GetMidNodeFlux(iVal, 1, idx_temp[4]) -
                                                  data_manager->GetMidNodeFlux(iVal, 1, idx_temp[1])) +
                                         coef3 * (data_manager->GetMidNodeFlux(iVal, 1, idx_temp[5]) -
                                                  data_manager->GetMidNodeFlux(iVal, 1, idx_temp[0]));
                    }
                    // k direction
                    if (grid->GetDim() == 3)
                    {
                        idx_temp[0] = idx_proxy->GetIdx(i, j, k - 3);
                        idx_temp[1] = idx_proxy->GetIdx(i, j, k - 2);
                        idx_temp[2] = idx_proxy->GetIdx(i, j, k - 1);
                        idx_temp[3] = idx_proxy->GetIdx(i, j, k);
                        idx_temp[4] = idx_proxy->GetIdx(i, j, k + 1);
                        idx_temp[5] = idx_proxy->GetIdx(i, j, k + 2);
                        for (int iVal = 0; iVal < equ_num; ++iVal)
                        {
                            res_tmp[iVal] -= coef1 * (data_manager->GetMidNodeFlux(iVal, 2, idx_temp[3]) -
                                                      data_manager->GetMidNodeFlux(iVal, 2, idx_temp[2])) +
                                             coef2 * (data_manager->GetMidNodeFlux(iVal, 2, idx_temp[4]) -
                                                      data_manager->GetMidNodeFlux(iVal, 2, idx_temp[1])) +
                                             coef3 * (data_manager->GetMidNodeFlux(iVal, 2, idx_temp[5]) -
                                                      data_manager->GetMidNodeFlux(iVal, 2, idx_temp[0]));
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