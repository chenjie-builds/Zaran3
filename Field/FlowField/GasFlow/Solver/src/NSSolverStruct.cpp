#include "NSSolverStruct.h"
#include "MathBasic.h"
#include "FlowSolverStructPara.h"
#include "Log.h"
namespace zaran
{
    NSSolverStruct::NSSolverStruct(int index, string name, FlowSolverPara *para, GridStruct *grid,
                                   DataManagerNSStruct *data_manager)
        : NSSolver(index, name, para, grid, data_manager)
    {
        int ni, nj, nk;
        m_idx_proxy = new StructIdxProxy(grid);
        ni = m_idx_proxy->GetNi();
        nj = m_idx_proxy->GetNj();
        nk = m_idx_proxy->GetNk();
        int node_num = ni * nj * nk;
        m_node_metrics = new Metrics(node_num);
        m_metrics_half_i = new Metrics(node_num);
        m_metrics_half_j = new Metrics(node_num);
        m_metrics_half_k = new Metrics(node_num);
    }
    NSSolverStruct::~NSSolverStruct()
    {
        if (m_idx_proxy)
        {
            delete m_idx_proxy;
            m_idx_proxy = nullptr;
        }
        if (m_node_metrics)
        {
            delete m_node_metrics;
            m_node_metrics = nullptr;
        }
        if (m_metrics_half_i)
        {
            delete m_metrics_half_i;
            m_metrics_half_i = nullptr;
        }
        if (m_metrics_half_j)
        {
            delete m_metrics_half_j;
            m_metrics_half_j = nullptr;
        }
        if (m_metrics_half_k)
        {
            delete m_metrics_half_k;
            m_metrics_half_k = nullptr;
        }
    }
    DataManagerNSStruct *NSSolverStruct::GetDataManager()
    {
        return static_cast<DataManagerNSStruct *>(NSSolver::GetDataManager());
    }
    Metrics *NSSolverStruct::GetNodeMetrics()
    {
        return m_node_metrics;
    }
    StructIdxProxy *NSSolverStruct::GetIdxProxy()
    {
        return m_idx_proxy;
    }
    FlowSolverStructPara *NSSolverStruct::GetPara()
    {
        return static_cast<FlowSolverStructPara *>(NSSolver::GetPara());
    }
    void NSSolverStruct::InitFieldFarFlow()
    {
        auto grid = GetGrid();
        auto data_manager = GetDataManager();
        auto ni = m_idx_proxy->GetNi();
        auto nj = m_idx_proxy->GetNj();
        auto nk = m_idx_proxy->GetNk();
        FlowSolverPara *para = GetPara();
        double prim_far[5];
        prim_far[0] = para->GetInflowDensity();
        prim_far[1] = para->GetInflowVelocityX();
        prim_far[2] = para->GetInflowVelocityY();
        prim_far[3] = para->GetInflowVelocityZ();
        prim_far[4] = para->GetInflowPressure();
        for (int idx_eq = 0; idx_eq < para->GetEqNum(); ++idx_eq)
        {
            for (int k = 0; k < nk; ++k)
            {
                for (int j = 0; j < nj; ++j)
                {
                    for (int i = 0; i < ni; ++i)
                    {
                        int idx = m_idx_proxy->GetIdx(i, j, k);
                        for (int idx_eq = 0; idx_eq < para->GetEqNum(); ++idx_eq)
                        {
                            data_manager->SetPrim(idx_eq, idx, prim_far[idx_eq]);
                        }
                    }
                }
            }
        }
    }
    void NSSolverStruct::InitFieldFarFieldNoVelocity()
    {
        auto grid = GetGrid();
        auto node = grid->GetNode();
        auto data_manager = GetDataManager();
        auto ni = m_idx_proxy->GetNi();
        auto nj = m_idx_proxy->GetNj();
        auto nk = m_idx_proxy->GetNk();
        FlowSolverPara *para = GetPara();
        double prim_far[5];
        prim_far[0] = para->GetInflowDensity();
        prim_far[1] = 0.0;
        prim_far[2] = 0.0;
        prim_far[3] = 0.0;
        prim_far[4] = para->GetInflowPressure();
        for (int k = 0; k < nk; ++k)
        {
            for (int j = 0; j < nj; ++j)
            {
                for (int i = 0; i < ni; ++i)
                {
                    int idx = m_idx_proxy->GetIdx(i, j, k);
                    for (int idx_eq = 0; idx_eq < para->GetEqNum(); ++idx_eq)
                    {
                        data_manager->SetPrim(idx_eq, idx, prim_far[idx_eq]);
                    }
                    auto x = node->GetCoord(i, j, k)[0];
                    if (x < 0.3001)
                    {
                        data_manager->SetPrim(0, idx, 1.0);
                        data_manager->SetPrim(1, idx, 0.0);
                        data_manager->SetPrim(2, idx, 0.0);
                        data_manager->SetPrim(3, idx, 0.0);
                        data_manager->SetPrim(4, idx, 1.0);
                    }
                    else
                    {
                        data_manager->SetPrim(0, idx, 0.125);
                        data_manager->SetPrim(1, idx, 0.0);
                        data_manager->SetPrim(2, idx, 0.0);
                        data_manager->SetPrim(3, idx, 0.0);
                        data_manager->SetPrim(4, idx, 0.1);
                    }
                    // data_manager->SetPrim(0, idx, pow(0.2 * i, 2) * pow(0.2 *
                    // j, 2.3));
                }
            }
        }
    }
    void NSSolverStruct::InitFieldBackup()
    {
    }
    void NSSolverStruct::InitFieldVortex()
    {
        auto grid = GetGrid();
        auto node = grid->GetNode();
        auto data_manager = GetDataManager();
        auto ni = m_idx_proxy->GetNi();
        auto nj = m_idx_proxy->GetNj();
        auto nk = m_idx_proxy->GetNk();
        FlowSolverPara *para = GetPara();
        double beta = 5.0;
        double x, y, z;
        double r2;
        double prim[5];
        double gamma = 1.4;
        for (int k = 0; k < nk; ++k)
        {
            for (int j = 0; j < nj; ++j)
            {
                for (int i = 0; i < ni; ++i)
                {
                    int idx = m_idx_proxy->GetIdx(i, j, k);
                    for (int idx_eq = 0; idx_eq < para->GetEqNum(); ++idx_eq)
                    {
                        data_manager->SetPrim(idx_eq, idx, prim[idx_eq]);
                    }
                    x = node->GetCoord(i, j, k)[0];
                    y = node->GetCoord(i, j, k)[1];
                    z = node->GetCoord(i, j, k)[2];
                    r2 = x * x + y * y + z * z;
                    prim[0] = pow(1.0 - (gamma - 1.0) * beta * beta * exp(1.0 - r2) / (8.0 * gamma * PI * PI), 1.0 / (gamma - 1.0));
                    prim[4] = pow(prim[0], gamma);
                    prim[1] = beta * exp(0.5 * (1.0 - r2)) / (2.0 * PI) * (-y);
                    prim[2] = beta * exp(0.5 * (1.0 - r2)) / (2.0 * PI) * (x);
                    prim[3] = 0.0;
                    for (int idx_eq = 0; idx_eq < para->GetEqNum(); ++idx_eq)
                    {
                        data_manager->SetPrim(idx_eq, idx, prim[idx_eq]);
                    }
                }
            }
        }
    }
    void NSSolverStruct::CalcMetricsS0()
    {
        auto para = GetPara();
        auto flux_diff_scheme = para->GetDifferenceScheme();
        if (flux_diff_scheme == DifferenceScheme::SecondOrder)
        {
            CalcMetricsS0_2nd();
        }
        else if (flux_diff_scheme == DifferenceScheme::SixthOrder)
        {
            CalcMetricsS0_6th();
        }
        else
        {
            Log::warn("FluxDifferenceScheme is not defined, use SecondOrder as default");
            CalcMetricsS0_2nd();
        }
    }

    void NSSolverStruct::CalcMetricsS0_2nd()
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
        // i+1/2,j+1/2,k+1/2处的坐标
        std::vector<std::vector<double>> coord_i(ni * nj * nk, std::vector<double>(3)),
            coord_j(ni * nj * nk, std::vector<double>(3)), coord_k(ni * nj * nk, std::vector<double>(3));
        // 计算之前，先把度量系数赋值为0
        for (int i = 0; i < ni; ++i)
        {
            for (int j = 0; j < nj; ++j)
            {
                for (int k = 0; k < nk; ++k)
                {
                    for (int iDim = 0; iDim < 4; ++iDim)
                    {
                        coef->GetX(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetY(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetZ(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetXi(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetEta(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetZeta(Idx(i, j, k))[iDim] = 0.0;
                        for (int jDim = 0; jDim < 3; ++jDim)
                        {
                            coef_mid[jDim]->GetX(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetY(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetZ(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetXi(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetEta(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetZeta(Idx(i, j, k))[iDim] = 0.0;
                        }
                    }
                }
            }
        }
        // 第一步：计算半点坐标
        for (int k = 0; k < nk - 1; ++k)
        {
            for (int j = 0; j < nj - 1; ++j)
            {
                for (int i = 0; i < ni - 1; ++i)
                {
                    for (int iDim = 0; iDim < 3; ++iDim)
                    {
                        int idx = idx_proxy->GetIdx(i, j, k);
                        coord_i[idx][iDim] = 0.5 * (node->GetCoord(i, j, k)[iDim] + node->GetCoord(i + 1, j, k)[iDim]);
                        coord_j[idx][iDim] = 0.5 * (node->GetCoord(i, j, k)[iDim] + node->GetCoord(i, j + 1, k)[iDim]);
                        coord_k[idx][iDim] = 0.5 * (node->GetCoord(i, j, k)[iDim] + node->GetCoord(i, j, k + 1)[iDim]);
                    }
                }
            }
        }
        // 第二步：根据半点坐标计算整点度量系数
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    int idx = idx_proxy->GetIdx(i, j, k);
                    auto coef_x = coef->GetX(idx);
                    auto coef_y = coef->GetY(idx);
                    auto coef_z = coef->GetZ(idx);
                    coef_x[0] = coord_i[Idx(i, j, k)][0] - coord_i[Idx(i - 1, j, k)][0];
                    coef_x[1] = coord_j[Idx(i, j, k)][0] - coord_j[Idx(i, j - 1, k)][0];
                    coef_x[2] = coord_k[Idx(i, j, k)][0] - coord_k[Idx(i, j, k - 1)][0];
                    coef_y[0] = coord_i[Idx(i, j, k)][1] - coord_i[Idx(i - 1, j, k)][1];
                    coef_y[1] = coord_j[Idx(i, j, k)][1] - coord_j[Idx(i, j - 1, k)][1];
                    coef_y[2] = coord_k[Idx(i, j, k)][1] - coord_k[Idx(i, j, k - 1)][1];
                    coef_z[0] = coord_i[Idx(i, j, k)][2] - coord_i[Idx(i - 1, j, k)][2];
                    coef_z[1] = coord_j[Idx(i, j, k)][2] - coord_j[Idx(i, j - 1, k)][2];
                    coef_z[2] = coord_k[Idx(i, j, k)][2] - coord_k[Idx(i, j, k - 1)][2];
                    coef_x[3] = 0.0;
                    coef_y[3] = 0.0;
                    if (grid->GetDim() == 2)
                    {
                        coef_x[2] = 0.0;
                        coef_y[2] = 0.0;
                        coef_z[0] = 0.0;
                        coef_z[1] = 0.0;
                        coef_z[2] = 1.0;
                    }
                    // 计算xi,eta,zeta
                    auto coef_xi = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, k));
                    coef_xi[0] = coef_y[1] * coef_z[2] - coef_z[1] * coef_y[2];
                    coef_xi[1] = coef_z[1] * coef_x[2] - coef_x[1] * coef_z[2];
                    coef_xi[2] = coef_x[1] * coef_y[2] - coef_y[1] * coef_x[2];
                    coef_xi[3] = 0.0;
                    auto coef_eta = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, k));
                    coef_eta[0] = coef_y[2] * coef_z[0] - coef_z[2] * coef_y[0];
                    coef_eta[1] = coef_z[2] * coef_x[0] - coef_x[2] * coef_z[0];
                    coef_eta[2] = coef_x[2] * coef_y[0] - coef_y[2] * coef_x[0];
                    coef_eta[3] = 0.0;
                    auto coef_zeta = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, k));
                    coef_zeta[0] = coef_y[0] * coef_z[1] - coef_z[0] * coef_y[1];
                    coef_zeta[1] = coef_z[0] * coef_x[1] - coef_x[0] * coef_z[1];
                    coef_zeta[2] = coef_x[0] * coef_y[1] - coef_y[0] * coef_x[1];
                    coef_zeta[3] = 0.0;
                }
            }
        }
    }

    void NSSolverStruct::CalcMetricsS0_6th()
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
        double inter_temp[6];
        // 用于度量系数计算的临时变量
        std::vector<std::vector<double>> coord_i(ni * nj * nk, std::vector<double>(3)),
            coord_j(ni * nj * nk, std::vector<double>(3)), coord_k(ni * nj * nk, std::vector<double>(3));
        // 计算之前，先把度量系数赋值为0
        for (int i = 0; i < ni; ++i)
        {
            for (int j = 0; j < nj; ++j)
            {
                for (int k = 0; k < nk; ++k)
                {
                    for (int iDim = 0; iDim < 4; ++iDim)
                    {
                        coef->GetX(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetY(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetZ(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetXi(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetEta(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetZeta(Idx(i, j, k))[iDim] = 0.0;
                        for (int jDim = 0; jDim < 3; ++jDim)
                        {
                            coef_mid[jDim]->GetX(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetY(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetZ(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetXi(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetEta(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetZeta(Idx(i, j, k))[iDim] = 0.0;
                        }
                    }
                }
            }
        }
        // 第一步：计算半点坐标
        // i+1/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 3; i < ni - 4; ++i)
                {
                    for (int iDim = 0; iDim < 3; ++iDim)
                    {
                        for (int iTemp = 0; iTemp < 6; iTemp++)
                        {
                            inter_temp[iTemp] = node->GetCoord(i + iTemp - 2, j, k)[iDim];
                        }
                        coord_i[Idx(i, j, k)][iDim] = MidNodeInter6th(inter_temp);
                    }
                }
            }
        }
        // i=1/2,3/2,5/2;ni-3/2,ni-5/2,ni-7/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int iDim = 0; iDim < 3; ++iDim)
                {
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(1 + iTemp, j, k)[iDim];
                    }
                    coord_i[Idx(0, j, k)][iDim] = MidNodeInter4thRight1(inter_temp);
                    coord_i[Idx(1, j, k)][iDim] = MidNodeInter4thRight2(inter_temp);
                    coord_i[Idx(2, j, k)][iDim] = MidNodeInter4th(inter_temp);
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(ni - 5 + iTemp, j, k)[iDim];
                    }
                    coord_i[Idx(ni - 2, j, k)][iDim] = MidNodeInter4thLeft1(inter_temp);
                    coord_i[Idx(ni - 3, j, k)][iDim] = MidNodeInter4thLeft2(inter_temp);
                    coord_i[Idx(ni - 4, j, k)][iDim] = MidNodeInter4th(inter_temp);
                }
            }
        }
        // j+1/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 3; j < nj - 4; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    for (int iDim = 0; iDim < 3; ++iDim)
                    {
                        for (int iTemp = 0; iTemp < 6; iTemp++)
                        {
                            inter_temp[iTemp] = node->GetCoord(i, j + iTemp - 2, k)[iDim];
                        }
                        coord_j[Idx(i, j, k)][iDim] = MidNodeInter6th(inter_temp);
                    }
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
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(i, 1 + iTemp, k)[iDim];
                    }
                    coord_j[Idx(i, 0, k)][iDim] = MidNodeInter4thRight1(inter_temp);
                    coord_j[Idx(i, 1, k)][iDim] = MidNodeInter4thRight2(inter_temp);
                    coord_j[Idx(i, 2, k)][iDim] = MidNodeInter4th(inter_temp);
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(i, nj - 5 + iTemp, k)[iDim];
                    }
                    coord_j[Idx(i, nj - 2, k)][iDim] = MidNodeInter4thLeft1(inter_temp);
                    coord_j[Idx(i, nj - 3, k)][iDim] = MidNodeInter4thLeft2(inter_temp);
                    coord_j[Idx(i, nj - 4, k)][iDim] = MidNodeInter4th(inter_temp);
                }
            }
        }
        // k+1/2
        for (int k = 3; k < nk - 4; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    for (int iDim = 0; iDim < 3; ++iDim)
                    {
                        for (int iTemp = 0; iTemp < 6; iTemp++)
                        {
                            inter_temp[iTemp] = node->GetCoord(i, j, k + iTemp - 2)[iDim];
                        }
                        coord_k[Idx(i, j, k)][iDim] = MidNodeInter6th(inter_temp);
                    }
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
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(i, j, 1 + iTemp)[iDim];
                    }
                    coord_k[Idx(i, j, 0)][iDim] = MidNodeInter4thRight1(inter_temp);
                    coord_k[Idx(i, j, 1)][iDim] = MidNodeInter4thRight2(inter_temp);
                    coord_k[Idx(i, j, 2)][iDim] = MidNodeInter4th(inter_temp);
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(i, j, nk - 5 + iTemp)[iDim];
                    }
                    coord_k[Idx(i, j, nk - 2)][iDim] = MidNodeInter4thLeft1(inter_temp);
                    coord_k[Idx(i, j, nk - 3)][iDim] = MidNodeInter4thLeft2(inter_temp);
                    coord_k[Idx(i, j, nk - 4)][iDim] = MidNodeInter4th(inter_temp);
                }
            }
        }
        // 第二步：根据半点坐标计算整点度量系数
        double diff_temp[6];
        // i direction
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 3; i < ni - 3; ++i)
                {
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_i[Idx(i + iTemp - 3, j, k)][0];
                    }
                    auto coef_x = coef->GetX(Idx(i, j, k));
                    coef_x[0] = NodeDifferece6th(diff_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_i[Idx(i + iTemp - 3, j, k)][1];
                    }
                    auto coef_y = coef->GetY(Idx(i, j, k));
                    coef_y[0] = NodeDifferece6th(diff_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_i[Idx(i + iTemp - 3, j, k)][2];
                    }
                    auto coef_z = coef->GetZ(Idx(i, j, k));
                    coef_z[0] = NodeDifferece6th(diff_temp);
                    if (grid->GetDim() == 2)
                    {
                        coef_z[0] = 0.0;
                    }
                }
            }
        }
        // j direction
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 3; j < nj - 3; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    auto coef_x = coef->GetX(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_j[Idx(i, j + iTemp - 3, k)][0];
                    }
                    coef_x[1] = NodeDifferece6th(diff_temp);
                    auto coef_y = coef->GetY(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_j[Idx(i, j + iTemp - 3, k)][1];
                    }
                    coef_y[1] = NodeDifferece6th(diff_temp);
                    auto coef_z = coef->GetZ(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_j[Idx(i, j + iTemp - 3, k)][2];
                    }
                    coef_z[1] = NodeDifferece6th(diff_temp);
                    if (grid->GetDim() == 2)
                    {
                        coef_z[1] = 0.0;
                    }
                }
            }
        }
        // k direction
        for (int k = 3; k < nk - 3; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    auto coef_x = coef->GetX(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_k[Idx(i, j, k + iTemp - 3)][0];
                    }
                    coef_x[2] = NodeDifferece6th(diff_temp);
                    auto coef_y = coef->GetY(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_k[Idx(i, j, k + iTemp - 3)][1];
                    }
                    coef_y[2] = NodeDifferece6th(diff_temp);
                    auto coef_z = coef->GetZ(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_k[Idx(i, j, k + iTemp - 3)][2];
                    }
                    coef_z[2] = NodeDifferece6th(diff_temp);
                    if (grid->GetDim() == 2)
                    {
                        coef_x[2] = 0.0;
                        coef_y[2] = 0.0;
                        coef_z[2] = 1.0;
                    }
                }
            }
        }
        // 计算边界整点的度量系数
        // i=1,2;ni-3,ni-2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_i[Idx(0 + iTemp, j, k)][0];
                }
                coef->GetX(Idx(1, j, k))[0] = NodeDifferece4thRight(diff_temp);
                coef->GetX(Idx(2, j, k))[0] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_i[Idx(0 + iTemp, j, k)][1];
                }
                coef->GetY(Idx(1, j, k))[0] = NodeDifferece4thRight(diff_temp);
                coef->GetY(Idx(2, j, k))[0] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_i[Idx(0 + iTemp, j, k)][2];
                }
                coef->GetZ(Idx(1, j, k))[0] = NodeDifferece4thRight(diff_temp);
                coef->GetZ(Idx(2, j, k))[0] = NodeDifferece4th(diff_temp);

                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_i[Idx(ni - 2 - iTemp, j, k)][0];
                }
                coef->GetX(Idx(ni - 2, j, k))[0] = NodeDifferece4thLeft(diff_temp);
                coef->GetX(Idx(ni - 3, j, k))[0] = -NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_i[Idx(ni - 2 - iTemp, j, k)][1];
                }
                coef->GetY(Idx(ni - 2, j, k))[0] = NodeDifferece4thLeft(diff_temp);
                coef->GetY(Idx(ni - 3, j, k))[0] = -NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_i[Idx(ni - 2 - iTemp, j, k)][2];
                }
                coef->GetZ(Idx(ni - 2, j, k))[0] = NodeDifferece4thLeft(diff_temp);
                coef->GetZ(Idx(ni - 3, j, k))[0] = -NodeDifferece4th(diff_temp);
            }
        }
        // j=1,2;nj-3,nj-2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_j[Idx(i, 0 + iTemp, k)][0];
                }
                coef->GetX(Idx(i, 1, k))[1] = NodeDifferece4thRight(diff_temp);
                coef->GetX(Idx(i, 2, k))[1] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_j[Idx(i, 0 + iTemp, k)][1];
                }
                coef->GetY(Idx(i, 1, k))[1] = NodeDifferece4thRight(diff_temp);
                coef->GetY(Idx(i, 2, k))[1] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_j[Idx(i, 0 + iTemp, k)][2];
                }
                coef->GetZ(Idx(i, 1, k))[1] = NodeDifferece4thRight(diff_temp);
                coef->GetZ(Idx(i, 2, k))[1] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_j[Idx(i, nj - 2 - iTemp, k)][0];
                }
                coef->GetX(Idx(i, nj - 2, k))[1] = NodeDifferece4thLeft(diff_temp);
                coef->GetX(Idx(i, nj - 3, k))[1] = -NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_j[Idx(i, nj - 2 - iTemp, k)][1];
                }
                coef->GetY(Idx(i, nj - 2, k))[1] = NodeDifferece4thLeft(diff_temp);
                coef->GetY(Idx(i, nj - 3, k))[1] = -NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_j[Idx(i, nj - 2 - iTemp, k)][2];
                }
                coef->GetZ(Idx(i, nj - 2, k))[1] = NodeDifferece4thLeft(diff_temp);
                coef->GetZ(Idx(i, nj - 3, k))[1] = -NodeDifferece4th(diff_temp);
            }
        }
        // k=1,2;nk-3,nk-2
        for (int j = 1; j < nj - 1; ++j)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_k[Idx(i, j, 0 + iTemp)][0];
                }
                coef->GetX(Idx(i, j, 1))[2] = NodeDifferece4thRight(diff_temp);
                coef->GetX(Idx(i, j, 2))[2] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_k[Idx(i, j, 0 + iTemp)][1];
                }
                coef->GetY(Idx(i, j, 1))[2] = NodeDifferece4thRight(diff_temp);
                coef->GetY(Idx(i, j, 2))[2] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_k[Idx(i, j, 0 + iTemp)][2];
                }
                coef->GetZ(Idx(i, j, 1))[2] = NodeDifferece4thRight(diff_temp);
                coef->GetZ(Idx(i, j, 2))[2] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_k[Idx(i, j, nk - 2 - iTemp)][0];
                }
                coef->GetX(Idx(i, j, nk - 2))[2] = NodeDifferece4thLeft(diff_temp);
                coef->GetX(Idx(i, j, nk - 3))[2] = -NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_k[Idx(i, j, nk - 2 - iTemp)][1];
                }
                coef->GetY(Idx(i, j, nk - 2))[2] = NodeDifferece4thLeft(diff_temp);
                coef->GetY(Idx(i, j, nk - 3))[2] = -NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_k[Idx(i, j, nk - 2 - iTemp)][2];
                }
                coef->GetZ(Idx(i, j, nk - 2))[2] = NodeDifferece4thLeft(diff_temp);
                coef->GetZ(Idx(i, j, nk - 3))[2] = -NodeDifferece4th(diff_temp);
                if (grid->GetDim() == 2)
                {
                    coef->GetX(Idx(i, j, 1))[2] = 0.0;
                    coef->GetX(Idx(i, j, 2))[2] = 0.0;
                    coef->GetX(Idx(i, j, nk - 2))[2] = 0.0;
                    coef->GetX(Idx(i, j, nk - 3))[2] = 0.0;
                    coef->GetY(Idx(i, j, 1))[2] = 0.0;
                    coef->GetY(Idx(i, j, 2))[2] = 0.0;
                    coef->GetY(Idx(i, j, nk - 2))[2] = 0.0;
                    coef->GetY(Idx(i, j, nk - 3))[2] = 0.0;
                    coef->GetZ(Idx(i, j, 1))[0] = coef->GetZ(Idx(i, j, 1))[1] = 0.0;
                    coef->GetZ(Idx(i, j, 2))[0] = coef->GetZ(Idx(i, j, 2))[1] = 0.0;
                    coef->GetZ(Idx(i, j, nk - 2))[0] = coef->GetZ(Idx(i, j, nk - 2))[1] = 0.0;
                    coef->GetZ(Idx(i, j, nk - 3))[0] = coef->GetZ(Idx(i, j, nk - 3))[1] = 0.0;
                    coef->GetZ(Idx(i, j, 1))[2] = 1.0;
                    coef->GetZ(Idx(i, j, 2))[2] = 1.0;
                    coef->GetZ(Idx(i, j, nk - 2))[2] = 1.0;
                    coef->GetZ(Idx(i, j, nk - 3))[2] = 1.0;
                }
            }
        }
        // 计算整点的度量系数
        for (int k = 1; k <= nk - 2; ++k)
        {
            for (int j = 1; j <= nj - 2; ++j)
            {
                for (int i = 1; i <= ni - 2; ++i)
                {
                    // 计算xi,eta,zeta
                    auto coef_x = GetNodeMetrics()->GetX(idx_proxy->GetIdx(i, j, k));
                    auto coef_y = GetNodeMetrics()->GetY(idx_proxy->GetIdx(i, j, k));
                    auto coef_z = GetNodeMetrics()->GetZ(idx_proxy->GetIdx(i, j, k));
                    auto coef_xi = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, k));
                    coef_xi[0] = coef_y[1] * coef_z[2] - coef_z[1] * coef_y[2];
                    coef_xi[1] = coef_z[1] * coef_x[2] - coef_x[1] * coef_z[2];
                    coef_xi[2] = coef_x[1] * coef_y[2] - coef_y[1] * coef_x[2];
                    coef_xi[3] = 0.0;
                    auto coef_eta = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, k));
                    coef_eta[0] = coef_y[2] * coef_z[0] - coef_z[2] * coef_y[0];
                    coef_eta[1] = coef_z[2] * coef_x[0] - coef_x[2] * coef_z[0];
                    coef_eta[2] = coef_x[2] * coef_y[0] - coef_y[2] * coef_x[0];
                    coef_eta[3] = 0.0;
                    auto coef_zeta = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, k));
                    coef_zeta[0] = coef_y[0] * coef_z[1] - coef_z[0] * coef_y[1];
                    coef_zeta[1] = coef_z[0] * coef_x[1] - coef_x[0] * coef_z[1];
                    coef_zeta[2] = coef_x[0] * coef_y[1] - coef_y[0] * coef_x[1];
                    coef_zeta[3] = 0.0;
                }
            }
        }
    }

    void NSSolverStruct::CalcCoordTransCoef()
    {
        CalcMetrics();
        CalcJacobian();
    }
    void NSSolverStruct::CalcMetrics()
    {
        if (GetGrid()->GetDim() == 2)
        {
            Log::warn("2D case, use Original as default");
            CalcMetricsS0();
            return;
        }
        auto para = GetPara();
        auto metrics_type = para->GetMetricsScheme();
        if (metrics_type == MetricsScheme::S0)
        {
            CalcMetricsS0();
        }
        else if (metrics_type == MetricsScheme::S1)
        {
            CalcMetricsS1();
        }
        else if (metrics_type == MetricsScheme::S2)
        {
            CalcMetricsS2();
        }
        else if (metrics_type == MetricsScheme::S3)
        {
            CalcMetricsS3();
        }
        else
        {
            Log::warn("MetricsType is not defined, use Originnal as default");
            CalcMetricsS0();
        }
    }
    void NSSolverStruct::CalcJacobian()
    {
        auto para = GetPara();
        auto jacobian_scheme = para->GetJacobianScheme();
        if (jacobian_scheme == JacobianScheme::V1)
        {
            CalcJacobianV1();
        }
        else if (jacobian_scheme == JacobianScheme::V2)
        {
            CalcJacobianV2();
        }
        else if (jacobian_scheme == JacobianScheme::V3)
        {
            CalcJacobianV3();
        }
        else
        {
            Log::warn("JacobianScheme is not defined, use V1 as default");
            CalcJacobianV1();
        }
    }
    void NSSolverStruct::CalcJacobianV1()
    {
        auto grid = GetGrid();
        auto node = grid->GetNode();
        int ni = grid->GetNi();
        int nj = grid->GetNj();
        int nk = grid->GetNk();
        auto idx_proxy = GetIdxProxy();
        auto node_metrics = GetNodeMetrics();
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    auto idx = idx_proxy->GetIdx(i, j, k);
                    auto &jacobian = node_metrics->GetJacobian(idx);
                    auto coef_x = node_metrics->GetX(idx);
                    auto coef_y = node_metrics->GetY(idx);
                    auto coef_z = node_metrics->GetZ(idx);
                    jacobian = coef_x[0] * coef_y[1] * coef_z[2] - coef_x[0] * coef_y[2] * coef_z[1] +
                               coef_x[1] * coef_y[2] * coef_z[0] - coef_x[1] * coef_y[0] * coef_z[2] +
                               coef_x[2] * coef_y[0] * coef_z[1] - coef_x[2] * coef_y[1] * coef_z[0];
                    jacobian = 1.0 / jacobian;
                }
            }
        }
    }
    void NSSolverStruct::CalcJacobianV2()
    {
        auto grid = GetGrid();
        auto node = grid->GetNode();
        int ni = grid->GetNi();
        int nj = grid->GetNj();
        int nk = grid->GetNk();
        auto idx_proxy = GetIdxProxy();
        auto node_metrics = GetNodeMetrics();
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    auto idx = idx_proxy->GetIdx(i, j, k);
                    auto &jacobian = node_metrics->GetJacobian(idx);
                    auto coef_x = node_metrics->GetX(idx);
                    auto coef_y = node_metrics->GetY(idx);
                    auto coef_z = node_metrics->GetZ(idx);
                    auto coef_xi = node_metrics->GetXi(idx);
                    auto coef_eta = node_metrics->GetEta(idx);
                    auto coef_zeta = node_metrics->GetZeta(idx);
                    jacobian = coef_x[0] * coef_xi[0] + coef_y[0] * coef_xi[1] + coef_z[0] * coef_xi[2] +
                               coef_x[1] * coef_eta[0] + coef_y[1] * coef_eta[1] + coef_z[1] * coef_eta[2] +
                               coef_x[2] * coef_zeta[0] + coef_y[2] * coef_zeta[1] + coef_z[2] * coef_zeta[2];
                    jacobian = jacobian / 3.0;
                    jacobian = 1.0 / jacobian;
                }
            }
        }
    }

    void NSSolverStruct::CalcJacobianV3()
    {
        auto para = GetPara();
        auto flux_diff_scheme = para->GetDifferenceScheme();
        if (flux_diff_scheme == DifferenceScheme::SecondOrder)
        {
            CalcJacobianV3_2nd();
        }
        else if (flux_diff_scheme == DifferenceScheme::SixthOrder)
        {
            CalcJacobianV3_6th();
        }
        else
        {
            Log::warn("FluxDifferenceScheme is not defined, use SecondOrder as default");
            CalcJacobianV3_2nd();
        }
    }
    void NSSolverStruct::CalcJacobianV3_2nd()
    {
        auto grid = GetGrid();
        auto node = grid->GetNode();
        auto coef = GetNodeMetrics();
        auto idx_proxy = GetIdxProxy();
        auto Idx = [&](int i, int j, int k)
        {
            return idx_proxy->GetIdx(i, j, k);
        };
        int ni = grid->GetNi();
        int nj = grid->GetNj();
        int nk = grid->GetNk();
        double inter_temp[2];
        // i+1/2,j+1/2,k+1/2处的临时变量
        std::vector<double> temp_i(ni * nj * nk, 0.0), temp_j(ni * nj * nk, 0.0), temp_k(ni * nj * nk, 0.0);
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    int idx0 = Idx(i, j, k);
                    int idx1 = Idx(i + 1, j, k);
                    auto coord0 = node->GetCoord(i, j, k);
                    auto coord1 = node->GetCoord(i + 1, j, k);
                    inter_temp[0] = coord0[0] * coef->GetXi(idx0)[0] + coord0[1] * coef->GetXi(idx0)[1] + coord0[2] * coef->GetXi(idx0)[2];
                    inter_temp[1] = coord1[0] * coef->GetXi(idx1)[0] + coord1[1] * coef->GetXi(idx1)[1] + coord1[2] * coef->GetXi(idx1)[2];
                    temp_i[idx0] = 0.5 * (inter_temp[0] + inter_temp[1]);
                    idx0 = Idx(i, j, k);
                    idx1 = Idx(i, j + 1, k);
                    coord0 = node->GetCoord(i, j, k);
                    coord1 = node->GetCoord(i, j + 1, k);
                    inter_temp[0] = coord0[0] * coef->GetEta(idx0)[0] + coord0[1] * coef->GetEta(idx0)[1] + coord0[2] * coef->GetEta(idx0)[2];
                    inter_temp[1] = coord1[0] * coef->GetEta(idx1)[0] + coord1[1] * coef->GetEta(idx1)[1] + coord1[2] * coef->GetEta(idx1)[2];
                    temp_j[idx0] = 0.5 * (inter_temp[0] + inter_temp[1]);
                    idx0 = Idx(i, j, k);
                    idx1 = Idx(i, j, k + 1);
                    coord0 = node->GetCoord(i, j, k);
                    coord1 = node->GetCoord(i, j, k + 1);
                    inter_temp[0] = coord0[0] * coef->GetZeta(idx0)[0] + coord0[1] * coef->GetZeta(idx0)[1] + coord0[2] * coef->GetZeta(idx0)[2];
                    inter_temp[1] = coord1[0] * coef->GetZeta(idx1)[0] + coord1[1] * coef->GetZeta(idx1)[1] + coord1[2] * coef->GetZeta(idx1)[2];
                    temp_k[idx0] = 0.5 * (inter_temp[0] + inter_temp[1]);
                }
            }
        }
        // i=1/2,(ni-1)-1/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                int idx0 = Idx(1, j, k);
                int idx1 = Idx(2, j, k);
                auto coord0 = node->GetCoord(1, j, k);
                auto coord1 = node->GetCoord(2, j, k);
                inter_temp[0] = coord0[0] * coef->GetXi(idx0)[0] + coord0[1] * coef->GetXi(idx0)[1] + coord0[2] * coef->GetXi(idx0)[2];
                inter_temp[1] = coord1[0] * coef->GetXi(idx1)[0] + coord1[1] * coef->GetXi(idx1)[1] + coord1[2] * coef->GetXi(idx1)[2];
                temp_i[idx0] = 0.5 * (3.0 * inter_temp[0] - 1.0 * inter_temp[1]);
                idx0 = Idx(ni - 2, j, k);
                idx1 = Idx(ni - 3, j, k);
                coord0 = node->GetCoord(ni - 2, j, k);
                coord1 = node->GetCoord(ni - 3, j, k);
                inter_temp[0] = coord0[0] * coef->GetXi(idx0)[0] + coord0[1] * coef->GetXi(idx0)[1] + coord0[2] * coef->GetXi(idx0)[2];
                inter_temp[1] = coord1[0] * coef->GetXi(idx1)[0] + coord1[1] * coef->GetXi(idx1)[1] + coord1[2] * coef->GetXi(idx1)[2];
                temp_i[idx0] = 0.5 * (3.0 * inter_temp[0] - 1.0 * inter_temp[1]);
            }
        }
        // j=1/2,(nj-1)-1/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                int idx0 = Idx(i, 1, k);
                int idx1 = Idx(i, 2, k);
                auto coord0 = node->GetCoord(i, 1, k);
                auto coord1 = node->GetCoord(i, 2, k);
                inter_temp[0] = coord0[0] * coef->GetEta(idx0)[0] + coord0[1] * coef->GetEta(idx0)[1] + coord0[2] * coef->GetEta(idx0)[2];
                inter_temp[1] = coord1[0] * coef->GetEta(idx1)[0] + coord1[1] * coef->GetEta(idx1)[1] + coord1[2] * coef->GetEta(idx1)[2];
                temp_j[idx0] = 0.5 * (3.0 * inter_temp[0] - 1.0 * inter_temp[1]);
                idx0 = Idx(i, nj - 2, k);
                idx1 = Idx(i, nj - 3, k);
                coord0 = node->GetCoord(i, nj - 2, k);
                coord1 = node->GetCoord(i, nj - 3, k);
                inter_temp[0] = coord0[0] * coef->GetEta(idx0)[0] + coord0[1] * coef->GetEta(idx0)[1] + coord0[2] * coef->GetEta(idx0)[2];
                inter_temp[1] = coord1[0] * coef->GetEta(idx1)[0] + coord1[1] * coef->GetEta(idx1)[1] + coord1[2] * coef->GetEta(idx1)[2];
                temp_j[idx0] = 0.5 * (3.0 * inter_temp[0] - 1.0 * inter_temp[1]);
            }
        }
        // k=1/2,(nk-1)-1/2
        for (int j = 1; j < nj - 1; ++j)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                int idx0 = Idx(i, j, 1);
                int idx1 = Idx(i, j, 2);
                auto coord0 = node->GetCoord(i, j, 1);
                auto coord1 = node->GetCoord(i, j, 2);
                inter_temp[0] = coord0[0] * coef->GetZeta(idx0)[0] + coord0[1] * coef->GetZeta(idx0)[1] + coord0[2] * coef->GetZeta(idx0)[2];
                inter_temp[1] = coord1[0] * coef->GetZeta(idx1)[0] + coord1[1] * coef->GetZeta(idx1)[1] + coord1[2] * coef->GetZeta(idx1)[2];
                temp_k[idx0] = 0.5 * (3.0 * inter_temp[0] - 1.0 * inter_temp[1]);
                idx0 = Idx(i, j, nk - 2);
                idx1 = Idx(i, j, nk - 3);
                coord0 = node->GetCoord(i, j, nk - 2);
                coord1 = node->GetCoord(i, j, nk - 3);
                inter_temp[0] = coord0[0] * coef->GetZeta(idx0)[0] + coord0[1] * coef->GetZeta(idx0)[1] + coord0[2] * coef->GetZeta(idx0)[2];
                inter_temp[1] = coord1[0] * coef->GetZeta(idx1)[0] + coord1[1] * coef->GetZeta(idx1)[1] + coord1[2] * coef->GetZeta(idx1)[2];
                temp_k[idx0] = 0.5 * (3.0 * inter_temp[0] - 1.0 * inter_temp[1]);
            }
        }
        // 计算整点的Jacobian
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    int idx = Idx(i, j, k);
                    auto &jacobian = coef->GetJacobian(idx);
                    jacobian = temp_i[Idx(i, j, k)] - temp_i[Idx(i - 1, j, k)] +
                               temp_j[Idx(i, j, k)] - temp_j[Idx(i, j - 1, k)] +
                               temp_k[Idx(i, j, k)] - temp_k[Idx(i, j, k - 1)];
                    jacobian = jacobian / 3.0;
                    jacobian = 1.0 / jacobian;
                }
            }
        }
    }
    void NSSolverStruct::CalcJacobianV3_6th()
    {
        auto grid = GetGrid();
        auto node = grid->GetNode();
        auto coef = GetNodeMetrics();
        auto idx_proxy = GetIdxProxy();
        auto Idx = [&](int i, int j, int k)
        {
            return idx_proxy->GetIdx(i, j, k);
        };
        int ni = grid->GetNi();
        int nj = grid->GetNj();
        int nk = grid->GetNk();
        double inter_temp[6];
        // i+1/2,j+1/2,k+1/2处的临时变量
        std::vector<double> temp_i(ni * nj * nk, 0.0), temp_j(ni * nj * nk, 0.0), temp_k(ni * nj * nk, 0.0);
        // i+1/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 3; i < ni - 4; ++i)
                {
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(i + iTemp - 2, j, k)[0] * coef->GetXi(Idx(i + iTemp - 2, j, k))[0] +
                                            node->GetCoord(i + iTemp - 2, j, k)[1] * coef->GetXi(Idx(i + iTemp - 2, j, k))[1] +
                                            node->GetCoord(i + iTemp - 2, j, k)[2] * coef->GetXi(Idx(i + iTemp - 2, j, k))[2];
                    }
                    temp_i[Idx(i, j, k)] = MidNodeInter6th(inter_temp);
                }
            }
        }
        // i=1/2,3/2,5/2;ni-3/2,ni-5/2,ni-7/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = node->GetCoord(1 + iTemp, j, k)[0] * coef->GetXi(Idx(1 + iTemp, j, k))[0] +
                                        node->GetCoord(1 + iTemp, j, k)[1] * coef->GetXi(Idx(1 + iTemp, j, k))[1] +
                                        node->GetCoord(1 + iTemp, j, k)[2] * coef->GetXi(Idx(1 + iTemp, j, k))[2];
                }
                temp_i[Idx(0, j, k)] = MidNodeInter4thRight1(inter_temp);
                temp_i[Idx(1, j, k)] = MidNodeInter4thRight2(inter_temp);
                temp_i[Idx(2, j, k)] = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = node->GetCoord(ni - 5 + iTemp, j, k)[0] * coef->GetXi(Idx(ni - 5 + iTemp, j, k))[0] +
                                        node->GetCoord(ni - 5 + iTemp, j, k)[1] * coef->GetXi(Idx(ni - 5 + iTemp, j, k))[1] +
                                        node->GetCoord(ni - 5 + iTemp, j, k)[2] * coef->GetXi(Idx(ni - 5 + iTemp, j, k))[2];
                }
                temp_i[Idx(ni - 2, j, k)] = MidNodeInter4thLeft1(inter_temp);
                temp_i[Idx(ni - 3, j, k)] = MidNodeInter4thLeft2(inter_temp);
                temp_i[Idx(ni - 4, j, k)] = MidNodeInter4th(inter_temp);
            }
        }
        // j+1/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 3; j < nj - 4; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(i, j + iTemp - 2, k)[0] * coef->GetEta(Idx(i, j + iTemp - 2, k))[0] +
                                            node->GetCoord(i, j + iTemp - 2, k)[1] * coef->GetEta(Idx(i, j + iTemp - 2, k))[1] +
                                            node->GetCoord(i, j + iTemp - 2, k)[2] * coef->GetEta(Idx(i, j + iTemp - 2, k))[2];
                    }
                    temp_j[Idx(i, j, k)] = MidNodeInter6th(inter_temp);
                }
            }
        }
        // j=1/2,3/2,5/2;nj-3/2,nj-5/2,nj-7/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = node->GetCoord(i, 1 + iTemp, k)[0] * coef->GetEta(Idx(i, 1 + iTemp, k))[0] +
                                        node->GetCoord(i, 1 + iTemp, k)[1] * coef->GetEta(Idx(i, 1 + iTemp, k))[1] +
                                        node->GetCoord(i, 1 + iTemp, k)[2] * coef->GetEta(Idx(i, 1 + iTemp, k))[2];
                }
                temp_j[Idx(i, 0, k)] = MidNodeInter4thRight1(inter_temp);
                temp_j[Idx(i, 1, k)] = MidNodeInter4thRight2(inter_temp);
                temp_j[Idx(i, 2, k)] = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = node->GetCoord(i, nj - 5 + iTemp, k)[0] * coef->GetEta(Idx(i, nj - 5 + iTemp, k))[0] +
                                        node->GetCoord(i, nj - 5 + iTemp, k)[1] * coef->GetEta(Idx(i, nj - 5 + iTemp, k))[1] +
                                        node->GetCoord(i, nj - 5 + iTemp, k)[2] * coef->GetEta(Idx(i, nj - 5 + iTemp, k))[2];
                }
                temp_j[Idx(i, nj - 2, k)] = MidNodeInter4thLeft1(inter_temp);
                temp_j[Idx(i, nj - 3, k)] = MidNodeInter4thLeft2(inter_temp);
                temp_j[Idx(i, nj - 4, k)] = MidNodeInter4th(inter_temp);
            }
        }
        // k+1/2
        for (int k = 3; k < nk - 4; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(i, j, k + iTemp - 2)[0] * coef->GetZeta(Idx(i, j, k + iTemp - 2))[0] +
                                            node->GetCoord(i, j, k + iTemp - 2)[1] * coef->GetZeta(Idx(i, j, k + iTemp - 2))[1] +
                                            node->GetCoord(i, j, k + iTemp - 2)[2] * coef->GetZeta(Idx(i, j, k + iTemp - 2))[2];
                    }
                    temp_k[Idx(i, j, k)] = MidNodeInter6th(inter_temp);
                }
            }
        }
        // k=1/2,3/2,5/2;nk-3/2,nk-5/2,nk-7/2
        for (int j = 1; j < nj - 1; ++j)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = node->GetCoord(i, j, 1 + iTemp)[0] * coef->GetZeta(Idx(i, j, 1 + iTemp))[0] +
                                        node->GetCoord(i, j, 1 + iTemp)[1] * coef->GetZeta(Idx(i, j, 1 + iTemp))[1] +
                                        node->GetCoord(i, j, 1 + iTemp)[2] * coef->GetZeta(Idx(i, j, 1 + iTemp))[2];
                }
                temp_k[Idx(i, j, 0)] = MidNodeInter4thRight1(inter_temp);
                temp_k[Idx(i, j, 1)] = MidNodeInter4thRight2(inter_temp);
                temp_k[Idx(i, j, 2)] = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = node->GetCoord(i, j, nk - 5 + iTemp)[0] * coef->GetZeta(Idx(i, j, nk - 5 + iTemp))[0] +
                                        node->GetCoord(i, j, nk - 5 + iTemp)[1] * coef->GetZeta(Idx(i, j, nk - 5 + iTemp))[1] +
                                        node->GetCoord(i, j, nk - 5 + iTemp)[2] * coef->GetZeta(Idx(i, j, nk - 5 + iTemp))[2];
                }
                temp_k[Idx(i, j, nk - 2)] = MidNodeInter4thLeft1(inter_temp);
                temp_k[Idx(i, j, nk - 3)] = MidNodeInter4thLeft2(inter_temp);
                temp_k[Idx(i, j, nk - 4)] = MidNodeInter4th(inter_temp);
            }
        }
        // 计算整点的Jacobian
        // 全部置零
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    coef->GetJacobian(Idx(i, j, k)) = 0.0;
                }
            }
        }
        double diff_temp[6];
        // i direction
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 3; i < ni - 3; ++i)
                {
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = temp_i[Idx(i + iTemp - 3, j, k)];
                    }
                    coef->GetJacobian(Idx(i, j, k)) += NodeDifferece6th(diff_temp);
                }
            }
        }
        // j direction
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 3; j < nj - 3; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = temp_j[Idx(i, j + iTemp - 3, k)];
                    }
                    coef->GetJacobian(Idx(i, j, k)) += NodeDifferece6th(diff_temp);
                }
            }
        }
        // k direction
        for (int k = 3; k < nk - 3; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = temp_k[Idx(i, j, k + iTemp - 3)];
                    }
                    coef->GetJacobian(Idx(i, j, k)) += NodeDifferece6th(diff_temp);
                }
            }
        }
        // 计算边界整点的Jacobian
        // i=1,2,ni-3,ni-2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = temp_i[Idx(0 + iTemp, j, k)];
                }
                coef->GetJacobian(Idx(1, j, k)) += NodeDifferece4thRight(diff_temp);
                coef->GetJacobian(Idx(2, j, k)) += NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = temp_i[Idx(ni - 2 - iTemp, j, k)];
                }
                coef->GetJacobian(Idx(ni - 2, j, k)) += NodeDifferece4thLeft(diff_temp);
                coef->GetJacobian(Idx(ni - 3, j, k)) += -NodeDifferece4th(diff_temp);
            }
        }
        // j=1,2,nj-3,nj-2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = temp_j[Idx(i, 0 + iTemp, k)];
                }
                coef->GetJacobian(Idx(i, 1, k)) += NodeDifferece4thRight(diff_temp);
                coef->GetJacobian(Idx(i, 2, k)) += NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = temp_j[Idx(i, nj - 2 - iTemp, k)];
                }
                coef->GetJacobian(Idx(i, nj - 2, k)) += NodeDifferece4thLeft(diff_temp);
                coef->GetJacobian(Idx(i, nj - 3, k)) += -NodeDifferece4th(diff_temp);
            }
        }
        // k=1,2,nk-3,nk-2
        for (int j = 1; j < nj - 1; ++j)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = temp_k[Idx(i, j, 0 + iTemp)];
                }
                coef->GetJacobian(Idx(i, j, 1)) += NodeDifferece4thRight(diff_temp);
                coef->GetJacobian(Idx(i, j, 2)) += NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = temp_k[Idx(i, j, nk - 2 - iTemp)];
                }
                coef->GetJacobian(Idx(i, j, nk - 2)) += NodeDifferece4thLeft(diff_temp);
                coef->GetJacobian(Idx(i, j, nk - 3)) += -NodeDifferece4th(diff_temp);
            }
        }
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    coef->GetJacobian(Idx(i, j, k)) = coef->GetJacobian(Idx(i, j, k)) / 3.0;
                    coef->GetJacobian(Idx(i, j, k)) = 1.0 / coef->GetJacobian(Idx(i, j, k));
                }
            }
        }
    }
    void NSSolverStruct::CalcMetricsS1()
    {
        auto para = GetPara();
        auto flux_diff_scheme = para->GetDifferenceScheme();
        if (flux_diff_scheme == DifferenceScheme::SecondOrder)
        {
            CalcMetricsS1_2nd();
        }
        else if (flux_diff_scheme == DifferenceScheme::SixthOrder)
        {
            CalcMetricsS1_6th();
        }
        else
        {
            Log::warn("FluxDifferenceScheme is not defined, use SecondOrder as default");
            CalcMetricsS1_2nd();
        }
    }
    void NSSolverStruct::CalcMetricsS1_2nd()
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
        int idx_temp[3];
        // i+1/2,j+1/2,k+1/2处的坐标
        std::vector<std::vector<double>>
            coord_i(ni * nj * nk, std::vector<double>(3)),
            coord_j(ni * nj * nk, std::vector<double>(3)), coord_k(ni * nj * nk, std::vector<double>(3));
        // 计算之前，先把度量系数赋值为0
        for (int i = 0; i < ni; ++i)
        {
            for (int j = 0; j < nj; ++j)
            {
                for (int k = 0; k < nk; ++k)
                {
                    for (int iDim = 0; iDim < 4; ++iDim)
                    {
                        coef->GetX(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetY(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetZ(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetXi(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetEta(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetZeta(Idx(i, j, k))[iDim] = 0.0;
                        for (int jDim = 0; jDim < 3; ++jDim)
                        {
                            coef_mid[jDim]->GetX(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetY(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetZ(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetXi(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetEta(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetZeta(Idx(i, j, k))[iDim] = 0.0;
                        }
                    }
                }
            }
        }
        // 第一步：计算半点坐标
        for (int k = 0; k < nk - 1; ++k)
        {
            for (int j = 0; j < nj - 1; ++j)
            {
                for (int i = 0; i < ni - 1; ++i)
                {
                    for (int iDim = 0; iDim < 3; ++iDim)
                    {
                        coord_i[Idx(i, j, k)][iDim] =
                            0.5 * (node->GetCoord(i, j, k)[iDim] + node->GetCoord(i + 1, j, k)[iDim]);
                        coord_j[Idx(i, j, k)][iDim] =
                            0.5 * (node->GetCoord(i, j, k)[iDim] + node->GetCoord(i, j + 1, k)[iDim]);
                        coord_k[Idx(i, j, k)][iDim] =
                            0.5 * (node->GetCoord(i, j, k)[iDim] + node->GetCoord(i, j, k + 1)[iDim]);
                    }
                }
            }
        }
        // 第二步：根据半点坐标计算整点逆变换度量系数
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    int idx = idx_proxy->GetIdx(i, j, k);
                    auto coef_x = coef->GetX(idx);
                    auto coef_y = coef->GetY(idx);
                    auto coef_z = coef->GetZ(idx);
                    coef_x[0] = coord_i[Idx(i, j, k)][0] - coord_i[Idx(i - 1, j, k)][0];
                    coef_x[1] = coord_j[Idx(i, j, k)][0] - coord_j[Idx(i, j - 1, k)][0];
                    coef_x[2] = coord_k[Idx(i, j, k)][0] - coord_k[Idx(i, j, k - 1)][0];
                    coef_y[0] = coord_i[Idx(i, j, k)][1] - coord_i[Idx(i - 1, j, k)][1];
                    coef_y[1] = coord_j[Idx(i, j, k)][1] - coord_j[Idx(i, j - 1, k)][1];
                    coef_y[2] = coord_k[Idx(i, j, k)][1] - coord_k[Idx(i, j, k - 1)][1];
                    coef_z[0] = coord_i[Idx(i, j, k)][2] - coord_i[Idx(i - 1, j, k)][2];
                    coef_z[1] = coord_j[Idx(i, j, k)][2] - coord_j[Idx(i, j - 1, k)][2];
                    coef_z[2] = coord_k[Idx(i, j, k)][2] - coord_k[Idx(i, j, k - 1)][2];
                    coef_x[3] = 0.0;
                    coef_y[3] = 0.0;
                    coef_z[3] = 0.0;
                    if (grid->GetDim() == 2)
                    {
                        coef_z[0] = 0.0;
                        coef_z[1] = 0.0;
                        coef_z[2] = 1.0;
                    }
                }
            }
        }
        // 第三步：计算半点处逆变换度量系数和坐标乘积
        //  求出xi方向的临时变量：逆变换度量系数乘以坐标
        struct TempXi
        {
            double y_zeta_z, z_zeta_x, x_zeta_y, y_eta_z, z_eta_x, x_eta_y;
        };
        std::vector<TempXi> temp_xi(ni * nj * nk);
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    int idx = Idx(i, j, k);
                    auto &temp = temp_xi[idx];
                    temp.y_zeta_z = 0.5 * (coef->GetY(Idx(i + 1, j, k))[2] * node->GetCoord(i + 1, j, k)[2] +
                                           coef->GetY(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[2]);
                    temp.z_zeta_x = 0.5 * (coef->GetZ(Idx(i + 1, j, k))[2] * node->GetCoord(i + 1, j, k)[0] +
                                           coef->GetZ(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[0]);
                    temp.x_zeta_y = 0.5 * (coef->GetX(Idx(i + 1, j, k))[2] * node->GetCoord(i + 1, j, k)[1] +
                                           coef->GetX(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[1]);
                    temp.y_eta_z = 0.5 * (coef->GetY(Idx(i + 1, j, k))[1] * node->GetCoord(i + 1, j, k)[2] +
                                          coef->GetY(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[2]);
                    temp.z_eta_x = 0.5 * (coef->GetZ(Idx(i + 1, j, k))[1] * node->GetCoord(i + 1, j, k)[0] +
                                          coef->GetZ(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[0]);
                    temp.x_eta_y = 0.5 * (coef->GetX(Idx(i + 1, j, k))[1] * node->GetCoord(i + 1, j, k)[1] +
                                          coef->GetX(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[1]);
                }
            }
        }
        // i=1/2；(ni-1)-1/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                int idx = Idx(0, j, k);
                temp_xi[idx].y_zeta_z = 0.5 * (3.0 * coef->GetY(Idx(1, j, k))[2] * node->GetCoord(1, j, k)[2] -
                                               1.0 * coef->GetY(Idx(2, j, k))[2] * node->GetCoord(2, j, k)[2]);
                temp_xi[idx].z_zeta_x = 0.5 * (3.0 * coef->GetZ(Idx(1, j, k))[2] * node->GetCoord(1, j, k)[0] -
                                               1.0 * coef->GetZ(Idx(2, j, k))[2] * node->GetCoord(2, j, k)[0]);
                temp_xi[idx].x_zeta_y = 0.5 * (3.0 * coef->GetX(Idx(1, j, k))[2] * node->GetCoord(1, j, k)[1] -
                                               1.0 * coef->GetX(Idx(2, j, k))[2] * node->GetCoord(2, j, k)[1]);
                temp_xi[idx].y_eta_z = 0.5 * (3.0 * coef->GetY(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[2] -
                                              1.0 * coef->GetY(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[2]);
                temp_xi[idx].z_eta_x = 0.5 * (3.0 * coef->GetZ(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[0] -
                                              1.0 * coef->GetZ(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[0]);
                temp_xi[idx].x_eta_y = 0.5 * (3.0 * coef->GetX(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[1] -
                                              1.0 * coef->GetX(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[1]);
                idx = Idx(ni - 2, j, k);
                temp_xi[idx].y_zeta_z = 0.5 * (3.0 * coef->GetY(Idx(ni - 2, j, k))[2] * node->GetCoord(ni - 2, j, k)[2] -
                                               1.0 * coef->GetY(Idx(ni - 3, j, k))[2] * node->GetCoord(ni - 3, j, k)[2]);
                temp_xi[idx].z_zeta_x = 0.5 * (3.0 * coef->GetZ(Idx(ni - 2, j, k))[2] * node->GetCoord(ni - 2, j, k)[0] -
                                               1.0 * coef->GetZ(Idx(ni - 3, j, k))[2] * node->GetCoord(ni - 3, j, k)[0]);
                temp_xi[idx].x_zeta_y = 0.5 * (3.0 * coef->GetX(Idx(ni - 2, j, k))[2] * node->GetCoord(ni - 2, j, k)[1] -
                                               1.0 * coef->GetX(Idx(ni - 3, j, k))[2] * node->GetCoord(ni - 3, j, k)[1]);
                temp_xi[idx].y_eta_z = 0.5 * (3.0 * coef->GetY(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[2] -
                                              1.0 * coef->GetY(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[2]);
                temp_xi[idx].z_eta_x = 0.5 * (3.0 * coef->GetZ(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[0] -
                                              1.0 * coef->GetZ(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[0]);
                temp_xi[idx].x_eta_y = 0.5 * (3.0 * coef->GetX(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[1] -
                                              1.0 * coef->GetX(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[1]);
            }
        }
        // 求出eta方向的临时变量：逆变换度量系数乘以坐标
        struct TempEta
        {
            double y_zeta_z, z_zeta_x, x_zeta_y, y_xi_z, z_xi_x, x_xi_y;
        };
        std::vector<TempEta> temp_eta(ni * nj * nk);
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    int idx = Idx(i, j, k);
                    auto &temp = temp_eta[idx];
                    temp.y_zeta_z = 0.5 * (coef->GetY(Idx(i, j + 1, k))[2] * node->GetCoord(i, j + 1, k)[2] +
                                           coef->GetY(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[2]);
                    temp.z_zeta_x = 0.5 * (coef->GetZ(Idx(i, j + 1, k))[2] * node->GetCoord(i, j + 1, k)[0] +
                                           coef->GetZ(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[0]);
                    temp.x_zeta_y = 0.5 * (coef->GetX(Idx(i, j + 1, k))[2] * node->GetCoord(i, j + 1, k)[1] +
                                           coef->GetX(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[1]);
                    temp.y_xi_z = 0.5 * (coef->GetY(Idx(i, j + 1, k))[0] * node->GetCoord(i, j + 1, k)[2] +
                                         coef->GetY(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[2]);
                    temp.z_xi_x = 0.5 * (coef->GetZ(Idx(i, j + 1, k))[0] * node->GetCoord(i, j + 1, k)[0] +
                                         coef->GetZ(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[0]);
                    temp.x_xi_y = 0.5 * (coef->GetX(Idx(i, j + 1, k))[0] * node->GetCoord(i, j + 1, k)[1] +
                                         coef->GetX(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[1]);
                }
            }
        }
        // j=1/2；(nj-1)-1/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                int idx = Idx(i, 0, k);
                temp_eta[idx].y_zeta_z = 0.5 * (3.0 * coef->GetY(Idx(i, 1, k))[2] * node->GetCoord(i, 1, k)[2] -
                                                1.0 * coef->GetY(Idx(i, 2, k))[2] * node->GetCoord(i, 2, k)[2]);
                temp_eta[idx].z_zeta_x = 0.5 * (3.0 * coef->GetZ(Idx(i, 1, k))[2] * node->GetCoord(i, 1, k)[0] -
                                                1.0 * coef->GetZ(Idx(i, 2, k))[2] * node->GetCoord(i, 2, k)[0]);
                temp_eta[idx].x_zeta_y = 0.5 * (3.0 * coef->GetX(Idx(i, 1, k))[2] * node->GetCoord(i, 1, k)[1] -
                                                1.0 * coef->GetX(Idx(i, 2, k))[2] * node->GetCoord(i, 2, k)[1]);
                temp_eta[idx].y_xi_z = 0.5 * (3.0 * coef->GetY(Idx(i, 1, k))[0] * node->GetCoord(i, 1, k)[2] -
                                              1.0 * coef->GetY(Idx(i, 2, k))[0] * node->GetCoord(i, 2, k)[2]);
                temp_eta[idx].z_xi_x = 0.5 * (3.0 * coef->GetZ(Idx(i, 1, k))[0] * node->GetCoord(i, 1, k)[0] -
                                              1.0 * coef->GetZ(Idx(i, 2, k))[0] * node->GetCoord(i, 2, k)[0]);
                temp_eta[idx].x_xi_y = 0.5 * (3.0 * coef->GetX(Idx(i, 1, k))[0] * node->GetCoord(i, 1, k)[1] -
                                              1.0 * coef->GetX(Idx(i, 2, k))[0] * node->GetCoord(i, 2, k)[1]);
                idx = Idx(i, nj - 2, k);
                temp_eta[idx].y_zeta_z = 0.5 * (3.0 * coef->GetY(Idx(i, nj - 2, k))[2] * node->GetCoord(i, nj - 2, k)[2] -
                                                1.0 * coef->GetY(Idx(i, nj - 3, k))[2] * node->GetCoord(i, nj - 3, k)[2]);
                temp_eta[idx].z_zeta_x = 0.5 * (3.0 * coef->GetZ(Idx(i, nj - 2, k))[2] * node->GetCoord(i, nj - 2, k)[0] -
                                                1.0 * coef->GetZ(Idx(i, nj - 3, k))[2] * node->GetCoord(i, nj - 3, k)[0]);
                temp_eta[idx].x_zeta_y = 0.5 * (3.0 * coef->GetX(Idx(i, nj - 2, k))[2] * node->GetCoord(i, nj - 2, k)[1] -
                                                1.0 * coef->GetX(Idx(i, nj - 3, k))[2] * node->GetCoord(i, nj - 3, k)[1]);
                temp_eta[idx].y_xi_z = 0.5 * (3.0 * coef->GetY(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[2] -
                                              1.0 * coef->GetY(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[2]);
                temp_eta[idx].z_xi_x = 0.5 * (3.0 * coef->GetZ(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[0] -
                                              1.0 * coef->GetZ(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[0]);
                temp_eta[idx].x_xi_y = 0.5 * (3.0 * coef->GetX(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[1] -
                                              1.0 * coef->GetX(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[1]);
            }
        }

        // 求出zeta方向的临时变量：逆变换度量系数乘以坐标
        struct TempZeta
        {
            double y_eta_z, z_eta_x, x_eta_y, y_xi_z, z_xi_x, x_xi_y;
        };
        std::vector<TempZeta> temp_zeta(ni * nj * nk);
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    int idx = Idx(i, j, k);
                    auto &temp = temp_zeta[idx];
                    temp.y_eta_z = 0.5 * (coef->GetY(Idx(i, j, k + 1))[1] * node->GetCoord(i, j, k + 1)[2] +
                                          coef->GetY(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[2]);
                    temp.z_eta_x = 0.5 * (coef->GetZ(Idx(i, j, k + 1))[1] * node->GetCoord(i, j, k + 1)[0] +
                                          coef->GetZ(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[0]);
                    temp.x_eta_y = 0.5 * (coef->GetX(Idx(i, j, k + 1))[1] * node->GetCoord(i, j, k + 1)[1] +
                                          coef->GetX(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[1]);
                    temp.y_xi_z = 0.5 * (coef->GetY(Idx(i, j, k + 1))[0] * node->GetCoord(i, j, k + 1)[2] +
                                         coef->GetY(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[2]);
                    temp.z_xi_x = 0.5 * (coef->GetZ(Idx(i, j, k + 1))[0] * node->GetCoord(i, j, k + 1)[0] +
                                         coef->GetZ(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[0]);
                    temp.x_xi_y = 0.5 * (coef->GetX(Idx(i, j, k + 1))[0] * node->GetCoord(i, j, k + 1)[1] +
                                         coef->GetX(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[1]);
                }
            }
        }
        // k=1/2；(nk-1)-1/2
        for (int j = 1; j < nj - 1; ++j)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                int idx = Idx(i, j, 0);
                temp_zeta[idx].y_eta_z = 0.5 * (3.0 * coef->GetY(Idx(i, j, 1))[1] * node->GetCoord(i, j, 1)[2] -
                                                1.0 * coef->GetY(Idx(i, j, 2))[1] * node->GetCoord(i, j, 2)[2]);
                temp_zeta[idx].z_eta_x = 0.5 * (3.0 * coef->GetZ(Idx(i, j, 1))[1] * node->GetCoord(i, j, 1)[0] -
                                                1.0 * coef->GetZ(Idx(i, j, 2))[1] * node->GetCoord(i, j, 2)[0]);
                temp_zeta[idx].x_eta_y = 0.5 * (3.0 * coef->GetX(Idx(i, j, 1))[1] * node->GetCoord(i, j, 1)[1] -
                                                1.0 * coef->GetX(Idx(i, j, 2))[1] * node->GetCoord(i, j, 2)[1]);
                temp_zeta[idx].y_xi_z = 0.5 * (3.0 * coef->GetY(Idx(i, j, 1))[0] * node->GetCoord(i, j, 1)[2] -
                                               1.0 * coef->GetY(Idx(i, j, 2))[0] * node->GetCoord(i, j, 2)[2]);
                temp_zeta[idx].z_xi_x = 0.5 * (3.0 * coef->GetZ(Idx(i, j, 1))[0] * node->GetCoord(i, j, 1)[0] -
                                               1.0 * coef->GetZ(Idx(i, j, 2))[0] * node->GetCoord(i, j, 2)[0]);
                temp_zeta[idx].x_xi_y = 0.5 * (3.0 * coef->GetX(Idx(i, j, 1))[0] * node->GetCoord(i, j, 1)[1] -
                                               1.0 * coef->GetX(Idx(i, j, 2))[0] * node->GetCoord(i, j, 2)[1]);
                idx = Idx(i, j, nk - 2);
                temp_zeta[idx].y_eta_z = 0.5 * (3.0 * coef->GetY(Idx(i, j, nk - 2))[1] * node->GetCoord(i, j, nk - 2)[2] -
                                                1.0 * coef->GetY(Idx(i, j, nk - 3))[1] * node->GetCoord(i, j, nk - 3)[2]);
                temp_zeta[idx].z_eta_x = 0.5 * (3.0 * coef->GetZ(Idx(i, j, nk - 2))[1] * node->GetCoord(i, j, nk - 2)[0] -
                                                1.0 * coef->GetZ(Idx(i, j, nk - 3))[1] * node->GetCoord(i, j, nk - 3)[0]);
                temp_zeta[idx].x_eta_y = 0.5 * (3.0 * coef->GetX(Idx(i, j, nk - 2))[1] * node->GetCoord(i, j, nk - 2)[1] -
                                                1.0 * coef->GetX(Idx(i, j, nk - 3))[1] * node->GetCoord(i, j, nk - 3)[1]);
                temp_zeta[idx].y_xi_z = 0.5 * (3.0 * coef->GetY(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[2] -
                                               1.0 * coef->GetY(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[2]);
                temp_zeta[idx].z_xi_x = 0.5 * (3.0 * coef->GetZ(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[0] -
                                               1.0 * coef->GetZ(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[0]);
                temp_zeta[idx].x_xi_y = 0.5 * (3.0 * coef->GetX(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[1] -
                                               1.0 * coef->GetX(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[1]);
            }
        }

        // 第四步：根据半点坐标和半点逆变换度量系数使用守恒形式计算整点度量系数（CMM1)
        double temp[6][2];
        // 获取下标的lamda函数
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    // i direction
                    auto coef_xi = coef->GetXi(Idx(i, j, k));
                    auto coef_eta = coef->GetEta(Idx(i, j, k));
                    auto coef_zeta = coef->GetZeta(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 2; iTemp++)
                    {
                        int idx = Idx(i - 1 + iTemp, j, k);
                        temp[0][iTemp] = temp_xi[idx].y_zeta_z; // y_zeta*z
                        temp[1][iTemp] = temp_xi[idx].z_zeta_x; // z_zeta*x
                        temp[2][iTemp] = temp_xi[idx].x_zeta_y; // x_zeta*y
                        temp[3][iTemp] = -temp_xi[idx].y_eta_z; // -y_eta*z
                        temp[4][iTemp] = -temp_xi[idx].z_eta_x; // -z_eta*x
                        temp[5][iTemp] = -temp_xi[idx].x_eta_y; // -x_eta*y
                    }
                    coef_eta[0] += temp[0][1] - temp[0][0];
                    coef_eta[1] += temp[1][1] - temp[1][0];
                    coef_eta[2] += temp[2][1] - temp[2][0];
                    coef_zeta[0] += temp[3][1] - temp[3][0];
                    coef_zeta[1] += temp[4][1] - temp[4][0];
                    coef_zeta[2] += temp[5][1] - temp[5][0];
                    for (int iTemp = 0; iTemp < 2; iTemp++)
                    {
                        int idx = Idx(i, j - 1 + iTemp, k);
                        temp[0][iTemp] = -temp_eta[idx].y_zeta_z; // -y_zeta*z
                        temp[1][iTemp] = -temp_eta[idx].z_zeta_x; // -z_zeta*x
                        temp[2][iTemp] = -temp_eta[idx].x_zeta_y; // -x_zeta*y
                        temp[3][iTemp] = temp_eta[idx].y_xi_z;    // y_xi*z
                        temp[4][iTemp] = temp_eta[idx].z_xi_x;    // z_xi*x
                        temp[5][iTemp] = temp_eta[idx].x_xi_y;    // x_xi*y
                    }
                    coef_xi[0] += temp[0][1] - temp[0][0];
                    coef_xi[1] += temp[1][1] - temp[1][0];
                    coef_xi[2] += temp[2][1] - temp[2][0];
                    coef_zeta[0] += temp[3][1] - temp[3][0];
                    coef_zeta[1] += temp[4][1] - temp[4][0];
                    coef_zeta[2] += temp[5][1] - temp[5][0];
                    for (int iTemp = 0; iTemp < 2; iTemp++)
                    {
                        int idx = Idx(i, j, k - 1 + iTemp);
                        temp[0][iTemp] = temp_zeta[idx].y_eta_z; // y_eta*z
                        temp[1][iTemp] = temp_zeta[idx].z_eta_x; // z_eta*x
                        temp[2][iTemp] = temp_zeta[idx].x_eta_y; // x_eta*y
                        temp[3][iTemp] = -temp_zeta[idx].y_xi_z; // -y_xi*z
                        temp[4][iTemp] = -temp_zeta[idx].z_xi_x; // -z_xi*x
                        temp[5][iTemp] = -temp_zeta[idx].x_xi_y; // -x_xi*y
                    }
                    coef_xi[0] += temp[0][1] - temp[0][0];
                    coef_xi[1] += temp[1][1] - temp[1][0];
                    coef_xi[2] += temp[2][1] - temp[2][0];
                    coef_eta[0] += temp[3][1] - temp[3][0];
                    coef_eta[1] += temp[4][1] - temp[4][0];
                    coef_eta[2] += temp[5][1] - temp[5][0];
                    if (grid->GetDim() == 2)
                    {
                        coef_zeta[0] = 0.0;
                        coef_zeta[1] = 0.0;
                        coef_zeta[2] = 1.0;
                    }
                    coef_xi[3] = coef_eta[3] = coef_zeta[3] = 0.0;
                }
            }
        }
    }
    void NSSolverStruct::CalcMetricsS1_6th()
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
        double inter_temp[6];
        // i+1/2,j+1/2,k+1/2处的坐标
        std::vector<std::vector<double>> coord_i(ni * nj * nk, std::vector<double>(3)),
            coord_j(ni * nj * nk, std::vector<double>(3)), coord_k(ni * nj * nk, std::vector<double>(3));
        // 计算之前，先把度量系数赋值为0
        for (int i = 0; i < ni; ++i)
        {
            for (int j = 0; j < nj; ++j)
            {
                for (int k = 0; k < nk; ++k)
                {
                    for (int iDim = 0; iDim < 4; ++iDim)
                    {
                        coef->GetX(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetY(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetZ(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetXi(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetEta(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetZeta(Idx(i, j, k))[iDim] = 0.0;
                        for (int jDim = 0; jDim < 3; ++jDim)
                        {
                            coef_mid[jDim]->GetX(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetY(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetZ(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetXi(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetEta(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetZeta(Idx(i, j, k))[iDim] = 0.0;
                        }
                    }
                }
            }
        }
        // 第一步：计算半点坐标
        // i+1/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 3; i < ni - 4; ++i)
                {
                    for (int iDim = 0; iDim < 3; ++iDim)
                    {
                        for (int iTemp = 0; iTemp < 6; iTemp++)
                        {
                            inter_temp[iTemp] = node->GetCoord(i + iTemp - 2, j, k)[iDim];
                        }
                        coord_i[Idx(i, j, k)][iDim] = MidNodeInter6th(inter_temp);
                    }
                }
            }
        }
        // i=1/2,3/2,5/2;ni-3/2,ni-5/2,ni-7/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int iDim = 0; iDim < 3; ++iDim)
                {
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(1 + iTemp, j, k)[iDim];
                    }
                    coord_i[Idx(0, j, k)][iDim] = MidNodeInter4thRight1(inter_temp);
                    coord_i[Idx(1, j, k)][iDim] = MidNodeInter4thRight2(inter_temp);
                    coord_i[Idx(2, j, k)][iDim] = MidNodeInter4th(inter_temp);
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(ni - 5 + iTemp, j, k)[iDim];
                    }
                    coord_i[Idx(ni - 2, j, k)][iDim] = MidNodeInter4thLeft1(inter_temp);
                    coord_i[Idx(ni - 3, j, k)][iDim] = MidNodeInter4thLeft2(inter_temp);
                    coord_i[Idx(ni - 4, j, k)][iDim] = MidNodeInter4th(inter_temp);
                }
            }
        }
        // j+1/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 3; j < nj - 4; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    for (int iDim = 0; iDim < 3; ++iDim)
                    {
                        for (int iTemp = 0; iTemp < 6; iTemp++)
                        {
                            inter_temp[iTemp] = node->GetCoord(i, j + iTemp - 2, k)[iDim];
                        }
                        coord_j[Idx(i, j, k)][iDim] = MidNodeInter6th(inter_temp);
                    }
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
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(i, 1 + iTemp, k)[iDim];
                    }
                    coord_j[Idx(i, 0, k)][iDim] = MidNodeInter4thRight1(inter_temp);
                    coord_j[Idx(i, 1, k)][iDim] = MidNodeInter4thRight2(inter_temp);
                    coord_j[Idx(i, 2, k)][iDim] = MidNodeInter4th(inter_temp);
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(i, nj - 5 + iTemp, k)[iDim];
                    }
                    coord_j[Idx(i, nj - 2, k)][iDim] = MidNodeInter4thLeft1(inter_temp);
                    coord_j[Idx(i, nj - 3, k)][iDim] = MidNodeInter4thLeft2(inter_temp);
                    coord_j[Idx(i, nj - 4, k)][iDim] = MidNodeInter4th(inter_temp);
                }
            }
        }
        // k+1/2
        for (int k = 3; k < nk - 4; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    for (int iDim = 0; iDim < 3; ++iDim)
                    {
                        for (int iTemp = 0; iTemp < 6; iTemp++)
                        {
                            inter_temp[iTemp] = node->GetCoord(i, j, k + iTemp - 2)[iDim];
                        }
                        coord_k[Idx(i, j, k)][iDim] = MidNodeInter6th(inter_temp);
                    }
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
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(i, j, 1 + iTemp)[iDim];
                    }
                    coord_k[Idx(i, j, 0)][iDim] = MidNodeInter4thRight1(inter_temp);
                    coord_k[Idx(i, j, 1)][iDim] = MidNodeInter4thRight2(inter_temp);
                    coord_k[Idx(i, j, 2)][iDim] = MidNodeInter4th(inter_temp);
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(i, j, nk - 5 + iTemp)[iDim];
                    }
                    coord_k[Idx(i, j, nk - 2)][iDim] = MidNodeInter4thLeft1(inter_temp);
                    coord_k[Idx(i, j, nk - 3)][iDim] = MidNodeInter4thLeft2(inter_temp);
                    coord_k[Idx(i, j, nk - 4)][iDim] = MidNodeInter4th(inter_temp);
                }
            }
        }
        // 第二步：根据半点坐标计算整点逆变换度量系数
        double diff_temp[6];
        // i direction
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 3; i < ni - 3; ++i)
                {
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_i[Idx(i + iTemp - 3, j, k)][0];
                    }
                    auto coef_x = coef->GetX(Idx(i, j, k));
                    coef_x[0] = NodeDifferece6th(diff_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_i[Idx(i + iTemp - 3, j, k)][1];
                    }
                    auto coef_y = coef->GetY(Idx(i, j, k));
                    coef_y[0] = NodeDifferece6th(diff_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_i[Idx(i + iTemp - 3, j, k)][2];
                    }
                    auto coef_z = coef->GetZ(Idx(i, j, k));
                    coef_z[0] = NodeDifferece6th(diff_temp);
                }
            }
        }
        // j direction
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 3; j < nj - 3; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    auto coef_x = coef->GetX(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_j[Idx(i, j + iTemp - 3, k)][0];
                    }
                    coef_x[1] = NodeDifferece6th(diff_temp);
                    auto coef_y = coef->GetY(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_j[Idx(i, j + iTemp - 3, k)][1];
                    }
                    coef_y[1] = NodeDifferece6th(diff_temp);
                    auto coef_z = coef->GetZ(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_j[Idx(i, j + iTemp - 3, k)][2];
                    }
                    coef_z[1] = NodeDifferece6th(diff_temp);
                }
            }
        }
        // k direction
        for (int k = 3; k < nk - 3; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    auto coef_x = coef->GetX(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_k[Idx(i, j, k + iTemp - 3)][0];
                    }
                    coef_x[2] = NodeDifferece6th(diff_temp);
                    auto coef_y = coef->GetY(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_k[Idx(i, j, k + iTemp - 3)][1];
                    }
                    coef_y[2] = NodeDifferece6th(diff_temp);
                    auto coef_z = coef->GetZ(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_k[Idx(i, j, k + iTemp - 3)][2];
                    }
                    coef_z[2] = NodeDifferece6th(diff_temp);
                    if (grid->GetDim() == 2)
                    {
                        coef_x[2] = 0.0;
                        coef_y[2] = 0.0;
                        coef_z[2] = 1.0;
                    }
                }
            }
        }
        // 计算边界整点的度量系数
        // i=1,2;ni-3,ni-2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_i[Idx(0 + iTemp, j, k)][0];
                }
                coef->GetX(Idx(1, j, k))[0] = NodeDifferece4thRight(diff_temp);
                coef->GetX(Idx(2, j, k))[0] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_i[Idx(0 + iTemp, j, k)][1];
                }
                coef->GetY(Idx(1, j, k))[0] = NodeDifferece4thRight(diff_temp);
                coef->GetY(Idx(2, j, k))[0] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_i[Idx(0 + iTemp, j, k)][2];
                }
                coef->GetZ(Idx(1, j, k))[0] = NodeDifferece4thRight(diff_temp);
                coef->GetZ(Idx(2, j, k))[0] = NodeDifferece4th(diff_temp);

                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_i[Idx(ni - 2 - iTemp, j, k)][0];
                }
                coef->GetX(Idx(ni - 2, j, k))[0] = NodeDifferece4thLeft(diff_temp);
                coef->GetX(Idx(ni - 3, j, k))[0] = -NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_i[Idx(ni - 2 - iTemp, j, k)][1];
                }
                coef->GetY(Idx(ni - 2, j, k))[0] = NodeDifferece4thLeft(diff_temp);
                coef->GetY(Idx(ni - 3, j, k))[0] = -NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_i[Idx(ni - 2 - iTemp, j, k)][2];
                }
                coef->GetZ(Idx(ni - 2, j, k))[0] = NodeDifferece4thLeft(diff_temp);
                coef->GetZ(Idx(ni - 3, j, k))[0] = -NodeDifferece4th(diff_temp);
            }
        }
        // j=1,2;nj-3,nj-2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_j[Idx(i, 0 + iTemp, k)][0];
                }
                coef->GetX(Idx(i, 1, k))[1] = NodeDifferece4thRight(diff_temp);
                coef->GetX(Idx(i, 2, k))[1] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_j[Idx(i, 0 + iTemp, k)][1];
                }
                coef->GetY(Idx(i, 1, k))[1] = NodeDifferece4thRight(diff_temp);
                coef->GetY(Idx(i, 2, k))[1] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_j[Idx(i, 0 + iTemp, k)][2];
                }
                coef->GetZ(Idx(i, 1, k))[1] = NodeDifferece4thRight(diff_temp);
                coef->GetZ(Idx(i, 2, k))[1] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_j[Idx(i, nj - 2 - iTemp, k)][0];
                }
                coef->GetX(Idx(i, nj - 2, k))[1] = NodeDifferece4thLeft(diff_temp);
                coef->GetX(Idx(i, nj - 3, k))[1] = -NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_j[Idx(i, nj - 2 - iTemp, k)][1];
                }
                coef->GetY(Idx(i, nj - 2, k))[1] = NodeDifferece4thLeft(diff_temp);
                coef->GetY(Idx(i, nj - 3, k))[1] = -NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_j[Idx(i, nj - 2 - iTemp, k)][2];
                }
                coef->GetZ(Idx(i, nj - 2, k))[1] = NodeDifferece4thLeft(diff_temp);
                coef->GetZ(Idx(i, nj - 3, k))[1] = -NodeDifferece4th(diff_temp);
            }
        }
        // k=1,2;nk-3,nk-2
        for (int j = 1; j < nj - 1; ++j)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_k[Idx(i, j, 0 + iTemp)][0];
                }
                coef->GetX(Idx(i, j, 1))[2] = NodeDifferece4thRight(diff_temp);
                coef->GetX(Idx(i, j, 2))[2] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_k[Idx(i, j, 0 + iTemp)][1];
                }
                coef->GetY(Idx(i, j, 1))[2] = NodeDifferece4thRight(diff_temp);
                coef->GetY(Idx(i, j, 2))[2] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_k[Idx(i, j, 0 + iTemp)][2];
                }
                coef->GetZ(Idx(i, j, 1))[2] = NodeDifferece4thRight(diff_temp);
                coef->GetZ(Idx(i, j, 2))[2] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_k[Idx(i, j, nk - 2 - iTemp)][0];
                }
                coef->GetX(Idx(i, j, nk - 2))[2] = NodeDifferece4thLeft(diff_temp);
                coef->GetX(Idx(i, j, nk - 3))[2] = -NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_k[Idx(i, j, nk - 2 - iTemp)][1];
                }
                coef->GetY(Idx(i, j, nk - 2))[2] = NodeDifferece4thLeft(diff_temp);
                coef->GetY(Idx(i, j, nk - 3))[2] = -NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_k[Idx(i, j, nk - 2 - iTemp)][2];
                }
                coef->GetZ(Idx(i, j, nk - 2))[2] = NodeDifferece4thLeft(diff_temp);
                coef->GetZ(Idx(i, j, nk - 3))[2] = -NodeDifferece4th(diff_temp);
                if (grid->GetDim() == 2)
                {
                    coef->GetZ(Idx(i, j, 1))[0] = coef->GetZ(Idx(i, j, 1))[1] = 0.0;
                    coef->GetZ(Idx(i, j, 2))[0] = coef->GetZ(Idx(i, j, 2))[1] = 0.0;
                    coef->GetZ(Idx(i, j, nk - 2))[0] = coef->GetZ(Idx(i, j, nk - 2))[1] = 0.0;
                    coef->GetZ(Idx(i, j, nk - 3))[0] = coef->GetZ(Idx(i, j, nk - 3))[1] = 0.0;
                    coef->GetZ(Idx(i, j, 1))[2] = 1.0;
                    coef->GetZ(Idx(i, j, 2))[2] = 1.0;
                    coef->GetZ(Idx(i, j, nk - 2))[2] = 1.0;
                    coef->GetZ(Idx(i, j, nk - 3))[2] = 1.0;
                }
            }
        }

        // 第三步：计算半点处逆变换度量系数和坐标乘积
        //  求出xi方向的临时变量：逆变换度量系数乘以坐标
        struct TempXi
        {
            double y_zeta_z, z_zeta_x, x_zeta_y, y_eta_z, z_eta_x, x_eta_y;
        };
        std::vector<TempXi> temp_xi(ni * nj * nk);
        // i+1/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 3; i < ni - 4; ++i)
                {
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetY(Idx(i + iTemp - 2, j, k))[2] * node->GetCoord(i + iTemp - 2, j, k)[2];
                    }
                    temp_xi[Idx(i, j, k)].y_zeta_z = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZ(Idx(i + iTemp - 2, j, k))[2] * node->GetCoord(i + iTemp - 2, j, k)[0];
                    }
                    temp_xi[Idx(i, j, k)].z_zeta_x = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetX(Idx(i + iTemp - 2, j, k))[2] * node->GetCoord(i + iTemp - 2, j, k)[1];
                    }
                    temp_xi[Idx(i, j, k)].x_zeta_y = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetY(Idx(i + iTemp - 2, j, k))[1] * node->GetCoord(i + iTemp - 2, j, k)[2];
                    }
                    temp_xi[Idx(i, j, k)].y_eta_z = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZ(Idx(i + iTemp - 2, j, k))[1] * node->GetCoord(i + iTemp - 2, j, k)[0];
                    }
                    temp_xi[Idx(i, j, k)].z_eta_x = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetX(Idx(i + iTemp - 2, j, k))[1] * node->GetCoord(i + iTemp - 2, j, k)[1];
                    }
                    temp_xi[Idx(i, j, k)].x_eta_y = MidNodeInter6th(inter_temp);
                }
            }
        }
        // i=1/2,3/2,5/2;ni-3/2,ni-5/2,ni-7/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                // i=1/2，3/2，5/2
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(1 + iTemp, j, k))[2] * node->GetCoord(1 + iTemp, j, k)[2];
                }
                temp_xi[Idx(0, j, k)].y_zeta_z = MidNodeInter4thRight1(inter_temp);
                temp_xi[Idx(1, j, k)].y_zeta_z = MidNodeInter4thRight2(inter_temp);
                temp_xi[Idx(2, j, k)].y_zeta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(1 + iTemp, j, k))[2] * node->GetCoord(1 + iTemp, j, k)[0];
                }
                temp_xi[Idx(0, j, k)].z_zeta_x = MidNodeInter4thRight1(inter_temp);
                temp_xi[Idx(1, j, k)].z_zeta_x = MidNodeInter4thRight2(inter_temp);
                temp_xi[Idx(2, j, k)].z_zeta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(1 + iTemp, j, k))[2] * node->GetCoord(1 + iTemp, j, k)[1];
                }
                temp_xi[Idx(0, j, k)].x_zeta_y = MidNodeInter4thRight1(inter_temp);
                temp_xi[Idx(1, j, k)].x_zeta_y = MidNodeInter4thRight2(inter_temp);
                temp_xi[Idx(2, j, k)].x_zeta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(1 + iTemp, j, k))[1] * node->GetCoord(1 + iTemp, j, k)[2];
                }
                temp_xi[Idx(0, j, k)].y_eta_z = MidNodeInter4thRight1(inter_temp);
                temp_xi[Idx(1, j, k)].y_eta_z = MidNodeInter4thRight2(inter_temp);
                temp_xi[Idx(2, j, k)].y_eta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(1 + iTemp, j, k))[1] * node->GetCoord(1 + iTemp, j, k)[0];
                }
                temp_xi[Idx(0, j, k)].z_eta_x = MidNodeInter4thRight1(inter_temp);
                temp_xi[Idx(1, j, k)].z_eta_x = MidNodeInter4thRight2(inter_temp);
                temp_xi[Idx(2, j, k)].z_eta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(1 + iTemp, j, k))[1] * node->GetCoord(1 + iTemp, j, k)[1];
                }
                temp_xi[Idx(0, j, k)].x_eta_y = MidNodeInter4thRight1(inter_temp);
                temp_xi[Idx(1, j, k)].x_eta_y = MidNodeInter4thRight2(inter_temp);
                temp_xi[Idx(2, j, k)].x_eta_y = MidNodeInter4th(inter_temp);
                // i=ni-3/2,ni-5/2,ni-7/2
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(ni - 5 + iTemp, j, k))[2] * node->GetCoord(ni - 5 + iTemp, j, k)[2];
                }
                temp_xi[Idx(ni - 2, j, k)].y_zeta_z = MidNodeInter4thLeft1(inter_temp);
                temp_xi[Idx(ni - 3, j, k)].y_zeta_z = MidNodeInter4thLeft2(inter_temp);
                temp_xi[Idx(ni - 4, j, k)].y_zeta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(ni - 5 + iTemp, j, k))[2] * node->GetCoord(ni - 5 + iTemp, j, k)[0];
                }
                temp_xi[Idx(ni - 2, j, k)].z_zeta_x = MidNodeInter4thLeft1(inter_temp);
                temp_xi[Idx(ni - 3, j, k)].z_zeta_x = MidNodeInter4thLeft2(inter_temp);
                temp_xi[Idx(ni - 4, j, k)].z_zeta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(ni - 5 + iTemp, j, k))[2] * node->GetCoord(ni - 5 + iTemp, j, k)[1];
                }
                temp_xi[Idx(ni - 2, j, k)].x_zeta_y = MidNodeInter4thLeft1(inter_temp);
                temp_xi[Idx(ni - 3, j, k)].x_zeta_y = MidNodeInter4thLeft2(inter_temp);
                temp_xi[Idx(ni - 4, j, k)].x_zeta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(ni - 5 + iTemp, j, k))[1] * node->GetCoord(ni - 5 + iTemp, j, k)[2];
                }
                temp_xi[Idx(ni - 2, j, k)].y_eta_z = MidNodeInter4thLeft1(inter_temp);
                temp_xi[Idx(ni - 3, j, k)].y_eta_z = MidNodeInter4thLeft2(inter_temp);
                temp_xi[Idx(ni - 4, j, k)].y_eta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(ni - 5 + iTemp, j, k))[1] * node->GetCoord(ni - 5 + iTemp, j, k)[0];
                }
                temp_xi[Idx(ni - 2, j, k)].z_eta_x = MidNodeInter4thLeft1(inter_temp);
                temp_xi[Idx(ni - 3, j, k)].z_eta_x = MidNodeInter4thLeft2(inter_temp);
                temp_xi[Idx(ni - 4, j, k)].z_eta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(ni - 5 + iTemp, j, k))[1] * node->GetCoord(ni - 5 + iTemp, j, k)[1];
                }
                temp_xi[Idx(ni - 2, j, k)].x_eta_y = MidNodeInter4thLeft1(inter_temp);
                temp_xi[Idx(ni - 3, j, k)].x_eta_y = MidNodeInter4thLeft2(inter_temp);
                temp_xi[Idx(ni - 4, j, k)].x_eta_y = MidNodeInter4th(inter_temp);
            }
        }
        // 求出eta方向的临时变量：逆变换度量系数乘以坐标
        struct TempEta
        {
            double y_zeta_z, z_zeta_x, x_zeta_y, y_xi_z, z_xi_x, x_xi_y;
        };
        std::vector<TempEta> temp_eta(ni * nj * nk);
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 3; j < nj - 4; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetY(Idx(i, j + iTemp - 2, k))[2] * node->GetCoord(i, j + iTemp - 2, k)[2];
                    }
                    temp_eta[Idx(i, j, k)].y_zeta_z = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZ(Idx(i, j + iTemp - 2, k))[2] * node->GetCoord(i, j + iTemp - 2, k)[0];
                    }
                    temp_eta[Idx(i, j, k)].z_zeta_x = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetX(Idx(i, j + iTemp - 2, k))[2] * node->GetCoord(i, j + iTemp - 2, k)[1];
                    }
                    temp_eta[Idx(i, j, k)].x_zeta_y = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetY(Idx(i, j + iTemp - 2, k))[0] * node->GetCoord(i, j + iTemp - 2, k)[2];
                    }
                    temp_eta[Idx(i, j, k)].y_xi_z = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZ(Idx(i, j + iTemp - 2, k))[0] * node->GetCoord(i, j + iTemp - 2, k)[0];
                    }
                    temp_eta[Idx(i, j, k)].z_xi_x = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetX(Idx(i, j + iTemp - 2, k))[0] * node->GetCoord(i, j + iTemp - 2, k)[1];
                    }
                    temp_eta[Idx(i, j, k)].x_xi_y = MidNodeInter6th(inter_temp);
                }
            }
        }
        // j=1/2,3/2,5/2;nj-3/2,nj-5/2,nj-7/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                // j=1/2，3/2，5/2
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, 1 + iTemp, k))[2] * node->GetCoord(i, 1 + iTemp, k)[2];
                }
                temp_eta[Idx(i, 0, k)].y_zeta_z = MidNodeInter4thRight1(inter_temp);
                temp_eta[Idx(i, 1, k)].y_zeta_z = MidNodeInter4thRight2(inter_temp);
                temp_eta[Idx(i, 2, k)].y_zeta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, 1 + iTemp, k))[2] * node->GetCoord(i, 1 + iTemp, k)[0];
                }
                temp_eta[Idx(i, 0, k)].z_zeta_x = MidNodeInter4thRight1(inter_temp);
                temp_eta[Idx(i, 1, k)].z_zeta_x = MidNodeInter4thRight2(inter_temp);
                temp_eta[Idx(i, 2, k)].z_zeta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, 1 + iTemp, k))[2] * node->GetCoord(i, 1 + iTemp, k)[1];
                }
                temp_eta[Idx(i, 0, k)].x_zeta_y = MidNodeInter4thRight1(inter_temp);
                temp_eta[Idx(i, 1, k)].x_zeta_y = MidNodeInter4thRight2(inter_temp);
                temp_eta[Idx(i, 2, k)].x_zeta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, 1 + iTemp, k))[0] * node->GetCoord(i, 1 + iTemp, k)[2];
                }
                temp_eta[Idx(i, 0, k)].y_xi_z = MidNodeInter4thRight1(inter_temp);
                temp_eta[Idx(i, 1, k)].y_xi_z = MidNodeInter4thRight2(inter_temp);
                temp_eta[Idx(i, 2, k)].y_xi_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, 1 + iTemp, k))[0] * node->GetCoord(i, 1 + iTemp, k)[0];
                }
                temp_eta[Idx(i, 0, k)].z_xi_x = MidNodeInter4thRight1(inter_temp);
                temp_eta[Idx(i, 1, k)].z_xi_x = MidNodeInter4thRight2(inter_temp);
                temp_eta[Idx(i, 2, k)].z_xi_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, 1 + iTemp, k))[0] * node->GetCoord(i, 1 + iTemp, k)[1];
                }
                temp_eta[Idx(i, 0, k)].x_xi_y = MidNodeInter4thRight1(inter_temp);
                temp_eta[Idx(i, 1, k)].x_xi_y = MidNodeInter4thRight2(inter_temp);
                temp_eta[Idx(i, 2, k)].x_xi_y = MidNodeInter4th(inter_temp);
                // j=nj-3/2,nj-5/2,nj-7/2
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, nj - 5 + iTemp, k))[2] * node->GetCoord(i, nj - 5 + iTemp, k)[2];
                }
                temp_eta[Idx(i, nj - 2, k)].y_zeta_z = MidNodeInter4thLeft1(inter_temp);
                temp_eta[Idx(i, nj - 3, k)].y_zeta_z = MidNodeInter4thLeft2(inter_temp);
                temp_eta[Idx(i, nj - 4, k)].y_zeta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, nj - 5 + iTemp, k))[2] * node->GetCoord(i, nj - 5 + iTemp, k)[0];
                }
                temp_eta[Idx(i, nj - 2, k)].z_zeta_x = MidNodeInter4thLeft1(inter_temp);
                temp_eta[Idx(i, nj - 3, k)].z_zeta_x = MidNodeInter4thLeft2(inter_temp);
                temp_eta[Idx(i, nj - 4, k)].z_zeta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, nj - 5 + iTemp, k))[2] * node->GetCoord(i, nj - 5 + iTemp, k)[1];
                }
                temp_eta[Idx(i, nj - 2, k)].x_zeta_y = MidNodeInter4thLeft1(inter_temp);
                temp_eta[Idx(i, nj - 3, k)].x_zeta_y = MidNodeInter4thLeft2(inter_temp);
                temp_eta[Idx(i, nj - 4, k)].x_zeta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, nj - 5 + iTemp, k))[0] * node->GetCoord(i, nj - 5 + iTemp, k)[2];
                }
                temp_eta[Idx(i, nj - 2, k)].y_xi_z = MidNodeInter4thLeft1(inter_temp);
                temp_eta[Idx(i, nj - 3, k)].y_xi_z = MidNodeInter4thLeft2(inter_temp);
                temp_eta[Idx(i, nj - 4, k)].y_xi_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, nj - 5 + iTemp, k))[0] * node->GetCoord(i, nj - 5 + iTemp, k)[0];
                }
                temp_eta[Idx(i, nj - 2, k)].z_xi_x = MidNodeInter4thLeft1(inter_temp);
                temp_eta[Idx(i, nj - 3, k)].z_xi_x = MidNodeInter4thLeft2(inter_temp);
                temp_eta[Idx(i, nj - 4, k)].z_xi_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, nj - 5 + iTemp, k))[0] * node->GetCoord(i, nj - 5 + iTemp, k)[1];
                }
                temp_eta[Idx(i, nj - 2, k)].x_xi_y = MidNodeInter4thLeft1(inter_temp);
                temp_eta[Idx(i, nj - 3, k)].x_xi_y = MidNodeInter4thLeft2(inter_temp);
                temp_eta[Idx(i, nj - 4, k)].x_xi_y = MidNodeInter4th(inter_temp);
            }
        }

        // 求出zeta方向的临时变量：逆变换度量系数乘以坐标
        struct TempZeta
        {
            double y_eta_z, z_eta_x, x_eta_y, y_xi_z, z_xi_x, x_xi_y;
        };
        std::vector<TempZeta> temp_zeta(ni * nj * nk);
        for (int k = 3; k < nk - 4; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetY(Idx(i, j, k + iTemp - 2))[1] * node->GetCoord(i, j, k + iTemp - 2)[2];
                    }
                    temp_zeta[Idx(i, j, k)].y_eta_z = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZ(Idx(i, j, k + iTemp - 2))[1] * node->GetCoord(i, j, k + iTemp - 2)[0];
                    }
                    temp_zeta[Idx(i, j, k)].z_eta_x = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetX(Idx(i, j, k + iTemp - 2))[1] * node->GetCoord(i, j, k + iTemp - 2)[1];
                    }
                    temp_zeta[Idx(i, j, k)].x_eta_y = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetY(Idx(i, j, k + iTemp - 2))[0] * node->GetCoord(i, j, k + iTemp - 2)[2];
                    }
                    temp_zeta[Idx(i, j, k)].y_xi_z = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZ(Idx(i, j, k + iTemp - 2))[0] * node->GetCoord(i, j, k + iTemp - 2)[0];
                    }
                    temp_zeta[Idx(i, j, k)].z_xi_x = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetX(Idx(i, j, k + iTemp - 2))[0] * node->GetCoord(i, j, k + iTemp - 2)[1];
                    }
                    temp_zeta[Idx(i, j, k)].x_xi_y = MidNodeInter6th(inter_temp);
                }
            }
        }
        // k=1/2,3/2,5/2;nk-3/2,nk-5/2,nk-7/2
        for (int j = 1; j < nj - 1; ++j)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                // k=1/2，3/2，5/2
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, j, 1 + iTemp))[1] * node->GetCoord(i, j, 1 + iTemp)[2];
                }
                temp_zeta[Idx(i, j, 0)].y_eta_z = MidNodeInter4thRight1(inter_temp);
                temp_zeta[Idx(i, j, 1)].y_eta_z = MidNodeInter4thRight2(inter_temp);
                temp_zeta[Idx(i, j, 2)].y_eta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, j, 1 + iTemp))[1] * node->GetCoord(i, j, 1 + iTemp)[0];
                }
                temp_zeta[Idx(i, j, 0)].z_eta_x = MidNodeInter4thRight1(inter_temp);
                temp_zeta[Idx(i, j, 1)].z_eta_x = MidNodeInter4thRight2(inter_temp);
                temp_zeta[Idx(i, j, 2)].z_eta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, j, 1 + iTemp))[1] * node->GetCoord(i, j, 1 + iTemp)[1];
                }
                temp_zeta[Idx(i, j, 0)].x_eta_y = MidNodeInter4thRight1(inter_temp);
                temp_zeta[Idx(i, j, 1)].x_eta_y = MidNodeInter4thRight2(inter_temp);
                temp_zeta[Idx(i, j, 2)].x_eta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, j, 1 + iTemp))[0] * node->GetCoord(i, j, 1 + iTemp)[2];
                }
                temp_zeta[Idx(i, j, 0)].y_xi_z = MidNodeInter4thRight1(inter_temp);
                temp_zeta[Idx(i, j, 1)].y_xi_z = MidNodeInter4thRight2(inter_temp);
                temp_zeta[Idx(i, j, 2)].y_xi_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, j, 1 + iTemp))[0] * node->GetCoord(i, j, 1 + iTemp)[0];
                }
                temp_zeta[Idx(i, j, 0)].z_xi_x = MidNodeInter4thRight1(inter_temp);
                temp_zeta[Idx(i, j, 1)].z_xi_x = MidNodeInter4thRight2(inter_temp);
                temp_zeta[Idx(i, j, 2)].z_xi_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, j, 1 + iTemp))[0] * node->GetCoord(i, j, 1 + iTemp)[1];
                }
                temp_zeta[Idx(i, j, 0)].x_xi_y = MidNodeInter4thRight1(inter_temp);
                temp_zeta[Idx(i, j, 1)].x_xi_y = MidNodeInter4thRight2(inter_temp);
                temp_zeta[Idx(i, j, 2)].x_xi_y = MidNodeInter4th(inter_temp);
                // k=nk-3/2,nk-5/2,nk-7/2
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, j, nk - 5 + iTemp))[1] * node->GetCoord(i, j, nk - 5 + iTemp)[2];
                }
                temp_zeta[Idx(i, j, nk - 2)].y_eta_z = MidNodeInter4thLeft1(inter_temp);
                temp_zeta[Idx(i, j, nk - 3)].y_eta_z = MidNodeInter4thLeft2(inter_temp);
                temp_zeta[Idx(i, j, nk - 4)].y_eta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, j, nk - 5 + iTemp))[1] * node->GetCoord(i, j, nk - 5 + iTemp)[0];
                }
                temp_zeta[Idx(i, j, nk - 2)].z_eta_x = MidNodeInter4thLeft1(inter_temp);
                temp_zeta[Idx(i, j, nk - 3)].z_eta_x = MidNodeInter4thLeft2(inter_temp);
                temp_zeta[Idx(i, j, nk - 4)].z_eta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, j, nk - 5 + iTemp))[1] * node->GetCoord(i, j, nk - 5 + iTemp)[1];
                }
                temp_zeta[Idx(i, j, nk - 2)].x_eta_y = MidNodeInter4thLeft1(inter_temp);
                temp_zeta[Idx(i, j, nk - 3)].x_eta_y = MidNodeInter4thLeft2(inter_temp);
                temp_zeta[Idx(i, j, nk - 4)].x_eta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, j, nk - 5 + iTemp))[0] * node->GetCoord(i, j, nk - 5 + iTemp)[2];
                }
                temp_zeta[Idx(i, j, nk - 2)].y_xi_z = MidNodeInter4thLeft1(inter_temp);
                temp_zeta[Idx(i, j, nk - 3)].y_xi_z = MidNodeInter4thLeft2(inter_temp);
                temp_zeta[Idx(i, j, nk - 4)].y_xi_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, j, nk - 5 + iTemp))[0] * node->GetCoord(i, j, nk - 5 + iTemp)[0];
                }
                temp_zeta[Idx(i, j, nk - 2)].z_xi_x = MidNodeInter4thLeft1(inter_temp);
                temp_zeta[Idx(i, j, nk - 3)].z_xi_x = MidNodeInter4thLeft2(inter_temp);
                temp_zeta[Idx(i, j, nk - 4)].z_xi_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, j, nk - 5 + iTemp))[0] * node->GetCoord(i, j, nk - 5 + iTemp)[1];
                }
                temp_zeta[Idx(i, j, nk - 2)].x_xi_y = MidNodeInter4thLeft1(inter_temp);
                temp_zeta[Idx(i, j, nk - 3)].x_xi_y = MidNodeInter4thLeft2(inter_temp);
                temp_zeta[Idx(i, j, nk - 4)].x_xi_y = MidNodeInter4th(inter_temp);
            }
        }
        // 第五步：根据半点坐标和半点逆变换度量系数使用守恒形式计算整点度量系数（CMM1)
        // i direction
        double temp[6][6];
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 3; i < ni - 3; ++i)
                {
                    auto coef_xi = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, k));
                    auto coef_eta = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, k));
                    auto coef_zeta = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, k));
                    // i direction
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        int idx = Idx(i - 3 + iTemp, j, k);
                        temp[0][iTemp] = temp_xi[idx].y_zeta_z; // y_zeta*z
                        temp[1][iTemp] = temp_xi[idx].z_zeta_x; // z_zeta*x
                        temp[2][iTemp] = temp_xi[idx].x_zeta_y; // x_zeta*y
                        temp[3][iTemp] = -temp_xi[idx].y_eta_z; // -y_eta*z
                        temp[4][iTemp] = -temp_xi[idx].z_eta_x; // -z_eta*x
                        temp[5][iTemp] = -temp_xi[idx].x_eta_y; // -x_eta*y
                    }
                    coef_eta[0] += NodeDifferece6th(temp[0]);
                    coef_eta[1] += NodeDifferece6th(temp[1]);
                    coef_eta[2] += NodeDifferece6th(temp[2]);
                    coef_zeta[0] += NodeDifferece6th(temp[3]);
                    coef_zeta[1] += NodeDifferece6th(temp[4]);
                    coef_zeta[2] += NodeDifferece6th(temp[5]);
                }
            }
        }
        // j direction
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 3; j < nj - 3; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    auto coef_xi = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, k));
                    auto coef_eta = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, k));
                    auto coef_zeta = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        int idx = Idx(i, j - 3 + iTemp, k);
                        temp[0][iTemp] = -temp_eta[idx].y_zeta_z; // -y_zeta*z
                        temp[1][iTemp] = -temp_eta[idx].z_zeta_x; // -z_zeta*x
                        temp[2][iTemp] = -temp_eta[idx].x_zeta_y; // -x_zeta*y
                        temp[3][iTemp] = temp_eta[idx].y_xi_z;    // y_xi*z
                        temp[4][iTemp] = temp_eta[idx].z_xi_x;    // z_xi*x
                        temp[5][iTemp] = temp_eta[idx].x_xi_y;    // x_xi*y
                    }
                    coef_xi[0] += NodeDifferece6th(temp[0]);
                    coef_xi[1] += NodeDifferece6th(temp[1]);
                    coef_xi[2] += NodeDifferece6th(temp[2]);
                    coef_zeta[0] += NodeDifferece6th(temp[3]);
                    coef_zeta[1] += NodeDifferece6th(temp[4]);
                    coef_zeta[2] += NodeDifferece6th(temp[5]);
                }
            }
        }
        // k direction
        for (int k = 3; k < nk - 3; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    auto coef_xi = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, k));
                    auto coef_eta = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, k));
                    auto coef_zeta = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        int idx = Idx(i, j, k - 3 + iTemp);
                        temp[0][iTemp] = temp_zeta[idx].y_eta_z; // y_eta*z
                        temp[1][iTemp] = temp_zeta[idx].z_eta_x; // z_eta*x
                        temp[2][iTemp] = temp_zeta[idx].x_eta_y; // x_eta*y
                        temp[3][iTemp] = -temp_zeta[idx].y_xi_z; // -y_xi*z
                        temp[4][iTemp] = -temp_zeta[idx].z_xi_x; // -z_xi*x
                        temp[5][iTemp] = -temp_zeta[idx].x_xi_y; // -x_xi*y
                    }
                    coef_xi[0] += NodeDifferece6th(temp[0]);
                    coef_xi[1] += NodeDifferece6th(temp[1]);
                    coef_xi[2] += NodeDifferece6th(temp[2]);
                    coef_eta[0] += NodeDifferece6th(temp[3]);
                    coef_eta[1] += NodeDifferece6th(temp[4]);
                    coef_eta[2] += NodeDifferece6th(temp[5]);
                }
            }
        }
        // 第六步：计算边界点度量系数
        // i=1,2;ni-3,ni-2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                auto coef_xi1 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(1, j, k));
                auto coef_eta1 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(1, j, k));
                auto coef_zeta1 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(1, j, k));
                auto coef_xi2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(2, j, k));
                auto coef_eta2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(2, j, k));
                auto coef_zeta2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(2, j, k));
                for (int iTemp = 0; iTemp < 5; iTemp++)
                {
                    int idx = Idx(0 + iTemp, j, k);
                    temp[0][iTemp] = temp_xi[idx].y_zeta_z; // y_zeta*z
                    temp[1][iTemp] = temp_xi[idx].z_zeta_x; // z_zeta*x
                    temp[2][iTemp] = temp_xi[idx].x_zeta_y; // x_zeta*y
                    temp[3][iTemp] = -temp_xi[idx].y_eta_z; // -y_eta*z
                    temp[4][iTemp] = -temp_xi[idx].z_eta_x; // -z_eta*x
                    temp[5][iTemp] = -temp_xi[idx].x_eta_y; // -x_eta*y
                }
                coef_eta1[0] += NodeDifferece4thRight(temp[0]);
                coef_eta1[1] += NodeDifferece4thRight(temp[1]);
                coef_eta1[2] += NodeDifferece4thRight(temp[2]);
                coef_zeta1[0] += NodeDifferece4thRight(temp[3]);
                coef_zeta1[1] += NodeDifferece4thRight(temp[4]);
                coef_zeta1[2] += NodeDifferece4thRight(temp[5]);
                coef_eta2[0] += NodeDifferece4th(temp[0]);
                coef_eta2[1] += NodeDifferece4th(temp[1]);
                coef_eta2[2] += NodeDifferece4th(temp[2]);
                coef_zeta2[0] += NodeDifferece4th(temp[3]);
                coef_zeta2[1] += NodeDifferece4th(temp[4]);
                coef_zeta2[2] = NodeDifferece4th(temp[5]);
                auto coef_xi_n2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(ni - 2, j, k));
                auto coef_eta_n2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(ni - 2, j, k));
                auto coef_zeta_n2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(ni - 2, j, k));
                auto coef_xi_n3 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(ni - 3, j, k));
                auto coef_eta_n3 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(ni - 3, j, k));
                auto coef_zeta_n3 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(ni - 3, j, k));
                for (int iTemp = 0; iTemp < 5; iTemp++)
                {
                    int idx = Idx(ni - 2 - iTemp, j, k);
                    temp[0][iTemp] = temp_xi[idx].y_zeta_z; // y_zeta*z
                    temp[1][iTemp] = temp_xi[idx].z_zeta_x; // z_zeta*x
                    temp[2][iTemp] = temp_xi[idx].x_zeta_y; // x_zeta*y
                    temp[3][iTemp] = -temp_xi[idx].y_eta_z; // -y_eta*z
                    temp[4][iTemp] = -temp_xi[idx].z_eta_x; // -z_eta*x
                    temp[5][iTemp] = -temp_xi[idx].x_eta_y; // -x_eta*y
                }
                coef_eta_n2[0] += NodeDifferece4thLeft(temp[0]);
                coef_eta_n2[1] += NodeDifferece4thLeft(temp[1]);
                coef_eta_n2[2] += NodeDifferece4thLeft(temp[2]);
                coef_zeta_n2[0] += NodeDifferece4thLeft(temp[3]);
                coef_zeta_n2[1] += NodeDifferece4thLeft(temp[4]);
                coef_zeta_n2[2] += NodeDifferece4thLeft(temp[5]);
                coef_eta_n3[0] += -NodeDifferece4th(temp[0]);
                coef_eta_n3[1] += -NodeDifferece4th(temp[1]);
                coef_eta_n3[2] += -NodeDifferece4th(temp[2]);
                coef_zeta_n3[0] += -NodeDifferece4th(temp[3]);
                coef_zeta_n3[1] += -NodeDifferece4th(temp[4]);
                coef_zeta_n3[2] += -NodeDifferece4th(temp[5]);
            }
        }
        // j=1,2;nj-3,nj-2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                auto coef_xi1 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, 1, k));
                auto coef_eta1 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, 1, k));
                auto coef_zeta1 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, 1, k));
                auto coef_xi2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, 2, k));
                auto coef_eta2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, 2, k));
                auto coef_zeta2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, 2, k));
                for (int iTemp = 0; iTemp < 5; iTemp++)
                {
                    int idx = Idx(i, 0 + iTemp, k);
                    temp[0][iTemp] = -temp_eta[idx].y_zeta_z; // -y_zeta*z
                    temp[1][iTemp] = -temp_eta[idx].z_zeta_x; // -z_zeta*x
                    temp[2][iTemp] = -temp_eta[idx].x_zeta_y; // -x_zeta*y
                    temp[3][iTemp] = temp_eta[idx].y_xi_z;    // y_xi*z
                    temp[4][iTemp] = temp_eta[idx].z_xi_x;    // z_xi*x
                    temp[5][iTemp] = temp_eta[idx].x_xi_y;    // x_xi*y
                }
                coef_xi1[0] += NodeDifferece4thRight(temp[0]);
                coef_xi1[1] += NodeDifferece4thRight(temp[1]);
                coef_xi1[2] += NodeDifferece4thRight(temp[2]);
                coef_zeta1[0] += NodeDifferece4thRight(temp[3]);
                coef_zeta1[1] += NodeDifferece4thRight(temp[4]);
                coef_zeta1[2] += NodeDifferece4thRight(temp[5]);
                coef_xi2[0] += NodeDifferece4th(temp[0]);
                coef_xi2[1] += NodeDifferece4th(temp[1]);
                coef_xi2[2] += NodeDifferece4th(temp[2]);
                coef_zeta2[0] += NodeDifferece4th(temp[3]);
                coef_zeta2[1] += NodeDifferece4th(temp[4]);
                coef_zeta2[2] += NodeDifferece4th(temp[5]);
                auto coef_xi_n2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, nj - 2, k));
                auto coef_eta_n2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, nj - 2, k));
                auto coef_zeta_n2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, nj - 2, k));
                auto coef_xi_n3 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, nj - 3, k));
                auto coef_eta_n3 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, nj - 3, k));
                auto coef_zeta_n3 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, nj - 3, k));
                for (int iTemp = 0; iTemp < 5; iTemp++)
                {
                    int idx = Idx(i, nj - 2 - iTemp, k);
                    temp[0][iTemp] = -temp_eta[idx].y_zeta_z; // -y_zeta*z
                    temp[1][iTemp] = -temp_eta[idx].z_zeta_x; // -z_zeta*x
                    temp[2][iTemp] = -temp_eta[idx].x_zeta_y; // -x_zeta*y
                    temp[3][iTemp] = temp_eta[idx].y_xi_z;    // y_xi*z
                    temp[4][iTemp] = temp_eta[idx].z_xi_x;    // z_xi*x
                    temp[5][iTemp] = temp_eta[idx].x_xi_y;    // x_xi*y
                }
                coef_xi_n2[0] += NodeDifferece4thLeft(temp[0]);
                coef_xi_n2[1] += NodeDifferece4thLeft(temp[1]);
                coef_xi_n2[2] += NodeDifferece4thLeft(temp[2]);
                coef_zeta_n2[0] += NodeDifferece4thLeft(temp[3]);
                coef_zeta_n2[1] += NodeDifferece4thLeft(temp[4]);
                coef_zeta_n2[2] += NodeDifferece4thLeft(temp[5]);
                coef_xi_n3[0] += -NodeDifferece4th(temp[0]);
                coef_xi_n3[1] += -NodeDifferece4th(temp[1]);
                coef_xi_n3[2] += -NodeDifferece4th(temp[2]);
                coef_zeta_n3[0] += -NodeDifferece4th(temp[3]);
                coef_zeta_n3[1] += -NodeDifferece4th(temp[4]);
                coef_zeta_n3[2] += -NodeDifferece4th(temp[5]);
            }
        }
        // k=1,2;nk-3,nk-2
        for (int j = 1; j < nj - 1; ++j)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                auto coef_xi1 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, 1));
                auto coef_eta1 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, 1));
                auto coef_zeta1 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, 1));
                auto coef_xi2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, 2));
                auto coef_eta2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, 2));
                auto coef_zeta2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, 2));
                for (int iTemp = 0; iTemp < 5; iTemp++)
                {
                    int idx = Idx(i, j, 0 + iTemp);
                    temp[0][iTemp] = temp_zeta[idx].y_eta_z; // y_eta*z
                    temp[1][iTemp] = temp_zeta[idx].z_eta_x; // z_eta*x
                    temp[2][iTemp] = temp_zeta[idx].x_eta_y; // x_eta*y
                    temp[3][iTemp] = -temp_zeta[idx].y_xi_z; // -y_xi*z
                    temp[4][iTemp] = -temp_zeta[idx].z_xi_x; // -z_xi*x
                    temp[5][iTemp] = -temp_zeta[idx].x_xi_y; // -x_xi*y
                }
                coef_xi1[0] += NodeDifferece4thRight(temp[0]);
                coef_xi1[1] += NodeDifferece4thRight(temp[1]);
                coef_xi1[2] += NodeDifferece4thRight(temp[2]);
                coef_eta1[0] += NodeDifferece4thRight(temp[3]);
                coef_eta1[1] += NodeDifferece4thRight(temp[4]);
                coef_eta1[2] += NodeDifferece4thRight(temp[5]);
                coef_xi2[0] += NodeDifferece4th(temp[0]);
                coef_xi2[1] += NodeDifferece4th(temp[1]);
                coef_xi2[2] += NodeDifferece4th(temp[2]);
                coef_eta2[0] += NodeDifferece4th(temp[3]);
                coef_eta2[1] += NodeDifferece4th(temp[4]);
                coef_eta2[2] += NodeDifferece4th(temp[5]);
                auto coef_xi_n2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, nk - 2));
                auto coef_eta_n2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, nk - 2));
                auto coef_zeta_n2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, nk - 2));
                auto coef_xi_n3 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, nk - 3));
                auto coef_eta_n3 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, nk - 3));
                auto coef_zeta_n3 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, nk - 3));
                for (int iTemp = 0; iTemp < 5; iTemp++)
                {
                    int idx = Idx(i, j, nk - 2 - iTemp);
                    temp[0][iTemp] = temp_zeta[idx].y_eta_z; // y_eta*z
                    temp[1][iTemp] = temp_zeta[idx].z_eta_x; // z_eta*x
                    temp[2][iTemp] = temp_zeta[idx].x_eta_y; // x_eta*y
                    temp[3][iTemp] = -temp_zeta[idx].y_xi_z; // -y_xi*z
                    temp[4][iTemp] = -temp_zeta[idx].z_xi_x; // -z_xi*x
                    temp[5][iTemp] = -temp_zeta[idx].x_xi_y; // -x_xi*y
                }
                coef_xi_n2[0] += NodeDifferece4thLeft(temp[0]);
                coef_xi_n2[1] += NodeDifferece4thLeft(temp[1]);
                coef_xi_n2[2] += NodeDifferece4thLeft(temp[2]);
                coef_eta_n2[0] += NodeDifferece4thLeft(temp[3]);
                coef_eta_n2[1] += NodeDifferece4thLeft(temp[4]);
                coef_eta_n2[2] += NodeDifferece4thLeft(temp[5]);
                coef_xi_n3[0] += -NodeDifferece4th(temp[0]);
                coef_xi_n3[1] += -NodeDifferece4th(temp[1]);
                coef_xi_n3[2] += -NodeDifferece4th(temp[2]);
                coef_eta_n3[0] += -NodeDifferece4th(temp[3]);
                coef_eta_n3[1] += -NodeDifferece4th(temp[4]);
                coef_eta_n3[2] += -NodeDifferece4th(temp[5]);
            }
        }
    }
    void NSSolverStruct::CalcMetricsS2()
    {
        auto para = GetPara();
        auto flux_diff_scheme = para->GetDifferenceScheme();
        if (flux_diff_scheme == DifferenceScheme::SecondOrder)
        {
            CalcMetricsS2_2nd();
        }
        else if (flux_diff_scheme == DifferenceScheme::SixthOrder)
        {
            CalcMetricsS2_6th();
        }
        else
        {
            Log::warn("FluxDifferenceScheme is not defined, use SecondOrder as default");
            CalcMetricsS2_2nd();
        }
    }
    void NSSolverStruct::CalcMetricsS2_2nd()
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
        int idx_temp[3];
        // i+1/2,j+1/2,k+1/2处的坐标
        std::vector<std::vector<double>> coord_i(ni * nj * nk, std::vector<double>(3)),
            coord_j(ni * nj * nk, std::vector<double>(3)), coord_k(ni * nj * nk, std::vector<double>(3));
        // 计算之前，先把度量系数赋值为0
        for (int i = 0; i < ni; ++i)
        {
            for (int j = 0; j < nj; ++j)
            {
                for (int k = 0; k < nk; ++k)
                {
                    for (int iDim = 0; iDim < 4; ++iDim)
                    {
                        coef->GetX(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetY(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetZ(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetXi(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetEta(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetZeta(Idx(i, j, k))[iDim] = 0.0;
                        for (int jDim = 0; jDim < 3; ++jDim)
                        {
                            coef_mid[jDim]->GetX(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetY(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetZ(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetXi(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetEta(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetZeta(Idx(i, j, k))[iDim] = 0.0;
                        }
                    }
                }
            }
        }
        // 第一步：计算半点坐标
        for (int k = 0; k < nk - 1; ++k)
        {
            for (int j = 0; j < nj - 1; ++j)
            {
                for (int i = 0; i < ni - 1; ++i)
                {
                    for (int iDim = 0; iDim < 3; ++iDim)
                    {
                        coord_i[Idx(i, j, k)][iDim] =
                            0.5 * (node->GetCoord(i, j, k)[iDim] + node->GetCoord(i + 1, j, k)[iDim]);
                        coord_j[Idx(i, j, k)][iDim] =
                            0.5 * (node->GetCoord(i, j, k)[iDim] + node->GetCoord(i, j + 1, k)[iDim]);
                        coord_k[Idx(i, j, k)][iDim] =
                            0.5 * (node->GetCoord(i, j, k)[iDim] + node->GetCoord(i, j, k + 1)[iDim]);
                    }
                }
            }
        }
        // 第二步：根据半点坐标计算整点逆变换度量系数
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    int idx = idx_proxy->GetIdx(i, j, k);
                    auto coef_x = coef->GetX(idx);
                    auto coef_y = coef->GetY(idx);
                    auto coef_z = coef->GetZ(idx);
                    for (int iDim = 0; iDim < 3; ++iDim)
                    {
                        // i direction
                        coef_x[iDim] = coord_i[Idx(i, j, k)][iDim] - coord_i[Idx(i - 1, j, k)][iDim];
                        // j direction
                        coef_y[iDim] = coord_j[Idx(i, j, k)][iDim] - coord_j[Idx(i, j - 1, k)][iDim];
                        // k direction
                        coef_z[iDim] = coord_k[Idx(i, j, k)][iDim] - coord_k[Idx(i, j, k - 1)][iDim];
                    }
                    coef_x[3] = 0.0;
                    coef_y[3] = 0.0;
                    coef_z[3] = 0.0;
                    if (grid->GetDim() == 2)
                    {
                        coef_z[0] = 0.0;
                        coef_z[1] = 0.0;
                        coef_z[2] = 1.0;
                    }
                }
            }
        }
        struct TempXi
        {
            double z_zeta_y, x_zeta_z, y_zeta_x, z_eta_y, x_eta_z, y_eta_x;
        };
        std::vector<TempXi> temp_xi(ni * nj * nk);
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    int idx = Idx(i, j, k);
                    auto &temp = temp_xi[idx];
                    temp.z_zeta_y = 0.5 * (coef->GetZ(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[1] +
                                           coef->GetZ(Idx(i + 1, j, k))[2] * node->GetCoord(i + 1, j, k)[1]);
                    temp.x_zeta_z = 0.5 * (coef->GetX(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[2] +
                                           coef->GetX(Idx(i + 1, j, k))[2] * node->GetCoord(i + 1, j, k)[2]);
                    temp.y_zeta_x = 0.5 * (coef->GetY(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[0] +
                                           coef->GetY(Idx(i + 1, j, k))[2] * node->GetCoord(i + 1, j, k)[0]);
                    temp.z_eta_y = 0.5 * (coef->GetZ(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[1] +
                                          coef->GetZ(Idx(i + 1, j, k))[1] * node->GetCoord(i + 1, j, k)[1]);
                    temp.x_eta_z = 0.5 * (coef->GetX(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[2] +
                                          coef->GetX(Idx(i + 1, j, k))[1] * node->GetCoord(i + 1, j, k)[2]);
                    temp.y_eta_x = 0.5 * (coef->GetY(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[0] +
                                          coef->GetY(Idx(i + 1, j, k))[1] * node->GetCoord(i + 1, j, k)[0]);
                }
            }
        }
        // i=1/2;ni-3/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                temp_xi[Idx(0, j, k)].z_zeta_y = 0.5 * (3.0 * coef->GetZ(Idx(1, j, k))[2] * node->GetCoord(1, j, k)[1] - 1.0 * coef->GetZ(Idx(2, j, k))[2] * node->GetCoord(2, j, k)[1]);
                temp_xi[Idx(0, j, k)].x_zeta_z = 0.5 * (3.0 * coef->GetX(Idx(1, j, k))[2] * node->GetCoord(1, j, k)[2] - 1.0 * coef->GetX(Idx(2, j, k))[2] * node->GetCoord(2, j, k)[2]);
                temp_xi[Idx(0, j, k)].y_zeta_x = 0.5 * (3.0 * coef->GetY(Idx(1, j, k))[2] * node->GetCoord(1, j, k)[0] - 1.0 * coef->GetY(Idx(2, j, k))[2] * node->GetCoord(2, j, k)[0]);
                temp_xi[Idx(0, j, k)].z_eta_y = 0.5 * (3.0 * coef->GetZ(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[1] - 1.0 * coef->GetZ(Idx(2, j, k))[1] * node->GetCoord(2, j, k)[1]);
                temp_xi[Idx(0, j, k)].x_eta_z = 0.5 * (3.0 * coef->GetX(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[2] - 1.0 * coef->GetX(Idx(2, j, k))[1] * node->GetCoord(2, j, k)[2]);
                temp_xi[Idx(0, j, k)].y_eta_x = 0.5 * (3.0 * coef->GetY(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[0] - 1.0 * coef->GetY(Idx(2, j, k))[1] * node->GetCoord(2, j, k)[0]);
                temp_xi[Idx(ni - 2, j, k)].z_zeta_y = 0.5 * (3.0 * coef->GetZ(Idx(ni - 2, j, k))[2] * node->GetCoord(ni - 2, j, k)[1] - 1.0 * coef->GetZ(Idx(ni - 3, j, k))[2] * node->GetCoord(ni - 3, j, k)[1]);
                temp_xi[Idx(ni - 2, j, k)].x_zeta_z = 0.5 * (3.0 * coef->GetX(Idx(ni - 2, j, k))[2] * node->GetCoord(ni - 2, j, k)[2] - 1.0 * coef->GetX(Idx(ni - 3, j, k))[2] * node->GetCoord(ni - 3, j, k)[2]);
                temp_xi[Idx(ni - 2, j, k)].y_zeta_x = 0.5 * (3.0 * coef->GetY(Idx(ni - 2, j, k))[2] * node->GetCoord(ni - 2, j, k)[0] - 1.0 * coef->GetY(Idx(ni - 3, j, k))[2] * node->GetCoord(ni - 3, j, k)[0]);
                temp_xi[Idx(ni - 2, j, k)].z_eta_y = 0.5 * (3.0 * coef->GetZ(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[1] - 1.0 * coef->GetZ(Idx(ni - 3, j, k))[1] * node->GetCoord(ni - 3, j, k)[1]);
                temp_xi[Idx(ni - 2, j, k)].x_eta_z = 0.5 * (3.0 * coef->GetX(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[2] - 1.0 * coef->GetX(Idx(ni - 3, j, k))[1] * node->GetCoord(ni - 3, j, k)[2]);
                temp_xi[Idx(ni - 2, j, k)].y_eta_x = 0.5 * (3.0 * coef->GetY(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[0] - 1.0 * coef->GetY(Idx(ni - 3, j, k))[1] * node->GetCoord(ni - 3, j, k)[0]);
            }
        }
        // 求出eta方向的临时变量：逆变换度量系数乘以坐标
        struct TempEta
        {
            double z_zeta_y, x_zeta_z, y_zeta_x, z_xi_y, x_xi_z, y_xi_x;
        };
        std::vector<TempEta> temp_eta(ni * nj * nk);
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    int idx = Idx(i, j, k);
                    auto &temp = temp_eta[idx];
                    temp.z_zeta_y = 0.5 * (coef->GetZ(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[1] +
                                           coef->GetZ(Idx(i, j + 1, k))[2] * node->GetCoord(i, j + 1, k)[1]);
                    temp.x_zeta_z = 0.5 * (coef->GetX(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[2] +
                                           coef->GetX(Idx(i, j + 1, k))[2] * node->GetCoord(i, j + 1, k)[2]);
                    temp.y_zeta_x = 0.5 * (coef->GetY(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[0] +
                                           coef->GetY(Idx(i, j + 1, k))[2] * node->GetCoord(i, j + 1, k)[0]);
                    temp.z_xi_y = 0.5 * (coef->GetZ(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[1] +
                                         coef->GetZ(Idx(i, j + 1, k))[0] * node->GetCoord(i, j + 1, k)[1]);
                    temp.x_xi_z = 0.5 * (coef->GetX(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[2] +
                                         coef->GetX(Idx(i, j + 1, k))[0] * node->GetCoord(i, j + 1, k)[2]);
                    temp.y_xi_x = 0.5 * (coef->GetY(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[0] +
                                         coef->GetY(Idx(i, j + 1, k))[0] * node->GetCoord(i, j + 1, k)[0]);
                }
            }
        }
        // j=1/2;nj-3/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                temp_eta[Idx(i, 0, k)].z_zeta_y = 0.5 * (3.0 * coef->GetZ(Idx(i, 1, k))[2] * node->GetCoord(i, 1, k)[1] - 1.0 * coef->GetZ(Idx(i, 2, k))[2] * node->GetCoord(i, 2, k)[1]);
                temp_eta[Idx(i, 0, k)].x_zeta_z = 0.5 * (3.0 * coef->GetX(Idx(i, 1, k))[2] * node->GetCoord(i, 1, k)[2] - 1.0 * coef->GetX(Idx(i, 2, k))[2] * node->GetCoord(i, 2, k)[2]);
                temp_eta[Idx(i, 0, k)].y_zeta_x = 0.5 * (3.0 * coef->GetY(Idx(i, 1, k))[2] * node->GetCoord(i, 1, k)[0] - 1.0 * coef->GetY(Idx(i, 2, k))[2] * node->GetCoord(i, 2, k)[0]);
                temp_eta[Idx(i, 0, k)].z_xi_y = 0.5 * (3.0 * coef->GetZ(Idx(i, 1, k))[0] * node->GetCoord(i, 1, k)[1] - 1.0 * coef->GetZ(Idx(i, 2, k))[0] * node->GetCoord(i, 2, k)[1]);
                temp_eta[Idx(i, 0, k)].x_xi_z = 0.5 * (3.0 * coef->GetX(Idx(i, 1, k))[0] * node->GetCoord(i, 1, k)[2] - 1.0 * coef->GetX(Idx(i, 2, k))[0] * node->GetCoord(i, 2, k)[2]);
                temp_eta[Idx(i, 0, k)].y_xi_x = 0.5 * (3.0 * coef->GetY(Idx(i, 1, k))[0] * node->GetCoord(i, 1, k)[0] - 1.0 * coef->GetY(Idx(i, 2, k))[0] * node->GetCoord(i, 2, k)[0]);
                temp_eta[Idx(i, nj - 2, k)].z_zeta_y = 0.5 * (3.0 * coef->GetZ(Idx(i, nj - 2, k))[2] * node->GetCoord(i, nj - 2, k)[1] - 1.0 * coef->GetZ(Idx(i, nj - 3, k))[2] * node->GetCoord(i, nj - 3, k)[1]);
                temp_eta[Idx(i, nj - 2, k)].x_zeta_z = 0.5 * (3.0 * coef->GetX(Idx(i, nj - 2, k))[2] * node->GetCoord(i, nj - 2, k)[2] - 1.0 * coef->GetX(Idx(i, nj - 3, k))[2] * node->GetCoord(i, nj - 3, k)[2]);
                temp_eta[Idx(i, nj - 2, k)].y_zeta_x = 0.5 * (3.0 * coef->GetY(Idx(i, nj - 2, k))[2] * node->GetCoord(i, nj - 2, k)[0] - 1.0 * coef->GetY(Idx(i, nj - 3, k))[2] * node->GetCoord(i, nj - 3, k)[0]);
                temp_eta[Idx(i, nj - 2, k)].z_xi_y = 0.5 * (3.0 * coef->GetZ(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[1] - 1.0 * coef->GetZ(Idx(i, nj - 3, k))[0] * node->GetCoord(i, nj - 3, k)[1]);
                temp_eta[Idx(i, nj - 2, k)].x_xi_z = 0.5 * (3.0 * coef->GetX(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[2] - 1.0 * coef->GetX(Idx(i, nj - 3, k))[0] * node->GetCoord(i, nj - 3, k)[2]);
                temp_eta[Idx(i, nj - 2, k)].y_xi_x = 0.5 * (3.0 * coef->GetY(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[0] - 1.0 * coef->GetY(Idx(i, nj - 3, k))[0] * node->GetCoord(i, nj - 3, k)[0]);
            }
        }
        // 求出zeta方向的临时变量：逆变换度量系数乘以坐标
        struct TempZeta
        {
            double z_eta_y, x_eta_z, y_eta_x, z_xi_y, x_xi_z, y_xi_x;
        };
        std::vector<TempZeta> temp_zeta(ni * nj * nk);
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    int idx = Idx(i, j, k);
                    auto &temp = temp_zeta[idx];
                    temp.z_eta_y = 0.5 * (coef->GetZ(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[1] +
                                          coef->GetZ(Idx(i, j, k + 1))[1] * node->GetCoord(i, j, k + 1)[1]);
                    temp.x_eta_z = 0.5 * (coef->GetX(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[2] +
                                          coef->GetX(Idx(i, j, k + 1))[1] * node->GetCoord(i, j, k + 1)[2]);
                    temp.y_eta_x = 0.5 * (coef->GetY(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[0] +
                                          coef->GetY(Idx(i, j, k + 1))[1] * node->GetCoord(i, j, k + 1)[0]);
                    temp.z_xi_y = 0.5 * (coef->GetZ(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[1] +
                                         coef->GetZ(Idx(i, j, k + 1))[0] * node->GetCoord(i, j, k + 1)[1]);
                    temp.x_xi_z = 0.5 * (coef->GetX(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[2] +
                                         coef->GetX(Idx(i, j, k + 1))[0] * node->GetCoord(i, j, k + 1)[2]);
                    temp.y_xi_x = 0.5 * (coef->GetY(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[0] +
                                         coef->GetY(Idx(i, j, k + 1))[0] * node->GetCoord(i, j, k + 1)[0]);
                }
            }
        }
        // k=1/2;nk-3/2
        for (int j = 1; j < nj - 1; ++j)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                temp_zeta[Idx(i, j, 0)].z_eta_y = 0.5 * (3.0 * coef->GetZ(Idx(i, j, 1))[1] * node->GetCoord(i, j, 1)[1] - 1.0 * coef->GetZ(Idx(i, j, 2))[1] * node->GetCoord(i, j, 2)[1]);
                temp_zeta[Idx(i, j, 0)].x_eta_z = 0.5 * (3.0 * coef->GetX(Idx(i, j, 1))[1] * node->GetCoord(i, j, 1)[2] - 1.0 * coef->GetX(Idx(i, j, 2))[1] * node->GetCoord(i, j, 2)[2]);
                temp_zeta[Idx(i, j, 0)].y_eta_x = 0.5 * (3.0 * coef->GetY(Idx(i, j, 1))[1] * node->GetCoord(i, j, 1)[0] - 1.0 * coef->GetY(Idx(i, j, 2))[1] * node->GetCoord(i, j, 2)[0]);
                temp_zeta[Idx(i, j, 0)].z_xi_y = 0.5 * (3.0 * coef->GetZ(Idx(i, j, 1))[0] * node->GetCoord(i, j, 1)[1] - 1.0 * coef->GetZ(Idx(i, j, 2))[0] * node->GetCoord(i, j, 2)[1]);
                temp_zeta[Idx(i, j, 0)].x_xi_z = 0.5 * (3.0 * coef->GetX(Idx(i, j, 1))[0] * node->GetCoord(i, j, 1)[2] - 1.0 * coef->GetX(Idx(i, j, 2))[0] * node->GetCoord(i, j, 2)[2]);
                temp_zeta[Idx(i, j, 0)].y_xi_x = 0.5 * (3.0 * coef->GetY(Idx(i, j, 1))[0] * node->GetCoord(i, j, 1)[0] - 1.0 * coef->GetY(Idx(i, j, 2))[0] * node->GetCoord(i, j, 2)[0]);
                temp_zeta[Idx(i, j, nk - 2)].z_eta_y = 0.5 * (3.0 * coef->GetZ(Idx(i, j, nk - 2))[1] * node->GetCoord(i, j, nk - 2)[1] - 1.0 * coef->GetZ(Idx(i, j, nk - 3))[1] * node->GetCoord(i, j, nk - 3)[1]);
                temp_zeta[Idx(i, j, nk - 2)].x_eta_z = 0.5 * (3.0 * coef->GetX(Idx(i, j, nk - 2))[1] * node->GetCoord(i, j, nk - 2)[2] - 1.0 * coef->GetX(Idx(i, j, nk - 3))[1] * node->GetCoord(i, j, nk - 3)[2]);
                temp_zeta[Idx(i, j, nk - 2)].y_eta_x = 0.5 * (3.0 * coef->GetY(Idx(i, j, nk - 2))[1] * node->GetCoord(i, j, nk - 2)[0] - 1.0 * coef->GetY(Idx(i, j, nk - 3))[1] * node->GetCoord(i, j, nk - 3)[0]);
                temp_zeta[Idx(i, j, nk - 2)].z_xi_y = 0.5 * (3.0 * coef->GetZ(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[1] - 1.0 * coef->GetZ(Idx(i, j, nk - 3))[0] * node->GetCoord(i, j, nk - 3)[1]);
                temp_zeta[Idx(i, j, nk - 2)].x_xi_z = 0.5 * (3.0 * coef->GetX(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[2] - 1.0 * coef->GetX(Idx(i, j, nk - 3))[0] * node->GetCoord(i, j, nk - 3)[2]);
                temp_zeta[Idx(i, j, nk - 2)].y_xi_x = 0.5 * (3.0 * coef->GetY(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[0] - 1.0 * coef->GetY(Idx(i, j, nk - 3))[0] * node->GetCoord(i, j, nk - 3)[0]);
            }
        }
        // 第五步：根据半点坐标和半点逆变换度量系数使用守恒形式计算整点度量系数（CMM2)
        double temp[6][2];
        // 获取下标的lamda函数
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    // i direction
                    auto coef_xi = coef->GetXi(Idx(i, j, k));
                    auto coef_eta = coef->GetEta(Idx(i, j, k));
                    auto coef_zeta = coef->GetZeta(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 2; iTemp++)
                    {
                        int idx = Idx(i - 1 + iTemp, j, k);
                        temp[0][iTemp] = -temp_xi[idx].z_zeta_y; // -z_zeta*y
                        temp[1][iTemp] = -temp_xi[idx].x_zeta_z; // -x_zeta*z
                        temp[2][iTemp] = -temp_xi[idx].y_zeta_x; // -y_zeta*x
                        temp[3][iTemp] = temp_xi[idx].z_eta_y;   // z_eta*y
                        temp[4][iTemp] = temp_xi[idx].x_eta_z;   // x_eta*z
                        temp[5][iTemp] = temp_xi[idx].y_eta_x;   // y_eta*x
                    }
                    coef_eta[0] += temp[0][1] - temp[0][0];
                    coef_eta[1] += temp[1][1] - temp[1][0];
                    coef_eta[2] += temp[2][1] - temp[2][0];
                    coef_zeta[0] += temp[3][1] - temp[3][0];
                    coef_zeta[1] += temp[4][1] - temp[4][0];
                    coef_zeta[2] += temp[5][1] - temp[5][0];
                    for (int iTemp = 0; iTemp < 2; iTemp++)
                    {
                        int idx = Idx(i, j - 1 + iTemp, k);
                        temp[0][iTemp] = temp_eta[idx].z_zeta_y; // z_zeta*y
                        temp[1][iTemp] = temp_eta[idx].x_zeta_z; // x_zeta*z
                        temp[2][iTemp] = temp_eta[idx].y_zeta_x; // y_zeta*x
                        temp[3][iTemp] = -temp_eta[idx].z_xi_y;  // -z_xi*y
                        temp[4][iTemp] = -temp_eta[idx].x_xi_z;  // -x_xi*z
                        temp[5][iTemp] = -temp_eta[idx].y_xi_x;  // -y_xi*x
                    }
                    coef_xi[0] += temp[0][1] - temp[0][0];
                    coef_xi[1] += temp[1][1] - temp[1][0];
                    coef_xi[2] += temp[2][1] - temp[2][0];
                    coef_zeta[0] += temp[3][1] - temp[3][0];
                    coef_zeta[1] += temp[4][1] - temp[4][0];
                    coef_zeta[2] += temp[5][1] - temp[5][0];
                    for (int iTemp = 0; iTemp < 2; iTemp++)
                    {
                        int idx = Idx(i, j, k - 1 + iTemp);
                        temp[0][iTemp] = -temp_zeta[idx].z_eta_y; // -z_eta*y
                        temp[1][iTemp] = -temp_zeta[idx].x_eta_z; // -x_eta*z
                        temp[2][iTemp] = -temp_zeta[idx].y_eta_x; // -y_eta*x
                        temp[3][iTemp] = temp_zeta[idx].z_xi_y;   // z_xi*y
                        temp[4][iTemp] = temp_zeta[idx].x_xi_z;   // x_xi*z
                        temp[5][iTemp] = temp_zeta[idx].y_xi_x;   // y_xi*x
                    }
                    coef_xi[0] += temp[0][1] - temp[0][0];
                    coef_xi[1] += temp[1][1] - temp[1][0];
                    coef_xi[2] += temp[2][1] - temp[2][0];
                    coef_eta[0] += temp[3][1] - temp[3][0];
                    coef_eta[1] += temp[4][1] - temp[4][0];
                    coef_eta[2] += temp[5][1] - temp[5][0];
                    if (grid->GetDim() == 2)
                    {
                        coef_zeta[0] = 0.0;
                        coef_zeta[1] = 0.0;
                        coef_zeta[2] = 1.0;
                    }
                    coef_xi[3] = coef_eta[3] = coef_zeta[3] = 0.0;
                }
            }
        }
    }
    void NSSolverStruct::CalcMetricsS2_6th()
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
        double inter_temp[6];
        // i+1/2,j+1/2,k+1/2处的坐标
        std::vector<std::vector<double>> coord_i(ni * nj * nk, std::vector<double>(3)),
            coord_j(ni * nj * nk, std::vector<double>(3)), coord_k(ni * nj * nk, std::vector<double>(3));
        // 计算之前，先把度量系数赋值为0
        for (int i = 0; i < ni; ++i)
        {
            for (int j = 0; j < nj; ++j)
            {
                for (int k = 0; k < nk; ++k)
                {
                    for (int iDim = 0; iDim < 4; ++iDim)
                    {
                        coef->GetX(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetY(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetZ(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetXi(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetEta(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetZeta(Idx(i, j, k))[iDim] = 0.0;
                        for (int jDim = 0; jDim < 3; ++jDim)
                        {
                            coef_mid[jDim]->GetX(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetY(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetZ(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetXi(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetEta(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetZeta(Idx(i, j, k))[iDim] = 0.0;
                        }
                    }
                }
            }
        }
        // 第一步：计算半点坐标
        // i+1/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 3; i < ni - 4; ++i)
                {
                    for (int iDim = 0; iDim < 3; ++iDim)
                    {
                        for (int iTemp = 0; iTemp < 6; iTemp++)
                        {
                            inter_temp[iTemp] = node->GetCoord(i + iTemp - 2, j, k)[iDim];
                        }
                        coord_i[Idx(i, j, k)][iDim] = MidNodeInter6th(inter_temp);
                    }
                }
            }
        }
        // i=1/2,3/2,5/2;ni-3/2,ni-5/2,ni-7/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int iDim = 0; iDim < 3; ++iDim)
                {
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(1 + iTemp, j, k)[iDim];
                    }
                    coord_i[Idx(0, j, k)][iDim] = MidNodeInter4thRight1(inter_temp);
                    coord_i[Idx(1, j, k)][iDim] = MidNodeInter4thRight2(inter_temp);
                    coord_i[Idx(2, j, k)][iDim] = MidNodeInter4th(inter_temp);
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(ni - 5 + iTemp, j, k)[iDim];
                    }
                    coord_i[Idx(ni - 2, j, k)][iDim] = MidNodeInter4thLeft1(inter_temp);
                    coord_i[Idx(ni - 3, j, k)][iDim] = MidNodeInter4thLeft2(inter_temp);
                    coord_i[Idx(ni - 4, j, k)][iDim] = MidNodeInter4th(inter_temp);
                }
            }
        }
        // j+1/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 3; j < nj - 4; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    for (int iDim = 0; iDim < 3; ++iDim)
                    {
                        for (int iTemp = 0; iTemp < 6; iTemp++)
                        {
                            inter_temp[iTemp] = node->GetCoord(i, j + iTemp - 2, k)[iDim];
                        }
                        coord_j[Idx(i, j, k)][iDim] = MidNodeInter6th(inter_temp);
                    }
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
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(i, 1 + iTemp, k)[iDim];
                    }
                    coord_j[Idx(i, 0, k)][iDim] = MidNodeInter4thRight1(inter_temp);
                    coord_j[Idx(i, 1, k)][iDim] = MidNodeInter4thRight2(inter_temp);
                    coord_j[Idx(i, 2, k)][iDim] = MidNodeInter4th(inter_temp);
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(i, nj - 5 + iTemp, k)[iDim];
                    }
                    coord_j[Idx(i, nj - 2, k)][iDim] = MidNodeInter4thLeft1(inter_temp);
                    coord_j[Idx(i, nj - 3, k)][iDim] = MidNodeInter4thLeft2(inter_temp);
                    coord_j[Idx(i, nj - 4, k)][iDim] = MidNodeInter4th(inter_temp);
                }
            }
        }
        // k+1/2
        for (int k = 3; k < nk - 4; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    for (int iDim = 0; iDim < 3; ++iDim)
                    {
                        for (int iTemp = 0; iTemp < 6; iTemp++)
                        {
                            inter_temp[iTemp] = node->GetCoord(i, j, k + iTemp - 2)[iDim];
                        }
                        coord_k[Idx(i, j, k)][iDim] = MidNodeInter6th(inter_temp);
                    }
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
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(i, j, 1 + iTemp)[iDim];
                    }
                    coord_k[Idx(i, j, 0)][iDim] = MidNodeInter4thRight1(inter_temp);
                    coord_k[Idx(i, j, 1)][iDim] = MidNodeInter4thRight2(inter_temp);
                    coord_k[Idx(i, j, 2)][iDim] = MidNodeInter4th(inter_temp);
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(i, j, nk - 5 + iTemp)[iDim];
                    }
                    coord_k[Idx(i, j, nk - 2)][iDim] = MidNodeInter4thLeft1(inter_temp);
                    coord_k[Idx(i, j, nk - 3)][iDim] = MidNodeInter4thLeft2(inter_temp);
                    coord_k[Idx(i, j, nk - 4)][iDim] = MidNodeInter4th(inter_temp);
                }
            }
        }
        // 第二步：根据半点坐标计算整点逆变换度量系数
        double diff_temp[6];
        // i direction
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 3; i < ni - 3; ++i)
                {
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_i[Idx(i + iTemp - 3, j, k)][0];
                    }
                    auto coef_x = coef->GetX(Idx(i, j, k));
                    coef_x[0] = NodeDifferece6th(diff_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_i[Idx(i + iTemp - 3, j, k)][1];
                    }
                    auto coef_y = coef->GetY(Idx(i, j, k));
                    coef_y[0] = NodeDifferece6th(diff_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_i[Idx(i + iTemp - 3, j, k)][2];
                    }
                    auto coef_z = coef->GetZ(Idx(i, j, k));
                    coef_z[0] = NodeDifferece6th(diff_temp);
                }
            }
        }
        // j direction
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 3; j < nj - 3; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    auto coef_x = coef->GetX(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_j[Idx(i, j + iTemp - 3, k)][0];
                    }
                    coef_x[1] = NodeDifferece6th(diff_temp);
                    auto coef_y = coef->GetY(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_j[Idx(i, j + iTemp - 3, k)][1];
                    }
                    coef_y[1] = NodeDifferece6th(diff_temp);
                    auto coef_z = coef->GetZ(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_j[Idx(i, j + iTemp - 3, k)][2];
                    }
                    coef_z[1] = NodeDifferece6th(diff_temp);
                }
            }
        }
        // k direction
        for (int k = 3; k < nk - 3; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    auto coef_x = coef->GetX(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_k[Idx(i, j, k + iTemp - 3)][0];
                    }
                    coef_x[2] = NodeDifferece6th(diff_temp);
                    auto coef_y = coef->GetY(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_k[Idx(i, j, k + iTemp - 3)][1];
                    }
                    coef_y[2] = NodeDifferece6th(diff_temp);
                    auto coef_z = coef->GetZ(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_k[Idx(i, j, k + iTemp - 3)][2];
                    }
                    coef_z[2] = NodeDifferece6th(diff_temp);
                    if (grid->GetDim() == 2)
                    {
                        coef_x[2] = 0.0;
                        coef_y[2] = 0.0;
                        coef_z[2] = 1.0;
                    }
                }
            }
        }
        // 计算边界整点的度量系数
        // i=1,2;ni-3,ni-2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_i[Idx(0 + iTemp, j, k)][0];
                }
                coef->GetX(Idx(1, j, k))[0] = NodeDifferece4thRight(diff_temp);
                coef->GetX(Idx(2, j, k))[0] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_i[Idx(0 + iTemp, j, k)][1];
                }
                coef->GetY(Idx(1, j, k))[0] = NodeDifferece4thRight(diff_temp);
                coef->GetY(Idx(2, j, k))[0] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_i[Idx(0 + iTemp, j, k)][2];
                }
                coef->GetZ(Idx(1, j, k))[0] = NodeDifferece4thRight(diff_temp);
                coef->GetZ(Idx(2, j, k))[0] = NodeDifferece4th(diff_temp);

                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_i[Idx(ni - 2 - iTemp, j, k)][0];
                }
                coef->GetX(Idx(ni - 2, j, k))[0] = NodeDifferece4thLeft(diff_temp);
                coef->GetX(Idx(ni - 3, j, k))[0] = -NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_i[Idx(ni - 2 - iTemp, j, k)][1];
                }
                coef->GetY(Idx(ni - 2, j, k))[0] = NodeDifferece4thLeft(diff_temp);
                coef->GetY(Idx(ni - 3, j, k))[0] = -NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_i[Idx(ni - 2 - iTemp, j, k)][2];
                }
                coef->GetZ(Idx(ni - 2, j, k))[0] = NodeDifferece4thLeft(diff_temp);
                coef->GetZ(Idx(ni - 3, j, k))[0] = -NodeDifferece4th(diff_temp);
            }
        }
        // j=1,2;nj-3,nj-2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_j[Idx(i, 0 + iTemp, k)][0];
                }
                coef->GetX(Idx(i, 1, k))[1] = NodeDifferece4thRight(diff_temp);
                coef->GetX(Idx(i, 2, k))[1] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_j[Idx(i, 0 + iTemp, k)][1];
                }
                coef->GetY(Idx(i, 1, k))[1] = NodeDifferece4thRight(diff_temp);
                coef->GetY(Idx(i, 2, k))[1] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_j[Idx(i, 0 + iTemp, k)][2];
                }
                coef->GetZ(Idx(i, 1, k))[1] = NodeDifferece4thRight(diff_temp);
                coef->GetZ(Idx(i, 2, k))[1] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_j[Idx(i, nj - 2 - iTemp, k)][0];
                }
                coef->GetX(Idx(i, nj - 2, k))[1] = NodeDifferece4thLeft(diff_temp);
                coef->GetX(Idx(i, nj - 3, k))[1] = -NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_j[Idx(i, nj - 2 - iTemp, k)][1];
                }
                coef->GetY(Idx(i, nj - 2, k))[1] = NodeDifferece4thLeft(diff_temp);
                coef->GetY(Idx(i, nj - 3, k))[1] = -NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_j[Idx(i, nj - 2 - iTemp, k)][2];
                }
                coef->GetZ(Idx(i, nj - 2, k))[1] = NodeDifferece4thLeft(diff_temp);
                coef->GetZ(Idx(i, nj - 3, k))[1] = -NodeDifferece4th(diff_temp);
            }
        }
        // k=1,2;nk-3,nk-2
        for (int j = 1; j < nj - 1; ++j)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_k[Idx(i, j, 0 + iTemp)][0];
                }
                coef->GetX(Idx(i, j, 1))[2] = NodeDifferece4thRight(diff_temp);
                coef->GetX(Idx(i, j, 2))[2] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_k[Idx(i, j, 0 + iTemp)][1];
                }
                coef->GetY(Idx(i, j, 1))[2] = NodeDifferece4thRight(diff_temp);
                coef->GetY(Idx(i, j, 2))[2] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_k[Idx(i, j, 0 + iTemp)][2];
                }
                coef->GetZ(Idx(i, j, 1))[2] = NodeDifferece4thRight(diff_temp);
                coef->GetZ(Idx(i, j, 2))[2] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_k[Idx(i, j, nk - 2 - iTemp)][0];
                }
                coef->GetX(Idx(i, j, nk - 2))[2] = NodeDifferece4thLeft(diff_temp);
                coef->GetX(Idx(i, j, nk - 3))[2] = -NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_k[Idx(i, j, nk - 2 - iTemp)][1];
                }
                coef->GetY(Idx(i, j, nk - 2))[2] = NodeDifferece4thLeft(diff_temp);
                coef->GetY(Idx(i, j, nk - 3))[2] = -NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_k[Idx(i, j, nk - 2 - iTemp)][2];
                }
                coef->GetZ(Idx(i, j, nk - 2))[2] = NodeDifferece4thLeft(diff_temp);
                coef->GetZ(Idx(i, j, nk - 3))[2] = -NodeDifferece4th(diff_temp);
                if (grid->GetDim() == 2)
                {
                    coef->GetZ(Idx(i, j, 1))[0] = coef->GetZ(Idx(i, j, 1))[1] = 0.0;
                    coef->GetZ(Idx(i, j, 2))[0] = coef->GetZ(Idx(i, j, 2))[1] = 0.0;
                    coef->GetZ(Idx(i, j, nk - 2))[0] = coef->GetZ(Idx(i, j, nk - 2))[1] = 0.0;
                    coef->GetZ(Idx(i, j, nk - 3))[0] = coef->GetZ(Idx(i, j, nk - 3))[1] = 0.0;
                    coef->GetZ(Idx(i, j, 1))[2] = 1.0;
                    coef->GetZ(Idx(i, j, 2))[2] = 1.0;
                    coef->GetZ(Idx(i, j, nk - 2))[2] = 1.0;
                    coef->GetZ(Idx(i, j, nk - 3))[2] = 1.0;
                }
            }
        }
        // 第三步：计算半点处逆变换度量系数和坐标乘积
        //  求出xi方向的临时变量：逆变换度量系数乘以坐标
        struct TempXi
        {
            double z_zeta_y, x_zeta_z, y_zeta_x, z_eta_y, x_eta_z, y_eta_x;
        };
        std::vector<TempXi> temp_xi(ni * nj * nk);
        // i+1/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 3; i < ni - 4; ++i)
                {
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZ(Idx(i - 2 + iTemp, j, k))[2] * node->GetCoord(i - 2 + iTemp, j, k)[1];
                    }
                    temp_xi[Idx(i, j, k)].z_zeta_y = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetX(Idx(i - 2 + iTemp, j, k))[2] * node->GetCoord(i - 2 + iTemp, j, k)[2];
                    }
                    temp_xi[Idx(i, j, k)].x_zeta_z = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetY(Idx(i - 2 + iTemp, j, k))[2] * node->GetCoord(i - 2 + iTemp, j, k)[0];
                    }
                    temp_xi[Idx(i, j, k)].y_zeta_x = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZ(Idx(i - 2 + iTemp, j, k))[1] * node->GetCoord(i - 2 + iTemp, j, k)[1];
                    }
                    temp_xi[Idx(i, j, k)].z_eta_y = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetX(Idx(i - 2 + iTemp, j, k))[1] * node->GetCoord(i - 2 + iTemp, j, k)[2];
                    }
                    temp_xi[Idx(i, j, k)].x_eta_z = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetY(Idx(i - 2 + iTemp, j, k))[1] * node->GetCoord(i - 2 + iTemp, j, k)[0];
                    }
                    temp_xi[Idx(i, j, k)].y_eta_x = MidNodeInter6th(inter_temp);
                }
            }
        }
        // i=1/2,3/2,5/2;ni-3/2,ni-5/2,ni-7/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                // i=1/2，3/2，5/2
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(1 + iTemp, j, k))[2] * node->GetCoord(1 + iTemp, j, k)[1];
                }
                temp_xi[Idx(0, j, k)].z_zeta_y = MidNodeInter4thRight1(inter_temp);
                temp_xi[Idx(1, j, k)].z_zeta_y = MidNodeInter4thRight2(inter_temp);
                temp_xi[Idx(2, j, k)].z_zeta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(1 + iTemp, j, k))[2] * node->GetCoord(1 + iTemp, j, k)[2];
                }
                temp_xi[Idx(0, j, k)].x_zeta_z = MidNodeInter4thRight1(inter_temp);
                temp_xi[Idx(1, j, k)].x_zeta_z = MidNodeInter4thRight2(inter_temp);
                temp_xi[Idx(2, j, k)].x_zeta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(1 + iTemp, j, k))[2] * node->GetCoord(1 + iTemp, j, k)[0];
                }
                temp_xi[Idx(0, j, k)].y_zeta_x = MidNodeInter4thRight1(inter_temp);
                temp_xi[Idx(1, j, k)].y_zeta_x = MidNodeInter4thRight2(inter_temp);
                temp_xi[Idx(2, j, k)].y_zeta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(1 + iTemp, j, k))[1] * node->GetCoord(1 + iTemp, j, k)[1];
                }
                temp_xi[Idx(0, j, k)].z_eta_y = MidNodeInter4thRight1(inter_temp);
                temp_xi[Idx(1, j, k)].z_eta_y = MidNodeInter4thRight2(inter_temp);
                temp_xi[Idx(2, j, k)].z_eta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(1 + iTemp, j, k))[1] * node->GetCoord(1 + iTemp, j, k)[2];
                }
                temp_xi[Idx(0, j, k)].x_eta_z = MidNodeInter4thRight1(inter_temp);
                temp_xi[Idx(1, j, k)].x_eta_z = MidNodeInter4thRight2(inter_temp);
                temp_xi[Idx(2, j, k)].x_eta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(1 + iTemp, j, k))[1] * node->GetCoord(1 + iTemp, j, k)[0];
                }
                temp_xi[Idx(0, j, k)].y_eta_x = MidNodeInter4thRight1(inter_temp);
                temp_xi[Idx(1, j, k)].y_eta_x = MidNodeInter4thRight2(inter_temp);
                temp_xi[Idx(2, j, k)].y_eta_x = MidNodeInter4th(inter_temp);
                // i=ni-3/2,ni-5/2,ni-7/2
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(ni - 5 + iTemp, j, k))[2] * node->GetCoord(ni - 5 + iTemp, j, k)[1];
                }
                temp_xi[Idx(ni - 2, j, k)].z_zeta_y = MidNodeInter4thLeft1(inter_temp);
                temp_xi[Idx(ni - 3, j, k)].z_zeta_y = MidNodeInter4thLeft2(inter_temp);
                temp_xi[Idx(ni - 4, j, k)].z_zeta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(ni - 5 + iTemp, j, k))[2] * node->GetCoord(ni - 5 + iTemp, j, k)[2];
                }
                temp_xi[Idx(ni - 2, j, k)].x_zeta_z = MidNodeInter4thLeft1(inter_temp);
                temp_xi[Idx(ni - 3, j, k)].x_zeta_z = MidNodeInter4thLeft2(inter_temp);
                temp_xi[Idx(ni - 4, j, k)].x_zeta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(ni - 5 + iTemp, j, k))[2] * node->GetCoord(ni - 5 + iTemp, j, k)[0];
                }
                temp_xi[Idx(ni - 2, j, k)].y_zeta_x = MidNodeInter4thLeft1(inter_temp);
                temp_xi[Idx(ni - 3, j, k)].y_zeta_x = MidNodeInter4thLeft2(inter_temp);
                temp_xi[Idx(ni - 4, j, k)].y_zeta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(ni - 5 + iTemp, j, k))[1] * node->GetCoord(ni - 5 + iTemp, j, k)[1];
                }
                temp_xi[Idx(ni - 2, j, k)].z_eta_y = MidNodeInter4thLeft1(inter_temp);
                temp_xi[Idx(ni - 3, j, k)].z_eta_y = MidNodeInter4thLeft2(inter_temp);
                temp_xi[Idx(ni - 4, j, k)].z_eta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(ni - 5 + iTemp, j, k))[1] * node->GetCoord(ni - 5 + iTemp, j, k)[2];
                }
                temp_xi[Idx(ni - 2, j, k)].x_eta_z = MidNodeInter4thLeft1(inter_temp);
                temp_xi[Idx(ni - 3, j, k)].x_eta_z = MidNodeInter4thLeft2(inter_temp);
                temp_xi[Idx(ni - 4, j, k)].x_eta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(ni - 5 + iTemp, j, k))[1] * node->GetCoord(ni - 5 + iTemp, j, k)[0];
                }
                temp_xi[Idx(ni - 2, j, k)].y_eta_x = MidNodeInter4thLeft1(inter_temp);
                temp_xi[Idx(ni - 3, j, k)].y_eta_x = MidNodeInter4thLeft2(inter_temp);
                temp_xi[Idx(ni - 4, j, k)].y_eta_x = MidNodeInter4th(inter_temp);
            }
        }
        // 求出eta方向的临时变量：逆变换度量系数乘以坐标
        struct TempEta
        {
            double z_zeta_y, x_zeta_z, y_zeta_x, z_xi_y, x_xi_z, y_xi_x;
        };
        std::vector<TempEta> temp_eta(ni * nj * nk);
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 3; j < nj - 4; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZ(Idx(i, j - 2 + iTemp, k))[2] * node->GetCoord(i, j - 2 + iTemp, k)[1];
                    }
                    temp_eta[Idx(i, j, k)].z_zeta_y = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetX(Idx(i, j - 2 + iTemp, k))[2] * node->GetCoord(i, j - 2 + iTemp, k)[2];
                    }
                    temp_eta[Idx(i, j, k)].x_zeta_z = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetY(Idx(i, j - 2 + iTemp, k))[2] * node->GetCoord(i, j - 2 + iTemp, k)[0];
                    }
                    temp_eta[Idx(i, j, k)].y_zeta_x = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZ(Idx(i, j - 2 + iTemp, k))[0] * node->GetCoord(i, j - 2 + iTemp, k)[1];
                    }
                    temp_eta[Idx(i, j, k)].z_xi_y = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetX(Idx(i, j - 2 + iTemp, k))[0] * node->GetCoord(i, j - 2 + iTemp, k)[2];
                    }
                    temp_eta[Idx(i, j, k)].x_xi_z = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetY(Idx(i, j - 2 + iTemp, k))[0] * node->GetCoord(i, j - 2 + iTemp, k)[0];
                    }
                    temp_eta[Idx(i, j, k)].y_xi_x = MidNodeInter6th(inter_temp);
                }
            }
        }
        // j=1/2,3/2,5/2;nj-3/2,nj-5/2,nj-7/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                // j=1/2,3/2,5/2
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, 1 + iTemp, k))[2] * node->GetCoord(i, 1 + iTemp, k)[1];
                }
                temp_eta[Idx(i, 0, k)].z_zeta_y = MidNodeInter4thRight1(inter_temp);
                temp_eta[Idx(i, 1, k)].z_zeta_y = MidNodeInter4thRight2(inter_temp);
                temp_eta[Idx(i, 2, k)].z_zeta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, 1 + iTemp, k))[2] * node->GetCoord(i, 1 + iTemp, k)[2];
                }
                temp_eta[Idx(i, 0, k)].x_zeta_z = MidNodeInter4thRight1(inter_temp);
                temp_eta[Idx(i, 1, k)].x_zeta_z = MidNodeInter4thRight2(inter_temp);
                temp_eta[Idx(i, 2, k)].x_zeta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, 1 + iTemp, k))[2] * node->GetCoord(i, 1 + iTemp, k)[0];
                }
                temp_eta[Idx(i, 0, k)].y_zeta_x = MidNodeInter4thRight1(inter_temp);
                temp_eta[Idx(i, 1, k)].y_zeta_x = MidNodeInter4thRight2(inter_temp);
                temp_eta[Idx(i, 2, k)].y_zeta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, 1 + iTemp, k))[0] * node->GetCoord(i, 1 + iTemp, k)[1];
                }
                temp_eta[Idx(i, 0, k)].z_xi_y = MidNodeInter4thRight1(inter_temp);
                temp_eta[Idx(i, 1, k)].z_xi_y = MidNodeInter4thRight2(inter_temp);
                temp_eta[Idx(i, 2, k)].z_xi_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, 1 + iTemp, k))[0] * node->GetCoord(i, 1 + iTemp, k)[2];
                }
                temp_eta[Idx(i, 0, k)].x_xi_z = MidNodeInter4thRight1(inter_temp);
                temp_eta[Idx(i, 1, k)].x_xi_z = MidNodeInter4thRight2(inter_temp);
                temp_eta[Idx(i, 2, k)].x_xi_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, 1 + iTemp, k))[0] * node->GetCoord(i, 1 + iTemp, k)[0];
                }
                temp_eta[Idx(i, 0, k)].y_xi_x = MidNodeInter4thRight1(inter_temp);
                temp_eta[Idx(i, 1, k)].y_xi_x = MidNodeInter4thRight2(inter_temp);
                temp_eta[Idx(i, 2, k)].y_xi_x = MidNodeInter4th(inter_temp);
                // j=nj-3/2,nj-5/2,nj-7/2
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, nj - 5 + iTemp, k))[2] * node->GetCoord(i, nj - 5 + iTemp, k)[1];
                }
                temp_eta[Idx(i, nj - 2, k)].z_zeta_y = MidNodeInter4thLeft1(inter_temp);
                temp_eta[Idx(i, nj - 3, k)].z_zeta_y = MidNodeInter4thLeft2(inter_temp);
                temp_eta[Idx(i, nj - 4, k)].z_zeta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, nj - 5 + iTemp, k))[2] * node->GetCoord(i, nj - 5 + iTemp, k)[2];
                }
                temp_eta[Idx(i, nj - 2, k)].x_zeta_z = MidNodeInter4thLeft1(inter_temp);
                temp_eta[Idx(i, nj - 3, k)].x_zeta_z = MidNodeInter4thLeft2(inter_temp);
                temp_eta[Idx(i, nj - 4, k)].x_zeta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, nj - 5 + iTemp, k))[2] * node->GetCoord(i, nj - 5 + iTemp, k)[0];
                }
                temp_eta[Idx(i, nj - 2, k)].y_zeta_x = MidNodeInter4thLeft1(inter_temp);
                temp_eta[Idx(i, nj - 3, k)].y_zeta_x = MidNodeInter4thLeft2(inter_temp);
                temp_eta[Idx(i, nj - 4, k)].y_zeta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, nj - 5 + iTemp, k))[0] * node->GetCoord(i, nj - 5 + iTemp, k)[1];
                }
                temp_eta[Idx(i, nj - 2, k)].z_xi_y = MidNodeInter4thLeft1(inter_temp);
                temp_eta[Idx(i, nj - 3, k)].z_xi_y = MidNodeInter4thLeft2(inter_temp);
                temp_eta[Idx(i, nj - 4, k)].z_xi_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, nj - 5 + iTemp, k))[0] * node->GetCoord(i, nj - 5 + iTemp, k)[2];
                }
                temp_eta[Idx(i, nj - 2, k)].x_xi_z = MidNodeInter4thLeft1(inter_temp);
                temp_eta[Idx(i, nj - 3, k)].x_xi_z = MidNodeInter4thLeft2(inter_temp);
                temp_eta[Idx(i, nj - 4, k)].x_xi_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, nj - 5 + iTemp, k))[0] * node->GetCoord(i, nj - 5 + iTemp, k)[0];
                }
                temp_eta[Idx(i, nj - 2, k)].y_xi_x = MidNodeInter4thLeft1(inter_temp);
                temp_eta[Idx(i, nj - 3, k)].y_xi_x = MidNodeInter4thLeft2(inter_temp);
                temp_eta[Idx(i, nj - 4, k)].y_xi_x = MidNodeInter4th(inter_temp);
            }
        }
        // 求出zeta方向的临时变量：逆变换度量系数乘以坐标
        struct TempZeta
        {
            double z_eta_y, x_eta_z, y_eta_x, z_xi_y, x_xi_z, y_xi_x;
        };
        std::vector<TempZeta> temp_zeta(ni * nj * nk);
        for (int k = 3; k < nk - 4; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZ(Idx(i, j, k - 2 + iTemp))[1] * node->GetCoord(i, j, k - 2 + iTemp)[1];
                    }
                    temp_zeta[Idx(i, j, k)].z_eta_y = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetX(Idx(i, j, k - 2 + iTemp))[1] * node->GetCoord(i, j, k - 2 + iTemp)[2];
                    }
                    temp_zeta[Idx(i, j, k)].x_eta_z = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetY(Idx(i, j, k - 2 + iTemp))[1] * node->GetCoord(i, j, k - 2 + iTemp)[0];
                    }
                    temp_zeta[Idx(i, j, k)].y_eta_x = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZ(Idx(i, j, k - 2 + iTemp))[0] * node->GetCoord(i, j, k - 2 + iTemp)[1];
                    }
                    temp_zeta[Idx(i, j, k)].z_xi_y = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetX(Idx(i, j, k - 2 + iTemp))[0] * node->GetCoord(i, j, k - 2 + iTemp)[2];
                    }
                    temp_zeta[Idx(i, j, k)].x_xi_z = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetY(Idx(i, j, k - 2 + iTemp))[0] * node->GetCoord(i, j, k - 2 + iTemp)[0];
                    }
                    temp_zeta[Idx(i, j, k)].y_xi_x = MidNodeInter6th(inter_temp);
                }
            }
        }
        // k=1/2,3/2,5/2;nk-3/2,nk-5/2,nk-7/2
        for (int j = 1; j < nj - 1; ++j)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                // k=1/2,3/2,5/2
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, j, 1 + iTemp))[1] * node->GetCoord(i, j, 1 + iTemp)[1];
                }
                temp_zeta[Idx(i, j, 0)].z_eta_y = MidNodeInter4thRight1(inter_temp);
                temp_zeta[Idx(i, j, 1)].z_eta_y = MidNodeInter4thRight2(inter_temp);
                temp_zeta[Idx(i, j, 2)].z_eta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, j, 1 + iTemp))[1] * node->GetCoord(i, j, 1 + iTemp)[2];
                }
                temp_zeta[Idx(i, j, 0)].x_eta_z = MidNodeInter4thRight1(inter_temp);
                temp_zeta[Idx(i, j, 1)].x_eta_z = MidNodeInter4thRight2(inter_temp);
                temp_zeta[Idx(i, j, 2)].x_eta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, j, 1 + iTemp))[1] * node->GetCoord(i, j, 1 + iTemp)[0];
                }
                temp_zeta[Idx(i, j, 0)].y_eta_x = MidNodeInter4thRight1(inter_temp);
                temp_zeta[Idx(i, j, 1)].y_eta_x = MidNodeInter4thRight2(inter_temp);
                temp_zeta[Idx(i, j, 2)].y_eta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, j, 1 + iTemp))[0] * node->GetCoord(i, j, 1 + iTemp)[1];
                }
                temp_zeta[Idx(i, j, 0)].z_xi_y = MidNodeInter4thRight1(inter_temp);
                temp_zeta[Idx(i, j, 1)].z_xi_y = MidNodeInter4thRight2(inter_temp);
                temp_zeta[Idx(i, j, 2)].z_xi_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, j, 1 + iTemp))[0] * node->GetCoord(i, j, 1 + iTemp)[2];
                }
                temp_zeta[Idx(i, j, 0)].x_xi_z = MidNodeInter4thRight1(inter_temp);
                temp_zeta[Idx(i, j, 1)].x_xi_z = MidNodeInter4thRight2(inter_temp);
                temp_zeta[Idx(i, j, 2)].x_xi_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, j, 1 + iTemp))[0] * node->GetCoord(i, j, 1 + iTemp)[0];
                }
                temp_zeta[Idx(i, j, 0)].y_xi_x = MidNodeInter4thRight1(inter_temp);
                temp_zeta[Idx(i, j, 1)].y_xi_x = MidNodeInter4thRight2(inter_temp);
                temp_zeta[Idx(i, j, 2)].y_xi_x = MidNodeInter4th(inter_temp);
                // k=nk-3/2,nk-5/2,nk-7/2
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, j, nk - 5 + iTemp))[1] * node->GetCoord(i, j, nk - 5 + iTemp)[1];
                }
                temp_zeta[Idx(i, j, nk - 2)].z_eta_y = MidNodeInter4thLeft1(inter_temp);
                temp_zeta[Idx(i, j, nk - 3)].z_eta_y = MidNodeInter4thLeft2(inter_temp);
                temp_zeta[Idx(i, j, nk - 4)].z_eta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, j, nk - 5 + iTemp))[1] * node->GetCoord(i, j, nk - 5 + iTemp)[2];
                }
                temp_zeta[Idx(i, j, nk - 2)].x_eta_z = MidNodeInter4thLeft1(inter_temp);
                temp_zeta[Idx(i, j, nk - 3)].x_eta_z = MidNodeInter4thLeft2(inter_temp);
                temp_zeta[Idx(i, j, nk - 4)].x_eta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, j, nk - 5 + iTemp))[1] * node->GetCoord(i, j, nk - 5 + iTemp)[0];
                }
                temp_zeta[Idx(i, j, nk - 2)].y_eta_x = MidNodeInter4thLeft1(inter_temp);
                temp_zeta[Idx(i, j, nk - 3)].y_eta_x = MidNodeInter4thLeft2(inter_temp);
                temp_zeta[Idx(i, j, nk - 4)].y_eta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, j, nk - 5 + iTemp))[0] * node->GetCoord(i, j, nk - 5 + iTemp)[1];
                }
                temp_zeta[Idx(i, j, nk - 2)].z_xi_y = MidNodeInter4thLeft1(inter_temp);
                temp_zeta[Idx(i, j, nk - 3)].z_xi_y = MidNodeInter4thLeft2(inter_temp);
                temp_zeta[Idx(i, j, nk - 4)].z_xi_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, j, nk - 5 + iTemp))[0] * node->GetCoord(i, j, nk - 5 + iTemp)[2];
                }
                temp_zeta[Idx(i, j, nk - 2)].x_xi_z = MidNodeInter4thLeft1(inter_temp);
                temp_zeta[Idx(i, j, nk - 3)].x_xi_z = MidNodeInter4thLeft2(inter_temp);
                temp_zeta[Idx(i, j, nk - 4)].x_xi_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, j, nk - 5 + iTemp))[0] * node->GetCoord(i, j, nk - 5 + iTemp)[0];
                }
                temp_zeta[Idx(i, j, nk - 2)].y_xi_x = MidNodeInter4thLeft1(inter_temp);
                temp_zeta[Idx(i, j, nk - 3)].y_xi_x = MidNodeInter4thLeft2(inter_temp);
                temp_zeta[Idx(i, j, nk - 4)].y_xi_x = MidNodeInter4th(inter_temp);
            }
        }

        // 第五步：根据半点坐标和半点逆变换度量系数使用守恒形式计算整点度量系数（CMM1)
        // i direction
        double temp[6][6];
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 3; i < ni - 3; ++i)
                {
                    auto coef_xi = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, k));
                    auto coef_eta = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, k));
                    auto coef_zeta = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, k));
                    // i direction
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        int idx = Idx(i - 3 + iTemp, j, k);
                        temp[0][iTemp] = -temp_xi[idx].z_zeta_y; // -z_zeta*y
                        temp[1][iTemp] = -temp_xi[idx].x_zeta_z; // -x_zeta*z
                        temp[2][iTemp] = -temp_xi[idx].y_zeta_x; // -y_zeta*x
                        temp[3][iTemp] = temp_xi[idx].z_eta_y;   // z_eta*y
                        temp[4][iTemp] = temp_xi[idx].x_eta_z;   // x_eta*z
                        temp[5][iTemp] = temp_xi[idx].y_eta_x;   // y_eta*x
                    }
                    coef_eta[0] += NodeDifferece6th(temp[0]);
                    coef_eta[1] += NodeDifferece6th(temp[1]);
                    coef_eta[2] += NodeDifferece6th(temp[2]);
                    coef_zeta[0] += NodeDifferece6th(temp[3]);
                    coef_zeta[1] += NodeDifferece6th(temp[4]);
                    coef_zeta[2] += NodeDifferece6th(temp[5]);
                }
            }
        }
        // j direction
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 3; j < nj - 3; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    auto coef_xi = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, k));
                    auto coef_eta = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, k));
                    auto coef_zeta = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        int idx = Idx(i, j - 3 + iTemp, k);
                        temp[0][iTemp] = temp_eta[idx].z_zeta_y; // z_zeta*y
                        temp[1][iTemp] = temp_eta[idx].x_zeta_z; // x_zeta*z
                        temp[2][iTemp] = temp_eta[idx].y_zeta_x; // y_zeta*x
                        temp[3][iTemp] = -temp_eta[idx].z_xi_y;  // -z_xi*y
                        temp[4][iTemp] = -temp_eta[idx].x_xi_z;  // -x_xi*z
                        temp[5][iTemp] = -temp_eta[idx].y_xi_x;  // -y_xi*x
                    }
                    coef_xi[0] += NodeDifferece6th(temp[0]);
                    coef_xi[1] += NodeDifferece6th(temp[1]);
                    coef_xi[2] += NodeDifferece6th(temp[2]);
                    coef_zeta[0] += NodeDifferece6th(temp[3]);
                    coef_zeta[1] += NodeDifferece6th(temp[4]);
                    coef_zeta[2] += NodeDifferece6th(temp[5]);
                }
            }
        }
        // k direction
        for (int k = 3; k < nk - 3; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    auto coef_xi = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, k));
                    auto coef_eta = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, k));
                    auto coef_zeta = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        int idx = Idx(i, j, k - 3 + iTemp);
                        temp[0][iTemp] = -temp_zeta[idx].z_eta_y; // -z_eta*y
                        temp[1][iTemp] = -temp_zeta[idx].x_eta_z; // -x_eta*z
                        temp[2][iTemp] = -temp_zeta[idx].y_eta_x; // -y_eta*x
                        temp[3][iTemp] = temp_zeta[idx].z_xi_y;   // z_xi*y
                        temp[4][iTemp] = temp_zeta[idx].x_xi_z;   // x_xi*z
                        temp[5][iTemp] = temp_zeta[idx].y_xi_x;   // y_xi*x
                    }
                    coef_xi[0] += NodeDifferece6th(temp[0]);
                    coef_xi[1] += NodeDifferece6th(temp[1]);
                    coef_xi[2] += NodeDifferece6th(temp[2]);
                    coef_eta[0] += NodeDifferece6th(temp[3]);
                    coef_eta[1] += NodeDifferece6th(temp[4]);
                    coef_eta[2] += NodeDifferece6th(temp[5]);
                }
            }
        }
        // 第六步：计算边界点度量系数
        // i=1,2;ni-3,ni-2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                auto coef_xi1 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(1, j, k));
                auto coef_eta1 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(1, j, k));
                auto coef_zeta1 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(1, j, k));
                auto coef_xi2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(2, j, k));
                auto coef_eta2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(2, j, k));
                auto coef_zeta2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(2, j, k));
                for (int iTemp = 0; iTemp < 5; iTemp++)
                {
                    int idx = Idx(0 + iTemp, j, k);
                    temp[0][iTemp] = -temp_xi[idx].z_zeta_y; // -z_zeta*y
                    temp[1][iTemp] = -temp_xi[idx].x_zeta_z; // -x_zeta*z
                    temp[2][iTemp] = -temp_xi[idx].y_zeta_x; // -y_zeta*x
                    temp[3][iTemp] = temp_xi[idx].z_eta_y;   // z_eta*y
                    temp[4][iTemp] = temp_xi[idx].x_eta_z;   // x_eta*z
                    temp[5][iTemp] = temp_xi[idx].y_eta_x;   // y_eta*x
                }
                coef_eta1[0] += NodeDifferece4thRight(temp[0]);
                coef_eta1[1] += NodeDifferece4thRight(temp[1]);
                coef_eta1[2] += NodeDifferece4thRight(temp[2]);
                coef_zeta1[0] += NodeDifferece4thRight(temp[3]);
                coef_zeta1[1] += NodeDifferece4thRight(temp[4]);
                coef_zeta1[2] += NodeDifferece4thRight(temp[5]);
                coef_eta2[0] += NodeDifferece4th(temp[0]);
                coef_eta2[1] += NodeDifferece4th(temp[1]);
                coef_eta2[2] += NodeDifferece4th(temp[2]);
                coef_zeta2[0] += NodeDifferece4th(temp[3]);
                coef_zeta2[1] += NodeDifferece4th(temp[4]);
                coef_zeta2[2] = NodeDifferece4th(temp[5]);
                auto coef_xi_n2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(ni - 2, j, k));
                auto coef_eta_n2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(ni - 2, j, k));
                auto coef_zeta_n2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(ni - 2, j, k));
                auto coef_xi_n3 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(ni - 3, j, k));
                auto coef_eta_n3 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(ni - 3, j, k));
                auto coef_zeta_n3 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(ni - 3, j, k));
                for (int iTemp = 0; iTemp < 5; iTemp++)
                {
                    int idx = Idx(ni - 2 - iTemp, j, k);
                    temp[0][iTemp] = -temp_xi[idx].z_zeta_y; // -z_zeta*y
                    temp[1][iTemp] = -temp_xi[idx].x_zeta_z; // -x_zeta*z
                    temp[2][iTemp] = -temp_xi[idx].y_zeta_x; // -y_zeta*x
                    temp[3][iTemp] = temp_xi[idx].z_eta_y;   // z_eta*y
                    temp[4][iTemp] = temp_xi[idx].x_eta_z;   // x_eta*z
                    temp[5][iTemp] = temp_xi[idx].y_eta_x;   // y_eta*x
                }
                coef_eta_n2[0] += NodeDifferece4thLeft(temp[0]);
                coef_eta_n2[1] += NodeDifferece4thLeft(temp[1]);
                coef_eta_n2[2] += NodeDifferece4thLeft(temp[2]);
                coef_zeta_n2[0] += NodeDifferece4thLeft(temp[3]);
                coef_zeta_n2[1] += NodeDifferece4thLeft(temp[4]);
                coef_zeta_n2[2] += NodeDifferece4thLeft(temp[5]);
                coef_eta_n3[0] += -NodeDifferece4th(temp[0]);
                coef_eta_n3[1] += -NodeDifferece4th(temp[1]);
                coef_eta_n3[2] += -NodeDifferece4th(temp[2]);
                coef_zeta_n3[0] += -NodeDifferece4th(temp[3]);
                coef_zeta_n3[1] += -NodeDifferece4th(temp[4]);
                coef_zeta_n3[2] += -NodeDifferece4th(temp[5]);
            }
        }
        // j=1,2;nj-3,nj-2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                auto coef_xi1 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, 1, k));
                auto coef_eta1 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, 1, k));
                auto coef_zeta1 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, 1, k));
                auto coef_xi2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, 2, k));
                auto coef_eta2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, 2, k));
                auto coef_zeta2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, 2, k));
                for (int iTemp = 0; iTemp < 5; iTemp++)
                {
                    int idx = Idx(i, 0 + iTemp, k);
                    temp[0][iTemp] = temp_eta[idx].z_zeta_y; // z_zeta*y
                    temp[1][iTemp] = temp_eta[idx].x_zeta_z; // x_zeta*z
                    temp[2][iTemp] = temp_eta[idx].y_zeta_x; // y_zeta*x
                    temp[3][iTemp] = -temp_eta[idx].z_xi_y;  // -z_xi*y
                    temp[4][iTemp] = -temp_eta[idx].x_xi_z;  // -x_xi*z
                    temp[5][iTemp] = -temp_eta[idx].y_xi_x;  // -y_xi*x
                }
                coef_xi1[0] += NodeDifferece4thRight(temp[0]);
                coef_xi1[1] += NodeDifferece4thRight(temp[1]);
                coef_xi1[2] += NodeDifferece4thRight(temp[2]);
                coef_zeta1[0] += NodeDifferece4thRight(temp[3]);
                coef_zeta1[1] += NodeDifferece4thRight(temp[4]);
                coef_zeta1[2] += NodeDifferece4thRight(temp[5]);
                coef_xi2[0] += NodeDifferece4th(temp[0]);
                coef_xi2[1] += NodeDifferece4th(temp[1]);
                coef_xi2[2] += NodeDifferece4th(temp[2]);
                coef_zeta2[0] += NodeDifferece4th(temp[3]);
                coef_zeta2[1] += NodeDifferece4th(temp[4]);
                coef_zeta2[2] += NodeDifferece4th(temp[5]);
                auto coef_xi_n2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, nj - 2, k));
                auto coef_eta_n2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, nj - 2, k));
                auto coef_zeta_n2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, nj - 2, k));
                auto coef_xi_n3 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, nj - 3, k));
                auto coef_eta_n3 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, nj - 3, k));
                auto coef_zeta_n3 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, nj - 3, k));
                for (int iTemp = 0; iTemp < 5; iTemp++)
                {
                    int idx = Idx(i, nj - 2 - iTemp, k);
                    temp[0][iTemp] = temp_eta[idx].z_zeta_y; // z_zeta*y
                    temp[1][iTemp] = temp_eta[idx].x_zeta_z; // x_zeta*z
                    temp[2][iTemp] = temp_eta[idx].y_zeta_x; // y_zeta*x
                    temp[3][iTemp] = -temp_eta[idx].z_xi_y;  // -z_xi*y
                    temp[4][iTemp] = -temp_eta[idx].x_xi_z;  // -x_xi*z
                    temp[5][iTemp] = -temp_eta[idx].y_xi_x;  // -y_xi*x
                }
                coef_xi_n2[0] += NodeDifferece4thLeft(temp[0]);
                coef_xi_n2[1] += NodeDifferece4thLeft(temp[1]);
                coef_xi_n2[2] += NodeDifferece4thLeft(temp[2]);
                coef_zeta_n2[0] += NodeDifferece4thLeft(temp[3]);
                coef_zeta_n2[1] += NodeDifferece4thLeft(temp[4]);
                coef_zeta_n2[2] += NodeDifferece4thLeft(temp[5]);
                coef_xi_n3[0] += -NodeDifferece4th(temp[0]);
                coef_xi_n3[1] += -NodeDifferece4th(temp[1]);
                coef_xi_n3[2] += -NodeDifferece4th(temp[2]);
                coef_zeta_n3[0] += -NodeDifferece4th(temp[3]);
                coef_zeta_n3[1] += -NodeDifferece4th(temp[4]);
                coef_zeta_n3[2] += -NodeDifferece4th(temp[5]);
            }
        }
        // k=1,2;nk-3,nk-2
        for (int j = 1; j < nj - 1; ++j)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                auto coef_xi1 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, 1));
                auto coef_eta1 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, 1));
                auto coef_zeta1 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, 1));
                auto coef_xi2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, 2));
                auto coef_eta2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, 2));
                auto coef_zeta2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, 2));
                for (int iTemp = 0; iTemp < 5; iTemp++)
                {
                    int idx = Idx(i, j, 0 + iTemp);
                    temp[0][iTemp] = -temp_zeta[idx].z_eta_y; // -z_eta*y
                    temp[1][iTemp] = -temp_zeta[idx].x_eta_z; // -x_eta*z
                    temp[2][iTemp] = -temp_zeta[idx].y_eta_x; // -y_eta*x
                    temp[3][iTemp] = temp_zeta[idx].z_xi_y;   // z_xi*y
                    temp[4][iTemp] = temp_zeta[idx].x_xi_z;   // x_xi*z
                    temp[5][iTemp] = temp_zeta[idx].y_xi_x;   // y_xi*x
                }
                coef_xi1[0] += NodeDifferece4thRight(temp[0]);
                coef_xi1[1] += NodeDifferece4thRight(temp[1]);
                coef_xi1[2] += NodeDifferece4thRight(temp[2]);
                coef_eta1[0] += NodeDifferece4thRight(temp[3]);
                coef_eta1[1] += NodeDifferece4thRight(temp[4]);
                coef_eta1[2] += NodeDifferece4thRight(temp[5]);
                coef_xi2[0] += NodeDifferece4th(temp[0]);
                coef_xi2[1] += NodeDifferece4th(temp[1]);
                coef_xi2[2] += NodeDifferece4th(temp[2]);
                coef_eta2[0] += NodeDifferece4th(temp[3]);
                coef_eta2[1] += NodeDifferece4th(temp[4]);
                coef_eta2[2] += NodeDifferece4th(temp[5]);
                auto coef_xi_n2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, nk - 2));
                auto coef_eta_n2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, nk - 2));
                auto coef_zeta_n2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, nk - 2));
                auto coef_xi_n3 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, nk - 3));
                auto coef_eta_n3 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, nk - 3));
                auto coef_zeta_n3 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, nk - 3));
                for (int iTemp = 0; iTemp < 5; iTemp++)
                {
                    int idx = Idx(i, j, nk - 2 - iTemp);
                    temp[0][iTemp] = -temp_zeta[idx].z_eta_y; // -z_eta*y
                    temp[1][iTemp] = -temp_zeta[idx].x_eta_z; // -x_eta*z
                    temp[2][iTemp] = -temp_zeta[idx].y_eta_x; // -y_eta*x
                    temp[3][iTemp] = temp_zeta[idx].z_xi_y;   // z_xi*y
                    temp[4][iTemp] = temp_zeta[idx].x_xi_z;   // x_xi*z
                    temp[5][iTemp] = temp_zeta[idx].y_xi_x;   // y_xi*x
                }
                coef_xi_n2[0] += NodeDifferece4thLeft(temp[0]);
                coef_xi_n2[1] += NodeDifferece4thLeft(temp[1]);
                coef_xi_n2[2] += NodeDifferece4thLeft(temp[2]);
                coef_eta_n2[0] += NodeDifferece4thLeft(temp[3]);
                coef_eta_n2[1] += NodeDifferece4thLeft(temp[4]);
                coef_eta_n2[2] += NodeDifferece4thLeft(temp[5]);
                coef_xi_n3[0] += -NodeDifferece4th(temp[0]);
                coef_xi_n3[1] += -NodeDifferece4th(temp[1]);
                coef_xi_n3[2] += -NodeDifferece4th(temp[2]);
                coef_eta_n3[0] += -NodeDifferece4th(temp[3]);
                coef_eta_n3[1] += -NodeDifferece4th(temp[4]);
                coef_eta_n3[2] += -NodeDifferece4th(temp[5]);
            }
        }
    }
    void NSSolverStruct::CalcMetricsS3()
    {
        auto para = GetPara();
        auto flux_diff_scheme = para->GetDifferenceScheme();
        if (flux_diff_scheme == DifferenceScheme::SecondOrder)
        {
            CalcMetricsS3_2nd();
        }
        else if (flux_diff_scheme == DifferenceScheme::SixthOrder)
        {
            CalcMetricsS3_6th();
        }
        else
        {
            Log::warn("FluxDifferenceScheme is not defined, use SecondOrder as default");
            CalcMetricsS3_2nd();
        }
    }
    void NSSolverStruct::CalcMetricsS3_2nd()
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
        int idx_temp[3];
        // i+1/2,j+1/2,k+1/2处的坐标
        std::vector<std::vector<double>>
            coord_i(ni * nj * nk, std::vector<double>(3)),
            coord_j(ni * nj * nk, std::vector<double>(3)), coord_k(ni * nj * nk, std::vector<double>(3));
        // 计算之前，先把度量系数赋值为0
        for (int i = 0; i < ni; ++i)
        {
            for (int j = 0; j < nj; ++j)
            {
                for (int k = 0; k < nk; ++k)
                {
                    for (int iDim = 0; iDim < 4; ++iDim)
                    {
                        coef->GetX(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetY(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetZ(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetXi(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetEta(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetZeta(Idx(i, j, k))[iDim] = 0.0;
                        for (int jDim = 0; jDim < 3; ++jDim)
                        {
                            coef_mid[jDim]->GetX(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetY(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetZ(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetXi(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetEta(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetZeta(Idx(i, j, k))[iDim] = 0.0;
                        }
                    }
                }
            }
        }
        // 第一步：计算半点坐标
        for (int k = 0; k < nk - 1; ++k)
        {
            for (int j = 0; j < nj - 1; ++j)
            {
                for (int i = 0; i < ni - 1; ++i)
                {
                    for (int iDim = 0; iDim < 3; ++iDim)
                    {
                        coord_i[Idx(i, j, k)][iDim] =
                            0.5 * (node->GetCoord(i, j, k)[iDim] + node->GetCoord(i + 1, j, k)[iDim]);
                        coord_j[Idx(i, j, k)][iDim] =
                            0.5 * (node->GetCoord(i, j, k)[iDim] + node->GetCoord(i, j + 1, k)[iDim]);
                        coord_k[Idx(i, j, k)][iDim] =
                            0.5 * (node->GetCoord(i, j, k)[iDim] + node->GetCoord(i, j, k + 1)[iDim]);
                    }
                }
            }
        }
        // 第二步：根据半点坐标计算整点逆变换度量系数
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    int idx = idx_proxy->GetIdx(i, j, k);
                    auto coef_x = coef->GetX(idx);
                    auto coef_y = coef->GetY(idx);
                    auto coef_z = coef->GetZ(idx);
                    coef_x[0] = coord_i[Idx(i, j, k)][0] - coord_i[Idx(i - 1, j, k)][0];
                    coef_x[1] = coord_j[Idx(i, j, k)][0] - coord_j[Idx(i, j - 1, k)][0];
                    coef_x[2] = coord_k[Idx(i, j, k)][0] - coord_k[Idx(i, j, k - 1)][0];
                    coef_y[0] = coord_i[Idx(i, j, k)][1] - coord_i[Idx(i - 1, j, k)][1];
                    coef_y[1] = coord_j[Idx(i, j, k)][1] - coord_j[Idx(i, j - 1, k)][1];
                    coef_y[2] = coord_k[Idx(i, j, k)][1] - coord_k[Idx(i, j, k - 1)][1];
                    coef_z[0] = coord_i[Idx(i, j, k)][2] - coord_i[Idx(i - 1, j, k)][2];
                    coef_z[1] = coord_j[Idx(i, j, k)][2] - coord_j[Idx(i, j - 1, k)][2];
                    coef_z[2] = coord_k[Idx(i, j, k)][2] - coord_k[Idx(i, j, k - 1)][2];
                    coef_x[3] = 0.0;
                    coef_y[3] = 0.0;
                    coef_z[3] = 0.0;
                    if (grid->GetDim() == 2)
                    {
                        coef_z[0] = 0.0;
                        coef_z[1] = 0.0;
                        coef_z[2] = 1.0;
                    }
                }
            }
        }
        // 第三步：计算半点处逆变换度量系数和坐标乘积
        //  求出xi方向的临时变量：逆变换度量系数乘以坐标
        struct TempXi
        {
            double y_zeta_z, z_zeta_x, x_zeta_y, y_eta_z, z_eta_x, x_eta_y, z_zeta_y, x_zeta_z, y_zeta_x, z_eta_y, x_eta_z, y_eta_x;
        };
        std::vector<TempXi> temp_xi(ni * nj * nk);
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    int idx = Idx(i, j, k);
                    auto &temp = temp_xi[idx];
                    temp.y_zeta_z = 0.5 * (coef->GetY(Idx(i + 1, j, k))[2] * node->GetCoord(i + 1, j, k)[2] +
                                           coef->GetY(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[2]);
                    temp.z_zeta_x = 0.5 * (coef->GetZ(Idx(i + 1, j, k))[2] * node->GetCoord(i + 1, j, k)[0] +
                                           coef->GetZ(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[0]);
                    temp.x_zeta_y = 0.5 * (coef->GetX(Idx(i + 1, j, k))[2] * node->GetCoord(i + 1, j, k)[1] +
                                           coef->GetX(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[1]);
                    temp.y_eta_z = 0.5 * (coef->GetY(Idx(i + 1, j, k))[1] * node->GetCoord(i + 1, j, k)[2] +
                                          coef->GetY(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[2]);
                    temp.z_eta_x = 0.5 * (coef->GetZ(Idx(i + 1, j, k))[1] * node->GetCoord(i + 1, j, k)[0] +
                                          coef->GetZ(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[0]);
                    temp.x_eta_y = 0.5 * (coef->GetX(Idx(i + 1, j, k))[1] * node->GetCoord(i + 1, j, k)[1] +
                                          coef->GetX(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[1]);
                    temp.z_zeta_y = 0.5 * (coef->GetZ(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[1] +
                                           coef->GetZ(Idx(i + 1, j, k))[2] * node->GetCoord(i + 1, j, k)[1]);
                    temp.x_zeta_z = 0.5 * (coef->GetX(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[2] +
                                           coef->GetX(Idx(i + 1, j, k))[2] * node->GetCoord(i + 1, j, k)[2]);
                    temp.y_zeta_x = 0.5 * (coef->GetY(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[0] +
                                           coef->GetY(Idx(i + 1, j, k))[2] * node->GetCoord(i + 1, j, k)[0]);
                    temp.z_eta_y = 0.5 * (coef->GetZ(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[1] +
                                          coef->GetZ(Idx(i + 1, j, k))[1] * node->GetCoord(i + 1, j, k)[1]);
                    temp.x_eta_z = 0.5 * (coef->GetX(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[2] +
                                          coef->GetX(Idx(i + 1, j, k))[1] * node->GetCoord(i + 1, j, k)[2]);
                    temp.y_eta_x = 0.5 * (coef->GetY(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[0] +
                                          coef->GetY(Idx(i + 1, j, k))[1] * node->GetCoord(i + 1, j, k)[0]);
                }
            }
        }
        // i=1/2；(ni-1)-1/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                int idx = Idx(0, j, k);
                temp_xi[idx].y_zeta_z = 0.5 * (3.0 * coef->GetY(Idx(1, j, k))[2] * node->GetCoord(1, j, k)[2] -
                                               1.0 * coef->GetY(Idx(2, j, k))[2] * node->GetCoord(2, j, k)[2]);
                temp_xi[idx].z_zeta_x = 0.5 * (3.0 * coef->GetZ(Idx(1, j, k))[2] * node->GetCoord(1, j, k)[0] -
                                               1.0 * coef->GetZ(Idx(2, j, k))[2] * node->GetCoord(2, j, k)[0]);
                temp_xi[idx].x_zeta_y = 0.5 * (3.0 * coef->GetX(Idx(1, j, k))[2] * node->GetCoord(1, j, k)[1] -
                                               1.0 * coef->GetX(Idx(2, j, k))[2] * node->GetCoord(2, j, k)[1]);
                temp_xi[idx].y_eta_z = 0.5 * (3.0 * coef->GetY(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[2] -
                                              1.0 * coef->GetY(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[2]);
                temp_xi[idx].z_eta_x = 0.5 * (3.0 * coef->GetZ(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[0] -
                                              1.0 * coef->GetZ(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[0]);
                temp_xi[idx].x_eta_y = 0.5 * (3.0 * coef->GetX(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[1] -
                                              1.0 * coef->GetX(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[1]);
                temp_xi[Idx(0, j, k)].z_zeta_y = 0.5 * (3.0 * coef->GetZ(Idx(1, j, k))[2] * node->GetCoord(1, j, k)[1] - 1.0 * coef->GetZ(Idx(2, j, k))[2] * node->GetCoord(2, j, k)[1]);
                temp_xi[Idx(0, j, k)].x_zeta_z = 0.5 * (3.0 * coef->GetX(Idx(1, j, k))[2] * node->GetCoord(1, j, k)[2] - 1.0 * coef->GetX(Idx(2, j, k))[2] * node->GetCoord(2, j, k)[2]);
                temp_xi[Idx(0, j, k)].y_zeta_x = 0.5 * (3.0 * coef->GetY(Idx(1, j, k))[2] * node->GetCoord(1, j, k)[0] - 1.0 * coef->GetY(Idx(2, j, k))[2] * node->GetCoord(2, j, k)[0]);
                temp_xi[Idx(0, j, k)].z_eta_y = 0.5 * (3.0 * coef->GetZ(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[1] - 1.0 * coef->GetZ(Idx(2, j, k))[1] * node->GetCoord(2, j, k)[1]);
                temp_xi[Idx(0, j, k)].x_eta_z = 0.5 * (3.0 * coef->GetX(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[2] - 1.0 * coef->GetX(Idx(2, j, k))[1] * node->GetCoord(2, j, k)[2]);
                temp_xi[Idx(0, j, k)].y_eta_x = 0.5 * (3.0 * coef->GetY(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[0] - 1.0 * coef->GetY(Idx(2, j, k))[1] * node->GetCoord(2, j, k)[0]);

                idx = Idx(ni - 2, j, k);
                temp_xi[idx].y_zeta_z = 0.5 * (3.0 * coef->GetY(Idx(ni - 2, j, k))[2] * node->GetCoord(ni - 2, j, k)[2] -
                                               1.0 * coef->GetY(Idx(ni - 3, j, k))[2] * node->GetCoord(ni - 3, j, k)[2]);
                temp_xi[idx].z_zeta_x = 0.5 * (3.0 * coef->GetZ(Idx(ni - 2, j, k))[2] * node->GetCoord(ni - 2, j, k)[0] -
                                               1.0 * coef->GetZ(Idx(ni - 3, j, k))[2] * node->GetCoord(ni - 3, j, k)[0]);
                temp_xi[idx].x_zeta_y = 0.5 * (3.0 * coef->GetX(Idx(ni - 2, j, k))[2] * node->GetCoord(ni - 2, j, k)[1] -
                                               1.0 * coef->GetX(Idx(ni - 3, j, k))[2] * node->GetCoord(ni - 3, j, k)[1]);
                temp_xi[idx].y_eta_z = 0.5 * (3.0 * coef->GetY(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[2] -
                                              1.0 * coef->GetY(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[2]);
                temp_xi[idx].z_eta_x = 0.5 * (3.0 * coef->GetZ(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[0] -
                                              1.0 * coef->GetZ(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[0]);
                temp_xi[idx].x_eta_y = 0.5 * (3.0 * coef->GetX(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[1] -
                                              1.0 * coef->GetX(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[1]);
                temp_xi[Idx(ni - 2, j, k)].z_zeta_y = 0.5 * (3.0 * coef->GetZ(Idx(ni - 2, j, k))[2] * node->GetCoord(ni - 2, j, k)[1] - 1.0 * coef->GetZ(Idx(ni - 3, j, k))[2] * node->GetCoord(ni - 3, j, k)[1]);
                temp_xi[Idx(ni - 2, j, k)].x_zeta_z = 0.5 * (3.0 * coef->GetX(Idx(ni - 2, j, k))[2] * node->GetCoord(ni - 2, j, k)[2] - 1.0 * coef->GetX(Idx(ni - 3, j, k))[2] * node->GetCoord(ni - 3, j, k)[2]);
                temp_xi[Idx(ni - 2, j, k)].y_zeta_x = 0.5 * (3.0 * coef->GetY(Idx(ni - 2, j, k))[2] * node->GetCoord(ni - 2, j, k)[0] - 1.0 * coef->GetY(Idx(ni - 3, j, k))[2] * node->GetCoord(ni - 3, j, k)[0]);
                temp_xi[Idx(ni - 2, j, k)].z_eta_y = 0.5 * (3.0 * coef->GetZ(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[1] - 1.0 * coef->GetZ(Idx(ni - 3, j, k))[1] * node->GetCoord(ni - 3, j, k)[1]);
                temp_xi[Idx(ni - 2, j, k)].x_eta_z = 0.5 * (3.0 * coef->GetX(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[2] - 1.0 * coef->GetX(Idx(ni - 3, j, k))[1] * node->GetCoord(ni - 3, j, k)[2]);
                temp_xi[Idx(ni - 2, j, k)].y_eta_x = 0.5 * (3.0 * coef->GetY(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[0] - 1.0 * coef->GetY(Idx(ni - 3, j, k))[1] * node->GetCoord(ni - 3, j, k)[0]);
            }
        }
        // 求出eta方向的临时变量：逆变换度量系数乘以坐标
        struct TempEta
        {
            double y_zeta_z, z_zeta_x, x_zeta_y, y_xi_z, z_xi_x, x_xi_y, z_zeta_y, x_zeta_z, y_zeta_x, z_xi_y, x_xi_z, y_xi_x;
        };
        std::vector<TempEta> temp_eta(ni * nj * nk);
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    int idx = Idx(i, j, k);
                    auto &temp = temp_eta[idx];
                    temp.y_zeta_z = 0.5 * (coef->GetY(Idx(i, j + 1, k))[2] * node->GetCoord(i, j + 1, k)[2] +
                                           coef->GetY(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[2]);
                    temp.z_zeta_x = 0.5 * (coef->GetZ(Idx(i, j + 1, k))[2] * node->GetCoord(i, j + 1, k)[0] +
                                           coef->GetZ(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[0]);
                    temp.x_zeta_y = 0.5 * (coef->GetX(Idx(i, j + 1, k))[2] * node->GetCoord(i, j + 1, k)[1] +
                                           coef->GetX(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[1]);
                    temp.y_xi_z = 0.5 * (coef->GetY(Idx(i, j + 1, k))[0] * node->GetCoord(i, j + 1, k)[2] +
                                         coef->GetY(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[2]);
                    temp.z_xi_x = 0.5 * (coef->GetZ(Idx(i, j + 1, k))[0] * node->GetCoord(i, j + 1, k)[0] +
                                         coef->GetZ(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[0]);
                    temp.x_xi_y = 0.5 * (coef->GetX(Idx(i, j + 1, k))[0] * node->GetCoord(i, j + 1, k)[1] +
                                         coef->GetX(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[1]);
                    temp.z_zeta_y = 0.5 * (coef->GetZ(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[1] +
                                           coef->GetZ(Idx(i, j + 1, k))[2] * node->GetCoord(i, j + 1, k)[1]);
                    temp.x_zeta_z = 0.5 * (coef->GetX(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[2] +
                                           coef->GetX(Idx(i, j + 1, k))[2] * node->GetCoord(i, j + 1, k)[2]);
                    temp.y_zeta_x = 0.5 * (coef->GetY(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[0] +
                                           coef->GetY(Idx(i, j + 1, k))[2] * node->GetCoord(i, j + 1, k)[0]);
                    temp.z_xi_y = 0.5 * (coef->GetZ(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[1] +
                                         coef->GetZ(Idx(i, j + 1, k))[0] * node->GetCoord(i, j + 1, k)[1]);
                    temp.x_xi_z = 0.5 * (coef->GetX(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[2] +
                                         coef->GetX(Idx(i, j + 1, k))[0] * node->GetCoord(i, j + 1, k)[2]);
                    temp.y_xi_x = 0.5 * (coef->GetY(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[0] +
                                         coef->GetY(Idx(i, j + 1, k))[0] * node->GetCoord(i, j + 1, k)[0]);
                }
            }
        }
        // j=1/2；(nj-1)-1/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                int idx = Idx(i, 0, k);
                temp_eta[idx].y_zeta_z = 0.5 * (3.0 * coef->GetY(Idx(i, 1, k))[2] * node->GetCoord(i, 1, k)[2] -
                                                1.0 * coef->GetY(Idx(i, 2, k))[2] * node->GetCoord(i, 2, k)[2]);
                temp_eta[idx].z_zeta_x = 0.5 * (3.0 * coef->GetZ(Idx(i, 1, k))[2] * node->GetCoord(i, 1, k)[0] -
                                                1.0 * coef->GetZ(Idx(i, 2, k))[2] * node->GetCoord(i, 2, k)[0]);
                temp_eta[idx].x_zeta_y = 0.5 * (3.0 * coef->GetX(Idx(i, 1, k))[2] * node->GetCoord(i, 1, k)[1] -
                                                1.0 * coef->GetX(Idx(i, 2, k))[2] * node->GetCoord(i, 2, k)[1]);
                temp_eta[idx].y_xi_z = 0.5 * (3.0 * coef->GetY(Idx(i, 1, k))[0] * node->GetCoord(i, 1, k)[2] -
                                              1.0 * coef->GetY(Idx(i, 2, k))[0] * node->GetCoord(i, 2, k)[2]);
                temp_eta[idx].z_xi_x = 0.5 * (3.0 * coef->GetZ(Idx(i, 1, k))[0] * node->GetCoord(i, 1, k)[0] -
                                              1.0 * coef->GetZ(Idx(i, 2, k))[0] * node->GetCoord(i, 2, k)[0]);
                temp_eta[idx].x_xi_y = 0.5 * (3.0 * coef->GetX(Idx(i, 1, k))[0] * node->GetCoord(i, 1, k)[1] -
                                              1.0 * coef->GetX(Idx(i, 2, k))[0] * node->GetCoord(i, 2, k)[1]);
                temp_eta[Idx(i, 0, k)].z_zeta_y = 0.5 * (3.0 * coef->GetZ(Idx(i, 1, k))[2] * node->GetCoord(i, 1, k)[1] - 1.0 * coef->GetZ(Idx(i, 2, k))[2] * node->GetCoord(i, 2, k)[1]);
                temp_eta[Idx(i, 0, k)].x_zeta_z = 0.5 * (3.0 * coef->GetX(Idx(i, 1, k))[2] * node->GetCoord(i, 1, k)[2] - 1.0 * coef->GetX(Idx(i, 2, k))[2] * node->GetCoord(i, 2, k)[2]);
                temp_eta[Idx(i, 0, k)].y_zeta_x = 0.5 * (3.0 * coef->GetY(Idx(i, 1, k))[2] * node->GetCoord(i, 1, k)[0] - 1.0 * coef->GetY(Idx(i, 2, k))[2] * node->GetCoord(i, 2, k)[0]);
                temp_eta[Idx(i, 0, k)].z_xi_y = 0.5 * (3.0 * coef->GetZ(Idx(i, 1, k))[0] * node->GetCoord(i, 1, k)[1] - 1.0 * coef->GetZ(Idx(i, 2, k))[0] * node->GetCoord(i, 2, k)[1]);
                temp_eta[Idx(i, 0, k)].x_xi_z = 0.5 * (3.0 * coef->GetX(Idx(i, 1, k))[0] * node->GetCoord(i, 1, k)[2] - 1.0 * coef->GetX(Idx(i, 2, k))[0] * node->GetCoord(i, 2, k)[2]);
                temp_eta[Idx(i, 0, k)].y_xi_x = 0.5 * (3.0 * coef->GetY(Idx(i, 1, k))[0] * node->GetCoord(i, 1, k)[0] - 1.0 * coef->GetY(Idx(i, 2, k))[0] * node->GetCoord(i, 2, k)[0]);

                idx = Idx(i, nj - 2, k);
                temp_eta[idx].y_zeta_z = 0.5 * (3.0 * coef->GetY(Idx(i, nj - 2, k))[2] * node->GetCoord(i, nj - 2, k)[2] -
                                                1.0 * coef->GetY(Idx(i, nj - 3, k))[2] * node->GetCoord(i, nj - 3, k)[2]);
                temp_eta[idx].z_zeta_x = 0.5 * (3.0 * coef->GetZ(Idx(i, nj - 2, k))[2] * node->GetCoord(i, nj - 2, k)[0] -
                                                1.0 * coef->GetZ(Idx(i, nj - 3, k))[2] * node->GetCoord(i, nj - 3, k)[0]);
                temp_eta[idx].x_zeta_y = 0.5 * (3.0 * coef->GetX(Idx(i, nj - 2, k))[2] * node->GetCoord(i, nj - 2, k)[1] -
                                                1.0 * coef->GetX(Idx(i, nj - 3, k))[2] * node->GetCoord(i, nj - 3, k)[1]);
                temp_eta[idx].y_xi_z = 0.5 * (3.0 * coef->GetY(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[2] -
                                              1.0 * coef->GetY(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[2]);
                temp_eta[idx].z_xi_x = 0.5 * (3.0 * coef->GetZ(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[0] -
                                              1.0 * coef->GetZ(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[0]);
                temp_eta[idx].x_xi_y = 0.5 * (3.0 * coef->GetX(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[1] -
                                              1.0 * coef->GetX(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[1]);
                temp_eta[Idx(i, nj - 2, k)].z_zeta_y = 0.5 * (3.0 * coef->GetZ(Idx(i, nj - 2, k))[2] * node->GetCoord(i, nj - 2, k)[1] - 1.0 * coef->GetZ(Idx(i, nj - 3, k))[2] * node->GetCoord(i, nj - 3, k)[1]);
                temp_eta[Idx(i, nj - 2, k)].x_zeta_z = 0.5 * (3.0 * coef->GetX(Idx(i, nj - 2, k))[2] * node->GetCoord(i, nj - 2, k)[2] - 1.0 * coef->GetX(Idx(i, nj - 3, k))[2] * node->GetCoord(i, nj - 3, k)[2]);
                temp_eta[Idx(i, nj - 2, k)].y_zeta_x = 0.5 * (3.0 * coef->GetY(Idx(i, nj - 2, k))[2] * node->GetCoord(i, nj - 2, k)[0] - 1.0 * coef->GetY(Idx(i, nj - 3, k))[2] * node->GetCoord(i, nj - 3, k)[0]);
                temp_eta[Idx(i, nj - 2, k)].z_xi_y = 0.5 * (3.0 * coef->GetZ(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[1] - 1.0 * coef->GetZ(Idx(i, nj - 3, k))[0] * node->GetCoord(i, nj - 3, k)[1]);
                temp_eta[Idx(i, nj - 2, k)].x_xi_z = 0.5 * (3.0 * coef->GetX(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[2] - 1.0 * coef->GetX(Idx(i, nj - 3, k))[0] * node->GetCoord(i, nj - 3, k)[2]);
                temp_eta[Idx(i, nj - 2, k)].y_xi_x = 0.5 * (3.0 * coef->GetY(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[0] - 1.0 * coef->GetY(Idx(i, nj - 3, k))[0] * node->GetCoord(i, nj - 3, k)[0]);
            }
        }

        // 求出zeta方向的临时变量：逆变换度量系数乘以坐标
        struct TempZeta
        {
            double y_eta_z, z_eta_x, x_eta_y, y_xi_z, z_xi_x, x_xi_y, z_eta_y, x_eta_z, y_eta_x, z_xi_y, x_xi_z, y_xi_x;
        };
        std::vector<TempZeta> temp_zeta(ni * nj * nk);
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    int idx = Idx(i, j, k);
                    auto &temp = temp_zeta[idx];
                    temp.y_eta_z = 0.5 * (coef->GetY(Idx(i, j, k + 1))[1] * node->GetCoord(i, j, k + 1)[2] +
                                          coef->GetY(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[2]);
                    temp.z_eta_x = 0.5 * (coef->GetZ(Idx(i, j, k + 1))[1] * node->GetCoord(i, j, k + 1)[0] +
                                          coef->GetZ(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[0]);
                    temp.x_eta_y = 0.5 * (coef->GetX(Idx(i, j, k + 1))[1] * node->GetCoord(i, j, k + 1)[1] +
                                          coef->GetX(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[1]);
                    temp.y_xi_z = 0.5 * (coef->GetY(Idx(i, j, k + 1))[0] * node->GetCoord(i, j, k + 1)[2] +
                                         coef->GetY(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[2]);
                    temp.z_xi_x = 0.5 * (coef->GetZ(Idx(i, j, k + 1))[0] * node->GetCoord(i, j, k + 1)[0] +
                                         coef->GetZ(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[0]);
                    temp.x_xi_y = 0.5 * (coef->GetX(Idx(i, j, k + 1))[0] * node->GetCoord(i, j, k + 1)[1] +
                                         coef->GetX(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[1]);
                    temp.z_eta_y = 0.5 * (coef->GetZ(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[1] +
                                          coef->GetZ(Idx(i, j, k + 1))[1] * node->GetCoord(i, j, k + 1)[1]);
                    temp.x_eta_z = 0.5 * (coef->GetX(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[2] +
                                          coef->GetX(Idx(i, j, k + 1))[1] * node->GetCoord(i, j, k + 1)[2]);
                    temp.y_eta_x = 0.5 * (coef->GetY(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[0] +
                                          coef->GetY(Idx(i, j, k + 1))[1] * node->GetCoord(i, j, k + 1)[0]);
                    temp.z_xi_y = 0.5 * (coef->GetZ(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[1] +
                                         coef->GetZ(Idx(i, j, k + 1))[0] * node->GetCoord(i, j, k + 1)[1]);
                    temp.x_xi_z = 0.5 * (coef->GetX(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[2] +
                                         coef->GetX(Idx(i, j, k + 1))[0] * node->GetCoord(i, j, k + 1)[2]);
                    temp.y_xi_x = 0.5 * (coef->GetY(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[0] +
                                         coef->GetY(Idx(i, j, k + 1))[0] * node->GetCoord(i, j, k + 1)[0]);
                }
            }
        }
        // k=1/2；(nk-1)-1/2
        for (int j = 1; j < nj - 1; ++j)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                int idx = Idx(i, j, 0);
                temp_zeta[idx].y_eta_z = 0.5 * (3.0 * coef->GetY(Idx(i, j, 1))[1] * node->GetCoord(i, j, 1)[2] -
                                                1.0 * coef->GetY(Idx(i, j, 2))[1] * node->GetCoord(i, j, 2)[2]);
                temp_zeta[idx].z_eta_x = 0.5 * (3.0 * coef->GetZ(Idx(i, j, 1))[1] * node->GetCoord(i, j, 1)[0] -
                                                1.0 * coef->GetZ(Idx(i, j, 2))[1] * node->GetCoord(i, j, 2)[0]);
                temp_zeta[idx].x_eta_y = 0.5 * (3.0 * coef->GetX(Idx(i, j, 1))[1] * node->GetCoord(i, j, 1)[1] -
                                                1.0 * coef->GetX(Idx(i, j, 2))[1] * node->GetCoord(i, j, 2)[1]);
                temp_zeta[idx].y_xi_z = 0.5 * (3.0 * coef->GetY(Idx(i, j, 1))[0] * node->GetCoord(i, j, 1)[2] -
                                               1.0 * coef->GetY(Idx(i, j, 2))[0] * node->GetCoord(i, j, 2)[2]);
                temp_zeta[idx].z_xi_x = 0.5 * (3.0 * coef->GetZ(Idx(i, j, 1))[0] * node->GetCoord(i, j, 1)[0] -
                                               1.0 * coef->GetZ(Idx(i, j, 2))[0] * node->GetCoord(i, j, 2)[0]);
                temp_zeta[idx].x_xi_y = 0.5 * (3.0 * coef->GetX(Idx(i, j, 1))[0] * node->GetCoord(i, j, 1)[1] -
                                               1.0 * coef->GetX(Idx(i, j, 2))[0] * node->GetCoord(i, j, 2)[1]);
                temp_zeta[Idx(i, j, 0)].z_eta_y = 0.5 * (3.0 * coef->GetZ(Idx(i, j, 1))[1] * node->GetCoord(i, j, 1)[1] - 1.0 * coef->GetZ(Idx(i, j, 2))[1] * node->GetCoord(i, j, 2)[1]);
                temp_zeta[Idx(i, j, 0)].x_eta_z = 0.5 * (3.0 * coef->GetX(Idx(i, j, 1))[1] * node->GetCoord(i, j, 1)[2] - 1.0 * coef->GetX(Idx(i, j, 2))[1] * node->GetCoord(i, j, 2)[2]);
                temp_zeta[Idx(i, j, 0)].y_eta_x = 0.5 * (3.0 * coef->GetY(Idx(i, j, 1))[1] * node->GetCoord(i, j, 1)[0] - 1.0 * coef->GetY(Idx(i, j, 2))[1] * node->GetCoord(i, j, 2)[0]);
                temp_zeta[Idx(i, j, 0)].z_xi_y = 0.5 * (3.0 * coef->GetZ(Idx(i, j, 1))[0] * node->GetCoord(i, j, 1)[1] - 1.0 * coef->GetZ(Idx(i, j, 2))[0] * node->GetCoord(i, j, 2)[1]);
                temp_zeta[Idx(i, j, 0)].x_xi_z = 0.5 * (3.0 * coef->GetX(Idx(i, j, 1))[0] * node->GetCoord(i, j, 1)[2] - 1.0 * coef->GetX(Idx(i, j, 2))[0] * node->GetCoord(i, j, 2)[2]);
                temp_zeta[Idx(i, j, 0)].y_xi_x = 0.5 * (3.0 * coef->GetY(Idx(i, j, 1))[0] * node->GetCoord(i, j, 1)[0] - 1.0 * coef->GetY(Idx(i, j, 2))[0] * node->GetCoord(i, j, 2)[0]);
                idx = Idx(i, j, nk - 2);
                temp_zeta[idx].y_eta_z = 0.5 * (3.0 * coef->GetY(Idx(i, j, nk - 2))[1] * node->GetCoord(i, j, nk - 2)[2] -
                                                1.0 * coef->GetY(Idx(i, j, nk - 3))[1] * node->GetCoord(i, j, nk - 3)[2]);
                temp_zeta[idx].z_eta_x = 0.5 * (3.0 * coef->GetZ(Idx(i, j, nk - 2))[1] * node->GetCoord(i, j, nk - 2)[0] -
                                                1.0 * coef->GetZ(Idx(i, j, nk - 3))[1] * node->GetCoord(i, j, nk - 3)[0]);
                temp_zeta[idx].x_eta_y = 0.5 * (3.0 * coef->GetX(Idx(i, j, nk - 2))[1] * node->GetCoord(i, j, nk - 2)[1] -
                                                1.0 * coef->GetX(Idx(i, j, nk - 3))[1] * node->GetCoord(i, j, nk - 3)[1]);
                temp_zeta[idx].y_xi_z = 0.5 * (3.0 * coef->GetY(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[2] -
                                               1.0 * coef->GetY(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[2]);
                temp_zeta[idx].z_xi_x = 0.5 * (3.0 * coef->GetZ(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[0] -
                                               1.0 * coef->GetZ(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[0]);
                temp_zeta[idx].x_xi_y = 0.5 * (3.0 * coef->GetX(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[1] -
                                               1.0 * coef->GetX(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[1]);
                temp_zeta[Idx(i, j, nk - 2)].z_eta_y = 0.5 * (3.0 * coef->GetZ(Idx(i, j, nk - 2))[1] * node->GetCoord(i, j, nk - 2)[1] - 1.0 * coef->GetZ(Idx(i, j, nk - 3))[1] * node->GetCoord(i, j, nk - 3)[1]);
                temp_zeta[Idx(i, j, nk - 2)].x_eta_z = 0.5 * (3.0 * coef->GetX(Idx(i, j, nk - 2))[1] * node->GetCoord(i, j, nk - 2)[2] - 1.0 * coef->GetX(Idx(i, j, nk - 3))[1] * node->GetCoord(i, j, nk - 3)[2]);
                temp_zeta[Idx(i, j, nk - 2)].y_eta_x = 0.5 * (3.0 * coef->GetY(Idx(i, j, nk - 2))[1] * node->GetCoord(i, j, nk - 2)[0] - 1.0 * coef->GetY(Idx(i, j, nk - 3))[1] * node->GetCoord(i, j, nk - 3)[0]);
                temp_zeta[Idx(i, j, nk - 2)].z_xi_y = 0.5 * (3.0 * coef->GetZ(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[1] - 1.0 * coef->GetZ(Idx(i, j, nk - 3))[0] * node->GetCoord(i, j, nk - 3)[1]);
                temp_zeta[Idx(i, j, nk - 2)].x_xi_z = 0.5 * (3.0 * coef->GetX(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[2] - 1.0 * coef->GetX(Idx(i, j, nk - 3))[0] * node->GetCoord(i, j, nk - 3)[2]);
                temp_zeta[Idx(i, j, nk - 2)].y_xi_x = 0.5 * (3.0 * coef->GetY(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[0] - 1.0 * coef->GetY(Idx(i, j, nk - 3))[0] * node->GetCoord(i, j, nk - 3)[0]);
            }
        }

        // 第四步：根据半点坐标和半点逆变换度量系数使用守恒形式计算整点度量系数（CMM1)
        double temp[6][2];
        // 获取下标的lamda函数
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    // i direction
                    auto coef_xi = coef->GetXi(Idx(i, j, k));
                    auto coef_eta = coef->GetEta(Idx(i, j, k));
                    auto coef_zeta = coef->GetZeta(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 2; iTemp++)
                    {
                        int idx = Idx(i - 1 + iTemp, j, k);
                        temp[0][iTemp] = 0.5 * (temp_xi[idx].y_zeta_z - temp_xi[idx].z_zeta_y); // y_zeta*z-z_zeta*y
                        temp[1][iTemp] = 0.5 * (temp_xi[idx].z_zeta_x - temp_xi[idx].x_zeta_z); // z_zeta*x-x_zeta*z
                        temp[2][iTemp] = 0.5 * (temp_xi[idx].x_zeta_y - temp_xi[idx].y_zeta_x); // x_zeta*y-y_zeta*x
                        temp[3][iTemp] = 0.5 * (temp_xi[idx].z_eta_y - temp_xi[idx].y_eta_z);   // z_eta*y-y_eta*z
                        temp[4][iTemp] = 0.5 * (temp_xi[idx].x_eta_z - temp_xi[idx].z_eta_x);   // x_eta*z-z_eta*x
                        temp[5][iTemp] = 0.5 * (temp_xi[idx].y_eta_x - temp_xi[idx].x_eta_y);   // y_eta*x-x_eta*y
                    }
                    coef_eta[0] += temp[0][1] - temp[0][0];
                    coef_eta[1] += temp[1][1] - temp[1][0];
                    coef_eta[2] += temp[2][1] - temp[2][0];
                    coef_zeta[0] += temp[3][1] - temp[3][0];
                    coef_zeta[1] += temp[4][1] - temp[4][0];
                    coef_zeta[2] += temp[5][1] - temp[5][0];
                    for (int iTemp = 0; iTemp < 2; iTemp++)
                    {
                        int idx = Idx(i, j - 1 + iTemp, k);
                        temp[0][iTemp] = 0.5 * (temp_eta[idx].z_zeta_y - temp_eta[idx].y_zeta_z); // z_zeta*y-y_zeta*z
                        temp[1][iTemp] = 0.5 * (temp_eta[idx].x_zeta_z - temp_eta[idx].z_zeta_x); // x_zeta*z-z_zeta*x
                        temp[2][iTemp] = 0.5 * (temp_eta[idx].y_zeta_x - temp_eta[idx].x_zeta_y); //  y_zeta*x-x_zeta*y
                        temp[3][iTemp] = 0.5 * (temp_eta[idx].y_xi_z - temp_eta[idx].z_xi_y);     // y_xi*z-z_xi*y
                        temp[4][iTemp] = 0.5 * (temp_eta[idx].z_xi_x - temp_eta[idx].x_xi_z);     // z_xi*x-x_xi*z
                        temp[5][iTemp] = 0.5 * (temp_eta[idx].x_xi_y - temp_eta[idx].y_xi_x);     // x_xi*y-y_xi*x
                    }
                    coef_xi[0] += temp[0][1] - temp[0][0];
                    coef_xi[1] += temp[1][1] - temp[1][0];
                    coef_xi[2] += temp[2][1] - temp[2][0];
                    coef_zeta[0] += temp[3][1] - temp[3][0];
                    coef_zeta[1] += temp[4][1] - temp[4][0];
                    coef_zeta[2] += temp[5][1] - temp[5][0];
                    for (int iTemp = 0; iTemp < 2; iTemp++)
                    {
                        int idx = Idx(i, j, k - 1 + iTemp);
                        temp[0][iTemp] = 0.5 * (temp_zeta[idx].y_eta_z - temp_zeta[idx].z_eta_y); // y_eta*z-z_eta*y
                        temp[1][iTemp] = 0.5 * (temp_zeta[idx].z_eta_x - temp_zeta[idx].x_eta_z); // z_eta*x-x_eta*z
                        temp[2][iTemp] = 0.5 * (temp_zeta[idx].x_eta_y - temp_zeta[idx].y_eta_x); // x_eta*y-y_eta*x
                        temp[3][iTemp] = 0.5 * (temp_zeta[idx].z_xi_y - temp_zeta[idx].y_xi_z);   // z_xi*y-y_xi*z
                        temp[4][iTemp] = 0.5 * (temp_zeta[idx].x_xi_z - temp_zeta[idx].z_xi_x);   // x_xi*z-z_xi*x
                        temp[5][iTemp] = 0.5 * (temp_zeta[idx].y_xi_x - temp_zeta[idx].x_xi_y);   //  y_xi*x-x_xi*y
                    }
                    coef_xi[0] += temp[0][1] - temp[0][0];
                    coef_xi[1] += temp[1][1] - temp[1][0];
                    coef_xi[2] += temp[2][1] - temp[2][0];
                    coef_eta[0] += temp[3][1] - temp[3][0];
                    coef_eta[1] += temp[4][1] - temp[4][0];
                    coef_eta[2] += temp[5][1] - temp[5][0];
                    if (grid->GetDim() == 2)
                    {
                        coef_zeta[0] = 0.0;
                        coef_zeta[1] = 0.0;
                        coef_zeta[2] = 1.0;
                    }
                    coef_xi[3] = coef_eta[3] = coef_zeta[3] = 0.0;
                }
            }
        }
    }
    void NSSolverStruct::CalcMetricsS3_6th()
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
        double inter_temp[6];
        // i+1/2,j+1/2,k+1/2处的坐标
        std::vector<std::vector<double>> coord_i(ni * nj * nk, std::vector<double>(3)),
            coord_j(ni * nj * nk, std::vector<double>(3)), coord_k(ni * nj * nk, std::vector<double>(3));
        // 计算之前，先把度量系数赋值为0
        for (int i = 0; i < ni; ++i)
        {
            for (int j = 0; j < nj; ++j)
            {
                for (int k = 0; k < nk; ++k)
                {
                    for (int iDim = 0; iDim < 4; ++iDim)
                    {
                        coef->GetX(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetY(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetZ(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetXi(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetEta(Idx(i, j, k))[iDim] = 0.0;
                        coef->GetZeta(Idx(i, j, k))[iDim] = 0.0;
                        for (int jDim = 0; jDim < 3; ++jDim)
                        {
                            coef_mid[jDim]->GetX(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetY(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetZ(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetXi(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetEta(Idx(i, j, k))[iDim] = 0.0;
                            coef_mid[jDim]->GetZeta(Idx(i, j, k))[iDim] = 0.0;
                        }
                    }
                }
            }
        }
        // 第一步：计算半点坐标
        // i+1/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 3; i < ni - 4; ++i)
                {
                    for (int iDim = 0; iDim < 3; ++iDim)
                    {
                        for (int iTemp = 0; iTemp < 6; iTemp++)
                        {
                            inter_temp[iTemp] = node->GetCoord(i + iTemp - 2, j, k)[iDim];
                        }
                        coord_i[Idx(i, j, k)][iDim] = MidNodeInter6th(inter_temp);
                    }
                }
            }
        }
        // i=1/2,3/2,5/2;ni-3/2,ni-5/2,ni-7/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int iDim = 0; iDim < 3; ++iDim)
                {
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(1 + iTemp, j, k)[iDim];
                    }
                    coord_i[Idx(0, j, k)][iDim] = MidNodeInter4thRight1(inter_temp);
                    coord_i[Idx(1, j, k)][iDim] = MidNodeInter4thRight2(inter_temp);
                    coord_i[Idx(2, j, k)][iDim] = MidNodeInter4th(inter_temp);
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(ni - 5 + iTemp, j, k)[iDim];
                    }
                    coord_i[Idx(ni - 2, j, k)][iDim] = MidNodeInter4thLeft1(inter_temp);
                    coord_i[Idx(ni - 3, j, k)][iDim] = MidNodeInter4thLeft2(inter_temp);
                    coord_i[Idx(ni - 4, j, k)][iDim] = MidNodeInter4th(inter_temp);
                }
            }
        }
        // j+1/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 3; j < nj - 4; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    for (int iDim = 0; iDim < 3; ++iDim)
                    {
                        for (int iTemp = 0; iTemp < 6; iTemp++)
                        {
                            inter_temp[iTemp] = node->GetCoord(i, j + iTemp - 2, k)[iDim];
                        }
                        coord_j[Idx(i, j, k)][iDim] = MidNodeInter6th(inter_temp);
                    }
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
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(i, 1 + iTemp, k)[iDim];
                    }
                    coord_j[Idx(i, 0, k)][iDim] = MidNodeInter4thRight1(inter_temp);
                    coord_j[Idx(i, 1, k)][iDim] = MidNodeInter4thRight2(inter_temp);
                    coord_j[Idx(i, 2, k)][iDim] = MidNodeInter4th(inter_temp);
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(i, nj - 5 + iTemp, k)[iDim];
                    }
                    coord_j[Idx(i, nj - 2, k)][iDim] = MidNodeInter4thLeft1(inter_temp);
                    coord_j[Idx(i, nj - 3, k)][iDim] = MidNodeInter4thLeft2(inter_temp);
                    coord_j[Idx(i, nj - 4, k)][iDim] = MidNodeInter4th(inter_temp);
                }
            }
        }
        // k+1/2
        for (int k = 3; k < nk - 4; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    for (int iDim = 0; iDim < 3; ++iDim)
                    {
                        for (int iTemp = 0; iTemp < 6; iTemp++)
                        {
                            inter_temp[iTemp] = node->GetCoord(i, j, k + iTemp - 2)[iDim];
                        }
                        coord_k[Idx(i, j, k)][iDim] = MidNodeInter6th(inter_temp);
                    }
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
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(i, j, 1 + iTemp)[iDim];
                    }
                    coord_k[Idx(i, j, 0)][iDim] = MidNodeInter4thRight1(inter_temp);
                    coord_k[Idx(i, j, 1)][iDim] = MidNodeInter4thRight2(inter_temp);
                    coord_k[Idx(i, j, 2)][iDim] = MidNodeInter4th(inter_temp);
                    for (int iTemp = 0; iTemp < 4; iTemp++)
                    {
                        inter_temp[iTemp] = node->GetCoord(i, j, nk - 5 + iTemp)[iDim];
                    }
                    coord_k[Idx(i, j, nk - 2)][iDim] = MidNodeInter4thLeft1(inter_temp);
                    coord_k[Idx(i, j, nk - 3)][iDim] = MidNodeInter4thLeft2(inter_temp);
                    coord_k[Idx(i, j, nk - 4)][iDim] = MidNodeInter4th(inter_temp);
                }
            }
        }
        // 第二步：根据半点坐标计算整点逆变换度量系数
        double diff_temp[6];
        // i direction
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 3; i < ni - 3; ++i)
                {
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_i[Idx(i + iTemp - 3, j, k)][0];
                    }
                    auto coef_x = coef->GetX(Idx(i, j, k));
                    coef_x[0] = NodeDifferece6th(diff_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_i[Idx(i + iTemp - 3, j, k)][1];
                    }
                    auto coef_y = coef->GetY(Idx(i, j, k));
                    coef_y[0] = NodeDifferece6th(diff_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_i[Idx(i + iTemp - 3, j, k)][2];
                    }
                    auto coef_z = coef->GetZ(Idx(i, j, k));
                    coef_z[0] = NodeDifferece6th(diff_temp);
                }
            }
        }
        // j direction
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 3; j < nj - 3; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    auto coef_x = coef->GetX(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_j[Idx(i, j + iTemp - 3, k)][0];
                    }
                    coef_x[1] = NodeDifferece6th(diff_temp);
                    auto coef_y = coef->GetY(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_j[Idx(i, j + iTemp - 3, k)][1];
                    }
                    coef_y[1] = NodeDifferece6th(diff_temp);
                    auto coef_z = coef->GetZ(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_j[Idx(i, j + iTemp - 3, k)][2];
                    }
                    coef_z[1] = NodeDifferece6th(diff_temp);
                }
            }
        }
        // k direction
        for (int k = 3; k < nk - 3; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    auto coef_x = coef->GetX(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_k[Idx(i, j, k + iTemp - 3)][0];
                    }
                    coef_x[2] = NodeDifferece6th(diff_temp);
                    auto coef_y = coef->GetY(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_k[Idx(i, j, k + iTemp - 3)][1];
                    }
                    coef_y[2] = NodeDifferece6th(diff_temp);
                    auto coef_z = coef->GetZ(Idx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        diff_temp[iTemp] = coord_k[Idx(i, j, k + iTemp - 3)][2];
                    }
                    coef_z[2] = NodeDifferece6th(diff_temp);
                    if (grid->GetDim() == 2)
                    {
                        coef_x[2] = 0.0;
                        coef_y[2] = 0.0;
                        coef_z[2] = 1.0;
                    }
                }
            }
        }
        // 计算边界整点的度量系数
        // i=1,2;ni-3,ni-2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_i[Idx(0 + iTemp, j, k)][0];
                }
                coef->GetX(Idx(1, j, k))[0] = NodeDifferece4thRight(diff_temp);
                coef->GetX(Idx(2, j, k))[0] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_i[Idx(0 + iTemp, j, k)][1];
                }
                coef->GetY(Idx(1, j, k))[0] = NodeDifferece4thRight(diff_temp);
                coef->GetY(Idx(2, j, k))[0] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_i[Idx(0 + iTemp, j, k)][2];
                }
                coef->GetZ(Idx(1, j, k))[0] = NodeDifferece4thRight(diff_temp);
                coef->GetZ(Idx(2, j, k))[0] = NodeDifferece4th(diff_temp);

                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_i[Idx(ni - 2 - iTemp, j, k)][0];
                }
                coef->GetX(Idx(ni - 2, j, k))[0] = NodeDifferece4thLeft(diff_temp);
                coef->GetX(Idx(ni - 3, j, k))[0] = -NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_i[Idx(ni - 2 - iTemp, j, k)][1];
                }
                coef->GetY(Idx(ni - 2, j, k))[0] = NodeDifferece4thLeft(diff_temp);
                coef->GetY(Idx(ni - 3, j, k))[0] = -NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_i[Idx(ni - 2 - iTemp, j, k)][2];
                }
                coef->GetZ(Idx(ni - 2, j, k))[0] = NodeDifferece4thLeft(diff_temp);
                coef->GetZ(Idx(ni - 3, j, k))[0] = -NodeDifferece4th(diff_temp);
            }
        }
        // j=1,2;nj-3,nj-2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_j[Idx(i, 0 + iTemp, k)][0];
                }
                coef->GetX(Idx(i, 1, k))[1] = NodeDifferece4thRight(diff_temp);
                coef->GetX(Idx(i, 2, k))[1] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_j[Idx(i, 0 + iTemp, k)][1];
                }
                coef->GetY(Idx(i, 1, k))[1] = NodeDifferece4thRight(diff_temp);
                coef->GetY(Idx(i, 2, k))[1] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_j[Idx(i, 0 + iTemp, k)][2];
                }
                coef->GetZ(Idx(i, 1, k))[1] = NodeDifferece4thRight(diff_temp);
                coef->GetZ(Idx(i, 2, k))[1] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_j[Idx(i, nj - 2 - iTemp, k)][0];
                }
                coef->GetX(Idx(i, nj - 2, k))[1] = NodeDifferece4thLeft(diff_temp);
                coef->GetX(Idx(i, nj - 3, k))[1] = -NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_j[Idx(i, nj - 2 - iTemp, k)][1];
                }
                coef->GetY(Idx(i, nj - 2, k))[1] = NodeDifferece4thLeft(diff_temp);
                coef->GetY(Idx(i, nj - 3, k))[1] = -NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_j[Idx(i, nj - 2 - iTemp, k)][2];
                }
                coef->GetZ(Idx(i, nj - 2, k))[1] = NodeDifferece4thLeft(diff_temp);
                coef->GetZ(Idx(i, nj - 3, k))[1] = -NodeDifferece4th(diff_temp);
            }
        }
        // k=1,2;nk-3,nk-2
        for (int j = 1; j < nj - 1; ++j)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_k[Idx(i, j, 0 + iTemp)][0];
                }
                coef->GetX(Idx(i, j, 1))[2] = NodeDifferece4thRight(diff_temp);
                coef->GetX(Idx(i, j, 2))[2] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_k[Idx(i, j, 0 + iTemp)][1];
                }
                coef->GetY(Idx(i, j, 1))[2] = NodeDifferece4thRight(diff_temp);
                coef->GetY(Idx(i, j, 2))[2] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_k[Idx(i, j, 0 + iTemp)][2];
                }
                coef->GetZ(Idx(i, j, 1))[2] = NodeDifferece4thRight(diff_temp);
                coef->GetZ(Idx(i, j, 2))[2] = NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_k[Idx(i, j, nk - 2 - iTemp)][0];
                }
                coef->GetX(Idx(i, j, nk - 2))[2] = NodeDifferece4thLeft(diff_temp);
                coef->GetX(Idx(i, j, nk - 3))[2] = -NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_k[Idx(i, j, nk - 2 - iTemp)][1];
                }
                coef->GetY(Idx(i, j, nk - 2))[2] = NodeDifferece4thLeft(diff_temp);
                coef->GetY(Idx(i, j, nk - 3))[2] = -NodeDifferece4th(diff_temp);
                for (int iTemp = 0; iTemp < 6; iTemp++)
                {
                    diff_temp[iTemp] = coord_k[Idx(i, j, nk - 2 - iTemp)][2];
                }
                coef->GetZ(Idx(i, j, nk - 2))[2] = NodeDifferece4thLeft(diff_temp);
                coef->GetZ(Idx(i, j, nk - 3))[2] = -NodeDifferece4th(diff_temp);
                if (grid->GetDim() == 2)
                {
                    coef->GetZ(Idx(i, j, 1))[0] = coef->GetZ(Idx(i, j, 1))[1] = 0.0;
                    coef->GetZ(Idx(i, j, 2))[0] = coef->GetZ(Idx(i, j, 2))[1] = 0.0;
                    coef->GetZ(Idx(i, j, nk - 2))[0] = coef->GetZ(Idx(i, j, nk - 2))[1] = 0.0;
                    coef->GetZ(Idx(i, j, nk - 3))[0] = coef->GetZ(Idx(i, j, nk - 3))[1] = 0.0;
                    coef->GetZ(Idx(i, j, 1))[2] = 1.0;
                    coef->GetZ(Idx(i, j, 2))[2] = 1.0;
                    coef->GetZ(Idx(i, j, nk - 2))[2] = 1.0;
                    coef->GetZ(Idx(i, j, nk - 3))[2] = 1.0;
                }
            }
        }

        // 第三步：计算半点处逆变换度量系数和坐标乘积
        //  求出xi方向的临时变量：逆变换度量系数乘以坐标
        struct TempXi
        {
            double y_zeta_z, z_zeta_x, x_zeta_y, y_eta_z, z_eta_x, x_eta_y, z_zeta_y, x_zeta_z, y_zeta_x, z_eta_y, x_eta_z, y_eta_x;
        };
        std::vector<TempXi> temp_xi(ni * nj * nk);
        // i+1/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 3; i < ni - 4; ++i)
                {
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetY(Idx(i + iTemp - 2, j, k))[2] * node->GetCoord(i + iTemp - 2, j, k)[2];
                    }
                    temp_xi[Idx(i, j, k)].y_zeta_z = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZ(Idx(i + iTemp - 2, j, k))[2] * node->GetCoord(i + iTemp - 2, j, k)[0];
                    }
                    temp_xi[Idx(i, j, k)].z_zeta_x = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetX(Idx(i + iTemp - 2, j, k))[2] * node->GetCoord(i + iTemp - 2, j, k)[1];
                    }
                    temp_xi[Idx(i, j, k)].x_zeta_y = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetY(Idx(i + iTemp - 2, j, k))[1] * node->GetCoord(i + iTemp - 2, j, k)[2];
                    }
                    temp_xi[Idx(i, j, k)].y_eta_z = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZ(Idx(i + iTemp - 2, j, k))[1] * node->GetCoord(i + iTemp - 2, j, k)[0];
                    }
                    temp_xi[Idx(i, j, k)].z_eta_x = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetX(Idx(i + iTemp - 2, j, k))[1] * node->GetCoord(i + iTemp - 2, j, k)[1];
                    }
                    temp_xi[Idx(i, j, k)].x_eta_y = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZ(Idx(i - 2 + iTemp, j, k))[2] * node->GetCoord(i - 2 + iTemp, j, k)[1];
                    }
                    temp_xi[Idx(i, j, k)].z_zeta_y = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetX(Idx(i - 2 + iTemp, j, k))[2] * node->GetCoord(i - 2 + iTemp, j, k)[2];
                    }
                    temp_xi[Idx(i, j, k)].x_zeta_z = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetY(Idx(i - 2 + iTemp, j, k))[2] * node->GetCoord(i - 2 + iTemp, j, k)[0];
                    }
                    temp_xi[Idx(i, j, k)].y_zeta_x = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZ(Idx(i - 2 + iTemp, j, k))[1] * node->GetCoord(i - 2 + iTemp, j, k)[1];
                    }
                    temp_xi[Idx(i, j, k)].z_eta_y = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetX(Idx(i - 2 + iTemp, j, k))[1] * node->GetCoord(i - 2 + iTemp, j, k)[2];
                    }
                    temp_xi[Idx(i, j, k)].x_eta_z = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetY(Idx(i - 2 + iTemp, j, k))[1] * node->GetCoord(i - 2 + iTemp, j, k)[0];
                    }
                    temp_xi[Idx(i, j, k)].y_eta_x = MidNodeInter6th(inter_temp);
                }
            }
        }
        // i=1/2,3/2,5/2;ni-3/2,ni-5/2,ni-7/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                // i=1/2，3/2，5/2
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(1 + iTemp, j, k))[2] * node->GetCoord(1 + iTemp, j, k)[2];
                }
                temp_xi[Idx(0, j, k)].y_zeta_z = MidNodeInter4thRight1(inter_temp);
                temp_xi[Idx(1, j, k)].y_zeta_z = MidNodeInter4thRight2(inter_temp);
                temp_xi[Idx(2, j, k)].y_zeta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(1 + iTemp, j, k))[2] * node->GetCoord(1 + iTemp, j, k)[0];
                }
                temp_xi[Idx(0, j, k)].z_zeta_x = MidNodeInter4thRight1(inter_temp);
                temp_xi[Idx(1, j, k)].z_zeta_x = MidNodeInter4thRight2(inter_temp);
                temp_xi[Idx(2, j, k)].z_zeta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(1 + iTemp, j, k))[2] * node->GetCoord(1 + iTemp, j, k)[1];
                }
                temp_xi[Idx(0, j, k)].x_zeta_y = MidNodeInter4thRight1(inter_temp);
                temp_xi[Idx(1, j, k)].x_zeta_y = MidNodeInter4thRight2(inter_temp);
                temp_xi[Idx(2, j, k)].x_zeta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(1 + iTemp, j, k))[1] * node->GetCoord(1 + iTemp, j, k)[2];
                }
                temp_xi[Idx(0, j, k)].y_eta_z = MidNodeInter4thRight1(inter_temp);
                temp_xi[Idx(1, j, k)].y_eta_z = MidNodeInter4thRight2(inter_temp);
                temp_xi[Idx(2, j, k)].y_eta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(1 + iTemp, j, k))[1] * node->GetCoord(1 + iTemp, j, k)[0];
                }
                temp_xi[Idx(0, j, k)].z_eta_x = MidNodeInter4thRight1(inter_temp);
                temp_xi[Idx(1, j, k)].z_eta_x = MidNodeInter4thRight2(inter_temp);
                temp_xi[Idx(2, j, k)].z_eta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(1 + iTemp, j, k))[1] * node->GetCoord(1 + iTemp, j, k)[1];
                }
                temp_xi[Idx(0, j, k)].x_eta_y = MidNodeInter4thRight1(inter_temp);
                temp_xi[Idx(1, j, k)].x_eta_y = MidNodeInter4thRight2(inter_temp);
                temp_xi[Idx(2, j, k)].x_eta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(1 + iTemp, j, k))[2] * node->GetCoord(1 + iTemp, j, k)[1];
                }
                temp_xi[Idx(0, j, k)].z_zeta_y = MidNodeInter4thRight1(inter_temp);
                temp_xi[Idx(1, j, k)].z_zeta_y = MidNodeInter4thRight2(inter_temp);
                temp_xi[Idx(2, j, k)].z_zeta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(1 + iTemp, j, k))[2] * node->GetCoord(1 + iTemp, j, k)[2];
                }
                temp_xi[Idx(0, j, k)].x_zeta_z = MidNodeInter4thRight1(inter_temp);
                temp_xi[Idx(1, j, k)].x_zeta_z = MidNodeInter4thRight2(inter_temp);
                temp_xi[Idx(2, j, k)].x_zeta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(1 + iTemp, j, k))[2] * node->GetCoord(1 + iTemp, j, k)[0];
                }
                temp_xi[Idx(0, j, k)].y_zeta_x = MidNodeInter4thRight1(inter_temp);
                temp_xi[Idx(1, j, k)].y_zeta_x = MidNodeInter4thRight2(inter_temp);
                temp_xi[Idx(2, j, k)].y_zeta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(1 + iTemp, j, k))[1] * node->GetCoord(1 + iTemp, j, k)[1];
                }
                temp_xi[Idx(0, j, k)].z_eta_y = MidNodeInter4thRight1(inter_temp);
                temp_xi[Idx(1, j, k)].z_eta_y = MidNodeInter4thRight2(inter_temp);
                temp_xi[Idx(2, j, k)].z_eta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(1 + iTemp, j, k))[1] * node->GetCoord(1 + iTemp, j, k)[2];
                }
                temp_xi[Idx(0, j, k)].x_eta_z = MidNodeInter4thRight1(inter_temp);
                temp_xi[Idx(1, j, k)].x_eta_z = MidNodeInter4thRight2(inter_temp);
                temp_xi[Idx(2, j, k)].x_eta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(1 + iTemp, j, k))[1] * node->GetCoord(1 + iTemp, j, k)[0];
                }
                temp_xi[Idx(0, j, k)].y_eta_x = MidNodeInter4thRight1(inter_temp);
                temp_xi[Idx(1, j, k)].y_eta_x = MidNodeInter4thRight2(inter_temp);
                temp_xi[Idx(2, j, k)].y_eta_x = MidNodeInter4th(inter_temp);
                // i=ni-3/2,ni-5/2,ni-7/2
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(ni - 5 + iTemp, j, k))[2] * node->GetCoord(ni - 5 + iTemp, j, k)[2];
                }
                temp_xi[Idx(ni - 2, j, k)].y_zeta_z = MidNodeInter4thLeft1(inter_temp);
                temp_xi[Idx(ni - 3, j, k)].y_zeta_z = MidNodeInter4thLeft2(inter_temp);
                temp_xi[Idx(ni - 4, j, k)].y_zeta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(ni - 5 + iTemp, j, k))[2] * node->GetCoord(ni - 5 + iTemp, j, k)[0];
                }
                temp_xi[Idx(ni - 2, j, k)].z_zeta_x = MidNodeInter4thLeft1(inter_temp);
                temp_xi[Idx(ni - 3, j, k)].z_zeta_x = MidNodeInter4thLeft2(inter_temp);
                temp_xi[Idx(ni - 4, j, k)].z_zeta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(ni - 5 + iTemp, j, k))[2] * node->GetCoord(ni - 5 + iTemp, j, k)[1];
                }
                temp_xi[Idx(ni - 2, j, k)].x_zeta_y = MidNodeInter4thLeft1(inter_temp);
                temp_xi[Idx(ni - 3, j, k)].x_zeta_y = MidNodeInter4thLeft2(inter_temp);
                temp_xi[Idx(ni - 4, j, k)].x_zeta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(ni - 5 + iTemp, j, k))[1] * node->GetCoord(ni - 5 + iTemp, j, k)[2];
                }
                temp_xi[Idx(ni - 2, j, k)].y_eta_z = MidNodeInter4thLeft1(inter_temp);
                temp_xi[Idx(ni - 3, j, k)].y_eta_z = MidNodeInter4thLeft2(inter_temp);
                temp_xi[Idx(ni - 4, j, k)].y_eta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(ni - 5 + iTemp, j, k))[1] * node->GetCoord(ni - 5 + iTemp, j, k)[0];
                }
                temp_xi[Idx(ni - 2, j, k)].z_eta_x = MidNodeInter4thLeft1(inter_temp);
                temp_xi[Idx(ni - 3, j, k)].z_eta_x = MidNodeInter4thLeft2(inter_temp);
                temp_xi[Idx(ni - 4, j, k)].z_eta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(ni - 5 + iTemp, j, k))[1] * node->GetCoord(ni - 5 + iTemp, j, k)[1];
                }
                temp_xi[Idx(ni - 2, j, k)].x_eta_y = MidNodeInter4thLeft1(inter_temp);
                temp_xi[Idx(ni - 3, j, k)].x_eta_y = MidNodeInter4thLeft2(inter_temp);
                temp_xi[Idx(ni - 4, j, k)].x_eta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(ni - 5 + iTemp, j, k))[2] * node->GetCoord(ni - 5 + iTemp, j, k)[1];
                }
                temp_xi[Idx(ni - 2, j, k)].z_zeta_y = MidNodeInter4thLeft1(inter_temp);
                temp_xi[Idx(ni - 3, j, k)].z_zeta_y = MidNodeInter4thLeft2(inter_temp);
                temp_xi[Idx(ni - 4, j, k)].z_zeta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(ni - 5 + iTemp, j, k))[2] * node->GetCoord(ni - 5 + iTemp, j, k)[2];
                }
                temp_xi[Idx(ni - 2, j, k)].x_zeta_z = MidNodeInter4thLeft1(inter_temp);
                temp_xi[Idx(ni - 3, j, k)].x_zeta_z = MidNodeInter4thLeft2(inter_temp);
                temp_xi[Idx(ni - 4, j, k)].x_zeta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(ni - 5 + iTemp, j, k))[2] * node->GetCoord(ni - 5 + iTemp, j, k)[0];
                }
                temp_xi[Idx(ni - 2, j, k)].y_zeta_x = MidNodeInter4thLeft1(inter_temp);
                temp_xi[Idx(ni - 3, j, k)].y_zeta_x = MidNodeInter4thLeft2(inter_temp);
                temp_xi[Idx(ni - 4, j, k)].y_zeta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(ni - 5 + iTemp, j, k))[1] * node->GetCoord(ni - 5 + iTemp, j, k)[1];
                }
                temp_xi[Idx(ni - 2, j, k)].z_eta_y = MidNodeInter4thLeft1(inter_temp);
                temp_xi[Idx(ni - 3, j, k)].z_eta_y = MidNodeInter4thLeft2(inter_temp);
                temp_xi[Idx(ni - 4, j, k)].z_eta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(ni - 5 + iTemp, j, k))[1] * node->GetCoord(ni - 5 + iTemp, j, k)[2];
                }
                temp_xi[Idx(ni - 2, j, k)].x_eta_z = MidNodeInter4thLeft1(inter_temp);
                temp_xi[Idx(ni - 3, j, k)].x_eta_z = MidNodeInter4thLeft2(inter_temp);
                temp_xi[Idx(ni - 4, j, k)].x_eta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(ni - 5 + iTemp, j, k))[1] * node->GetCoord(ni - 5 + iTemp, j, k)[0];
                }
                temp_xi[Idx(ni - 2, j, k)].y_eta_x = MidNodeInter4thLeft1(inter_temp);
                temp_xi[Idx(ni - 3, j, k)].y_eta_x = MidNodeInter4thLeft2(inter_temp);
                temp_xi[Idx(ni - 4, j, k)].y_eta_x = MidNodeInter4th(inter_temp);
            }
        }
        // 求出eta方向的临时变量：逆变换度量系数乘以坐标
        struct TempEta
        {
            double y_zeta_z, z_zeta_x, x_zeta_y, y_xi_z, z_xi_x, x_xi_y, z_zeta_y, x_zeta_z, y_zeta_x, z_xi_y, x_xi_z, y_xi_x;
        };
        std::vector<TempEta> temp_eta(ni * nj * nk);
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 3; j < nj - 4; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetY(Idx(i, j + iTemp - 2, k))[2] * node->GetCoord(i, j + iTemp - 2, k)[2];
                    }
                    temp_eta[Idx(i, j, k)].y_zeta_z = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZ(Idx(i, j + iTemp - 2, k))[2] * node->GetCoord(i, j + iTemp - 2, k)[0];
                    }
                    temp_eta[Idx(i, j, k)].z_zeta_x = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetX(Idx(i, j + iTemp - 2, k))[2] * node->GetCoord(i, j + iTemp - 2, k)[1];
                    }
                    temp_eta[Idx(i, j, k)].x_zeta_y = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetY(Idx(i, j + iTemp - 2, k))[0] * node->GetCoord(i, j + iTemp - 2, k)[2];
                    }
                    temp_eta[Idx(i, j, k)].y_xi_z = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZ(Idx(i, j + iTemp - 2, k))[0] * node->GetCoord(i, j + iTemp - 2, k)[0];
                    }
                    temp_eta[Idx(i, j, k)].z_xi_x = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetX(Idx(i, j + iTemp - 2, k))[0] * node->GetCoord(i, j + iTemp - 2, k)[1];
                    }
                    temp_eta[Idx(i, j, k)].x_xi_y = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZ(Idx(i, j - 2 + iTemp, k))[2] * node->GetCoord(i, j - 2 + iTemp, k)[1];
                    }
                    temp_eta[Idx(i, j, k)].z_zeta_y = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetX(Idx(i, j - 2 + iTemp, k))[2] * node->GetCoord(i, j - 2 + iTemp, k)[2];
                    }
                    temp_eta[Idx(i, j, k)].x_zeta_z = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetY(Idx(i, j - 2 + iTemp, k))[2] * node->GetCoord(i, j - 2 + iTemp, k)[0];
                    }
                    temp_eta[Idx(i, j, k)].y_zeta_x = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZ(Idx(i, j - 2 + iTemp, k))[0] * node->GetCoord(i, j - 2 + iTemp, k)[1];
                    }
                    temp_eta[Idx(i, j, k)].z_xi_y = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetX(Idx(i, j - 2 + iTemp, k))[0] * node->GetCoord(i, j - 2 + iTemp, k)[2];
                    }
                    temp_eta[Idx(i, j, k)].x_xi_z = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetY(Idx(i, j - 2 + iTemp, k))[0] * node->GetCoord(i, j - 2 + iTemp, k)[0];
                    }
                    temp_eta[Idx(i, j, k)].y_xi_x = MidNodeInter6th(inter_temp);
                }
            }
        }
        // j=1/2,3/2,5/2;nj-3/2,nj-5/2,nj-7/2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                // j=1/2，3/2，5/2
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, 1 + iTemp, k))[2] * node->GetCoord(i, 1 + iTemp, k)[2];
                }
                temp_eta[Idx(i, 0, k)].y_zeta_z = MidNodeInter4thRight1(inter_temp);
                temp_eta[Idx(i, 1, k)].y_zeta_z = MidNodeInter4thRight2(inter_temp);
                temp_eta[Idx(i, 2, k)].y_zeta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, 1 + iTemp, k))[2] * node->GetCoord(i, 1 + iTemp, k)[0];
                }
                temp_eta[Idx(i, 0, k)].z_zeta_x = MidNodeInter4thRight1(inter_temp);
                temp_eta[Idx(i, 1, k)].z_zeta_x = MidNodeInter4thRight2(inter_temp);
                temp_eta[Idx(i, 2, k)].z_zeta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, 1 + iTemp, k))[2] * node->GetCoord(i, 1 + iTemp, k)[1];
                }
                temp_eta[Idx(i, 0, k)].x_zeta_y = MidNodeInter4thRight1(inter_temp);
                temp_eta[Idx(i, 1, k)].x_zeta_y = MidNodeInter4thRight2(inter_temp);
                temp_eta[Idx(i, 2, k)].x_zeta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, 1 + iTemp, k))[0] * node->GetCoord(i, 1 + iTemp, k)[2];
                }
                temp_eta[Idx(i, 0, k)].y_xi_z = MidNodeInter4thRight1(inter_temp);
                temp_eta[Idx(i, 1, k)].y_xi_z = MidNodeInter4thRight2(inter_temp);
                temp_eta[Idx(i, 2, k)].y_xi_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, 1 + iTemp, k))[0] * node->GetCoord(i, 1 + iTemp, k)[0];
                }
                temp_eta[Idx(i, 0, k)].z_xi_x = MidNodeInter4thRight1(inter_temp);
                temp_eta[Idx(i, 1, k)].z_xi_x = MidNodeInter4thRight2(inter_temp);
                temp_eta[Idx(i, 2, k)].z_xi_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, 1 + iTemp, k))[0] * node->GetCoord(i, 1 + iTemp, k)[1];
                }
                temp_eta[Idx(i, 0, k)].x_xi_y = MidNodeInter4thRight1(inter_temp);
                temp_eta[Idx(i, 1, k)].x_xi_y = MidNodeInter4thRight2(inter_temp);
                temp_eta[Idx(i, 2, k)].x_xi_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, 1 + iTemp, k))[2] * node->GetCoord(i, 1 + iTemp, k)[1];
                }
                temp_eta[Idx(i, 0, k)].z_zeta_y = MidNodeInter4thRight1(inter_temp);
                temp_eta[Idx(i, 1, k)].z_zeta_y = MidNodeInter4thRight2(inter_temp);
                temp_eta[Idx(i, 2, k)].z_zeta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, 1 + iTemp, k))[2] * node->GetCoord(i, 1 + iTemp, k)[2];
                }
                temp_eta[Idx(i, 0, k)].x_zeta_z = MidNodeInter4thRight1(inter_temp);
                temp_eta[Idx(i, 1, k)].x_zeta_z = MidNodeInter4thRight2(inter_temp);
                temp_eta[Idx(i, 2, k)].x_zeta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, 1 + iTemp, k))[2] * node->GetCoord(i, 1 + iTemp, k)[0];
                }
                temp_eta[Idx(i, 0, k)].y_zeta_x = MidNodeInter4thRight1(inter_temp);
                temp_eta[Idx(i, 1, k)].y_zeta_x = MidNodeInter4thRight2(inter_temp);
                temp_eta[Idx(i, 2, k)].y_zeta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, 1 + iTemp, k))[0] * node->GetCoord(i, 1 + iTemp, k)[1];
                }
                temp_eta[Idx(i, 0, k)].z_xi_y = MidNodeInter4thRight1(inter_temp);
                temp_eta[Idx(i, 1, k)].z_xi_y = MidNodeInter4thRight2(inter_temp);
                temp_eta[Idx(i, 2, k)].z_xi_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, 1 + iTemp, k))[0] * node->GetCoord(i, 1 + iTemp, k)[2];
                }
                temp_eta[Idx(i, 0, k)].x_xi_z = MidNodeInter4thRight1(inter_temp);
                temp_eta[Idx(i, 1, k)].x_xi_z = MidNodeInter4thRight2(inter_temp);
                temp_eta[Idx(i, 2, k)].x_xi_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, 1 + iTemp, k))[0] * node->GetCoord(i, 1 + iTemp, k)[0];
                }
                temp_eta[Idx(i, 0, k)].y_xi_x = MidNodeInter4thRight1(inter_temp);
                temp_eta[Idx(i, 1, k)].y_xi_x = MidNodeInter4thRight2(inter_temp);
                temp_eta[Idx(i, 2, k)].y_xi_x = MidNodeInter4th(inter_temp);
                // j=nj-3/2,nj-5/2,nj-7/2
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, nj - 5 + iTemp, k))[2] * node->GetCoord(i, nj - 5 + iTemp, k)[2];
                }
                temp_eta[Idx(i, nj - 2, k)].y_zeta_z = MidNodeInter4thLeft1(inter_temp);
                temp_eta[Idx(i, nj - 3, k)].y_zeta_z = MidNodeInter4thLeft2(inter_temp);
                temp_eta[Idx(i, nj - 4, k)].y_zeta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, nj - 5 + iTemp, k))[2] * node->GetCoord(i, nj - 5 + iTemp, k)[0];
                }
                temp_eta[Idx(i, nj - 2, k)].z_zeta_x = MidNodeInter4thLeft1(inter_temp);
                temp_eta[Idx(i, nj - 3, k)].z_zeta_x = MidNodeInter4thLeft2(inter_temp);
                temp_eta[Idx(i, nj - 4, k)].z_zeta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, nj - 5 + iTemp, k))[2] * node->GetCoord(i, nj - 5 + iTemp, k)[1];
                }
                temp_eta[Idx(i, nj - 2, k)].x_zeta_y = MidNodeInter4thLeft1(inter_temp);
                temp_eta[Idx(i, nj - 3, k)].x_zeta_y = MidNodeInter4thLeft2(inter_temp);
                temp_eta[Idx(i, nj - 4, k)].x_zeta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, nj - 5 + iTemp, k))[0] * node->GetCoord(i, nj - 5 + iTemp, k)[2];
                }
                temp_eta[Idx(i, nj - 2, k)].y_xi_z = MidNodeInter4thLeft1(inter_temp);
                temp_eta[Idx(i, nj - 3, k)].y_xi_z = MidNodeInter4thLeft2(inter_temp);
                temp_eta[Idx(i, nj - 4, k)].y_xi_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, nj - 5 + iTemp, k))[0] * node->GetCoord(i, nj - 5 + iTemp, k)[0];
                }
                temp_eta[Idx(i, nj - 2, k)].z_xi_x = MidNodeInter4thLeft1(inter_temp);
                temp_eta[Idx(i, nj - 3, k)].z_xi_x = MidNodeInter4thLeft2(inter_temp);
                temp_eta[Idx(i, nj - 4, k)].z_xi_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, nj - 5 + iTemp, k))[0] * node->GetCoord(i, nj - 5 + iTemp, k)[1];
                }
                temp_eta[Idx(i, nj - 2, k)].x_xi_y = MidNodeInter4thLeft1(inter_temp);
                temp_eta[Idx(i, nj - 3, k)].x_xi_y = MidNodeInter4thLeft2(inter_temp);
                temp_eta[Idx(i, nj - 4, k)].x_xi_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, nj - 5 + iTemp, k))[2] * node->GetCoord(i, nj - 5 + iTemp, k)[1];
                }
                temp_eta[Idx(i, nj - 2, k)].z_zeta_y = MidNodeInter4thLeft1(inter_temp);
                temp_eta[Idx(i, nj - 3, k)].z_zeta_y = MidNodeInter4thLeft2(inter_temp);
                temp_eta[Idx(i, nj - 4, k)].z_zeta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, nj - 5 + iTemp, k))[2] * node->GetCoord(i, nj - 5 + iTemp, k)[2];
                }
                temp_eta[Idx(i, nj - 2, k)].x_zeta_z = MidNodeInter4thLeft1(inter_temp);
                temp_eta[Idx(i, nj - 3, k)].x_zeta_z = MidNodeInter4thLeft2(inter_temp);
                temp_eta[Idx(i, nj - 4, k)].x_zeta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, nj - 5 + iTemp, k))[2] * node->GetCoord(i, nj - 5 + iTemp, k)[0];
                }
                temp_eta[Idx(i, nj - 2, k)].y_zeta_x = MidNodeInter4thLeft1(inter_temp);
                temp_eta[Idx(i, nj - 3, k)].y_zeta_x = MidNodeInter4thLeft2(inter_temp);
                temp_eta[Idx(i, nj - 4, k)].y_zeta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, nj - 5 + iTemp, k))[0] * node->GetCoord(i, nj - 5 + iTemp, k)[1];
                }
                temp_eta[Idx(i, nj - 2, k)].z_xi_y = MidNodeInter4thLeft1(inter_temp);
                temp_eta[Idx(i, nj - 3, k)].z_xi_y = MidNodeInter4thLeft2(inter_temp);
                temp_eta[Idx(i, nj - 4, k)].z_xi_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, nj - 5 + iTemp, k))[0] * node->GetCoord(i, nj - 5 + iTemp, k)[2];
                }
                temp_eta[Idx(i, nj - 2, k)].x_xi_z = MidNodeInter4thLeft1(inter_temp);
                temp_eta[Idx(i, nj - 3, k)].x_xi_z = MidNodeInter4thLeft2(inter_temp);
                temp_eta[Idx(i, nj - 4, k)].x_xi_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, nj - 5 + iTemp, k))[0] * node->GetCoord(i, nj - 5 + iTemp, k)[0];
                }
                temp_eta[Idx(i, nj - 2, k)].y_xi_x = MidNodeInter4thLeft1(inter_temp);
                temp_eta[Idx(i, nj - 3, k)].y_xi_x = MidNodeInter4thLeft2(inter_temp);
                temp_eta[Idx(i, nj - 4, k)].y_xi_x = MidNodeInter4th(inter_temp);
            }
        }

        // 求出zeta方向的临时变量：逆变换度量系数乘以坐标
        struct TempZeta
        {
            double y_eta_z, z_eta_x, x_eta_y, y_xi_z, z_xi_x, x_xi_y, z_eta_y, x_eta_z, y_eta_x, z_xi_y, x_xi_z, y_xi_x;
        };
        std::vector<TempZeta> temp_zeta(ni * nj * nk);
        for (int k = 3; k < nk - 4; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetY(Idx(i, j, k + iTemp - 2))[1] * node->GetCoord(i, j, k + iTemp - 2)[2];
                    }
                    temp_zeta[Idx(i, j, k)].y_eta_z = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZ(Idx(i, j, k + iTemp - 2))[1] * node->GetCoord(i, j, k + iTemp - 2)[0];
                    }
                    temp_zeta[Idx(i, j, k)].z_eta_x = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetX(Idx(i, j, k + iTemp - 2))[1] * node->GetCoord(i, j, k + iTemp - 2)[1];
                    }
                    temp_zeta[Idx(i, j, k)].x_eta_y = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetY(Idx(i, j, k + iTemp - 2))[0] * node->GetCoord(i, j, k + iTemp - 2)[2];
                    }
                    temp_zeta[Idx(i, j, k)].y_xi_z = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZ(Idx(i, j, k + iTemp - 2))[0] * node->GetCoord(i, j, k + iTemp - 2)[0];
                    }
                    temp_zeta[Idx(i, j, k)].z_xi_x = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetX(Idx(i, j, k + iTemp - 2))[0] * node->GetCoord(i, j, k + iTemp - 2)[1];
                    }
                    temp_zeta[Idx(i, j, k)].x_xi_y = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZ(Idx(i, j, k - 2 + iTemp))[1] * node->GetCoord(i, j, k - 2 + iTemp)[1];
                    }
                    temp_zeta[Idx(i, j, k)].z_eta_y = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetX(Idx(i, j, k - 2 + iTemp))[1] * node->GetCoord(i, j, k - 2 + iTemp)[2];
                    }
                    temp_zeta[Idx(i, j, k)].x_eta_z = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetY(Idx(i, j, k - 2 + iTemp))[1] * node->GetCoord(i, j, k - 2 + iTemp)[0];
                    }
                    temp_zeta[Idx(i, j, k)].y_eta_x = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetZ(Idx(i, j, k - 2 + iTemp))[0] * node->GetCoord(i, j, k - 2 + iTemp)[1];
                    }
                    temp_zeta[Idx(i, j, k)].z_xi_y = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetX(Idx(i, j, k - 2 + iTemp))[0] * node->GetCoord(i, j, k - 2 + iTemp)[2];
                    }
                    temp_zeta[Idx(i, j, k)].x_xi_z = MidNodeInter6th(inter_temp);
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        inter_temp[iTemp] = coef->GetY(Idx(i, j, k - 2 + iTemp))[0] * node->GetCoord(i, j, k - 2 + iTemp)[0];
                    }
                    temp_zeta[Idx(i, j, k)].y_xi_x = MidNodeInter6th(inter_temp);
                }
            }
        }
        // k=1/2,3/2,5/2;nk-3/2,nk-5/2,nk-7/2
        for (int j = 1; j < nj - 1; ++j)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                // k=1/2，3/2，5/2
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, j, 1 + iTemp))[1] * node->GetCoord(i, j, 1 + iTemp)[2];
                }
                temp_zeta[Idx(i, j, 0)].y_eta_z = MidNodeInter4thRight1(inter_temp);
                temp_zeta[Idx(i, j, 1)].y_eta_z = MidNodeInter4thRight2(inter_temp);
                temp_zeta[Idx(i, j, 2)].y_eta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, j, 1 + iTemp))[1] * node->GetCoord(i, j, 1 + iTemp)[0];
                }
                temp_zeta[Idx(i, j, 0)].z_eta_x = MidNodeInter4thRight1(inter_temp);
                temp_zeta[Idx(i, j, 1)].z_eta_x = MidNodeInter4thRight2(inter_temp);
                temp_zeta[Idx(i, j, 2)].z_eta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, j, 1 + iTemp))[1] * node->GetCoord(i, j, 1 + iTemp)[1];
                }
                temp_zeta[Idx(i, j, 0)].x_eta_y = MidNodeInter4thRight1(inter_temp);
                temp_zeta[Idx(i, j, 1)].x_eta_y = MidNodeInter4thRight2(inter_temp);
                temp_zeta[Idx(i, j, 2)].x_eta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, j, 1 + iTemp))[0] * node->GetCoord(i, j, 1 + iTemp)[2];
                }
                temp_zeta[Idx(i, j, 0)].y_xi_z = MidNodeInter4thRight1(inter_temp);
                temp_zeta[Idx(i, j, 1)].y_xi_z = MidNodeInter4thRight2(inter_temp);
                temp_zeta[Idx(i, j, 2)].y_xi_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, j, 1 + iTemp))[0] * node->GetCoord(i, j, 1 + iTemp)[0];
                }
                temp_zeta[Idx(i, j, 0)].z_xi_x = MidNodeInter4thRight1(inter_temp);
                temp_zeta[Idx(i, j, 1)].z_xi_x = MidNodeInter4thRight2(inter_temp);
                temp_zeta[Idx(i, j, 2)].z_xi_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, j, 1 + iTemp))[0] * node->GetCoord(i, j, 1 + iTemp)[1];
                }
                temp_zeta[Idx(i, j, 0)].x_xi_y = MidNodeInter4thRight1(inter_temp);
                temp_zeta[Idx(i, j, 1)].x_xi_y = MidNodeInter4thRight2(inter_temp);
                temp_zeta[Idx(i, j, 2)].x_xi_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, j, 1 + iTemp))[1] * node->GetCoord(i, j, 1 + iTemp)[1];
                }
                temp_zeta[Idx(i, j, 0)].z_eta_y = MidNodeInter4thRight1(inter_temp);
                temp_zeta[Idx(i, j, 1)].z_eta_y = MidNodeInter4thRight2(inter_temp);
                temp_zeta[Idx(i, j, 2)].z_eta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, j, 1 + iTemp))[1] * node->GetCoord(i, j, 1 + iTemp)[2];
                }
                temp_zeta[Idx(i, j, 0)].x_eta_z = MidNodeInter4thRight1(inter_temp);
                temp_zeta[Idx(i, j, 1)].x_eta_z = MidNodeInter4thRight2(inter_temp);
                temp_zeta[Idx(i, j, 2)].x_eta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, j, 1 + iTemp))[1] * node->GetCoord(i, j, 1 + iTemp)[0];
                }
                temp_zeta[Idx(i, j, 0)].y_eta_x = MidNodeInter4thRight1(inter_temp);
                temp_zeta[Idx(i, j, 1)].y_eta_x = MidNodeInter4thRight2(inter_temp);
                temp_zeta[Idx(i, j, 2)].y_eta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, j, 1 + iTemp))[0] * node->GetCoord(i, j, 1 + iTemp)[1];
                }
                temp_zeta[Idx(i, j, 0)].z_xi_y = MidNodeInter4thRight1(inter_temp);
                temp_zeta[Idx(i, j, 1)].z_xi_y = MidNodeInter4thRight2(inter_temp);
                temp_zeta[Idx(i, j, 2)].z_xi_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, j, 1 + iTemp))[0] * node->GetCoord(i, j, 1 + iTemp)[2];
                }
                temp_zeta[Idx(i, j, 0)].x_xi_z = MidNodeInter4thRight1(inter_temp);
                temp_zeta[Idx(i, j, 1)].x_xi_z = MidNodeInter4thRight2(inter_temp);
                temp_zeta[Idx(i, j, 2)].x_xi_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, j, 1 + iTemp))[0] * node->GetCoord(i, j, 1 + iTemp)[0];
                }
                temp_zeta[Idx(i, j, 0)].y_xi_x = MidNodeInter4thRight1(inter_temp);
                temp_zeta[Idx(i, j, 1)].y_xi_x = MidNodeInter4thRight2(inter_temp);
                temp_zeta[Idx(i, j, 2)].y_xi_x = MidNodeInter4th(inter_temp);
                // k=nk-3/2,nk-5/2,nk-7/2
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, j, nk - 5 + iTemp))[1] * node->GetCoord(i, j, nk - 5 + iTemp)[2];
                }
                temp_zeta[Idx(i, j, nk - 2)].y_eta_z = MidNodeInter4thLeft1(inter_temp);
                temp_zeta[Idx(i, j, nk - 3)].y_eta_z = MidNodeInter4thLeft2(inter_temp);
                temp_zeta[Idx(i, j, nk - 4)].y_eta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, j, nk - 5 + iTemp))[1] * node->GetCoord(i, j, nk - 5 + iTemp)[0];
                }
                temp_zeta[Idx(i, j, nk - 2)].z_eta_x = MidNodeInter4thLeft1(inter_temp);
                temp_zeta[Idx(i, j, nk - 3)].z_eta_x = MidNodeInter4thLeft2(inter_temp);
                temp_zeta[Idx(i, j, nk - 4)].z_eta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, j, nk - 5 + iTemp))[1] * node->GetCoord(i, j, nk - 5 + iTemp)[1];
                }
                temp_zeta[Idx(i, j, nk - 2)].x_eta_y = MidNodeInter4thLeft1(inter_temp);
                temp_zeta[Idx(i, j, nk - 3)].x_eta_y = MidNodeInter4thLeft2(inter_temp);
                temp_zeta[Idx(i, j, nk - 4)].x_eta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, j, nk - 5 + iTemp))[0] * node->GetCoord(i, j, nk - 5 + iTemp)[2];
                }
                temp_zeta[Idx(i, j, nk - 2)].y_xi_z = MidNodeInter4thLeft1(inter_temp);
                temp_zeta[Idx(i, j, nk - 3)].y_xi_z = MidNodeInter4thLeft2(inter_temp);
                temp_zeta[Idx(i, j, nk - 4)].y_xi_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, j, nk - 5 + iTemp))[0] * node->GetCoord(i, j, nk - 5 + iTemp)[0];
                }
                temp_zeta[Idx(i, j, nk - 2)].z_xi_x = MidNodeInter4thLeft1(inter_temp);
                temp_zeta[Idx(i, j, nk - 3)].z_xi_x = MidNodeInter4thLeft2(inter_temp);
                temp_zeta[Idx(i, j, nk - 4)].z_xi_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, j, nk - 5 + iTemp))[0] * node->GetCoord(i, j, nk - 5 + iTemp)[1];
                }
                temp_zeta[Idx(i, j, nk - 2)].x_xi_y = MidNodeInter4thLeft1(inter_temp);
                temp_zeta[Idx(i, j, nk - 3)].x_xi_y = MidNodeInter4thLeft2(inter_temp);
                temp_zeta[Idx(i, j, nk - 4)].x_xi_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, j, nk - 5 + iTemp))[1] * node->GetCoord(i, j, nk - 5 + iTemp)[1];
                }
                temp_zeta[Idx(i, j, nk - 2)].z_eta_y = MidNodeInter4thLeft1(inter_temp);
                temp_zeta[Idx(i, j, nk - 3)].z_eta_y = MidNodeInter4thLeft2(inter_temp);
                temp_zeta[Idx(i, j, nk - 4)].z_eta_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, j, nk - 5 + iTemp))[1] * node->GetCoord(i, j, nk - 5 + iTemp)[2];
                }
                temp_zeta[Idx(i, j, nk - 2)].x_eta_z = MidNodeInter4thLeft1(inter_temp);
                temp_zeta[Idx(i, j, nk - 3)].x_eta_z = MidNodeInter4thLeft2(inter_temp);
                temp_zeta[Idx(i, j, nk - 4)].x_eta_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, j, nk - 5 + iTemp))[1] * node->GetCoord(i, j, nk - 5 + iTemp)[0];
                }
                temp_zeta[Idx(i, j, nk - 2)].y_eta_x = MidNodeInter4thLeft1(inter_temp);
                temp_zeta[Idx(i, j, nk - 3)].y_eta_x = MidNodeInter4thLeft2(inter_temp);
                temp_zeta[Idx(i, j, nk - 4)].y_eta_x = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetZ(Idx(i, j, nk - 5 + iTemp))[0] * node->GetCoord(i, j, nk - 5 + iTemp)[1];
                }
                temp_zeta[Idx(i, j, nk - 2)].z_xi_y = MidNodeInter4thLeft1(inter_temp);
                temp_zeta[Idx(i, j, nk - 3)].z_xi_y = MidNodeInter4thLeft2(inter_temp);
                temp_zeta[Idx(i, j, nk - 4)].z_xi_y = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetX(Idx(i, j, nk - 5 + iTemp))[0] * node->GetCoord(i, j, nk - 5 + iTemp)[2];
                }
                temp_zeta[Idx(i, j, nk - 2)].x_xi_z = MidNodeInter4thLeft1(inter_temp);
                temp_zeta[Idx(i, j, nk - 3)].x_xi_z = MidNodeInter4thLeft2(inter_temp);
                temp_zeta[Idx(i, j, nk - 4)].x_xi_z = MidNodeInter4th(inter_temp);
                for (int iTemp = 0; iTemp < 4; iTemp++)
                {
                    inter_temp[iTemp] = coef->GetY(Idx(i, j, nk - 5 + iTemp))[0] * node->GetCoord(i, j, nk - 5 + iTemp)[0];
                }
                temp_zeta[Idx(i, j, nk - 2)].y_xi_x = MidNodeInter4thLeft1(inter_temp);
                temp_zeta[Idx(i, j, nk - 3)].y_xi_x = MidNodeInter4thLeft2(inter_temp);
                temp_zeta[Idx(i, j, nk - 4)].y_xi_x = MidNodeInter4th(inter_temp);
            }
        }
        // 第五步：根据半点坐标和半点逆变换度量系数使用守恒形式计算整点度量系数（CMM1)
        // i direction
        double temp[6][6];
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 3; i < ni - 3; ++i)
                {
                    auto coef_xi = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, k));
                    auto coef_eta = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, k));
                    auto coef_zeta = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, k));
                    // i direction
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        int idx = Idx(i - 3 + iTemp, j, k);
                        temp[0][iTemp] = 0.5 * (temp_xi[idx].y_zeta_z - temp_xi[idx].z_zeta_y); // y_zeta*z-z_zeta*y
                        temp[1][iTemp] = 0.5 * (temp_xi[idx].z_zeta_x - temp_xi[idx].x_zeta_z); // z_zeta*x-x_zeta*z
                        temp[2][iTemp] = 0.5 * (temp_xi[idx].x_zeta_y - temp_xi[idx].y_zeta_x); // x_zeta*y-y_zeta*x
                        temp[3][iTemp] = 0.5 * (temp_xi[idx].z_eta_y - temp_xi[idx].y_eta_z);   // z_eta*y-y_eta*z
                        temp[4][iTemp] = 0.5 * (temp_xi[idx].x_eta_z - temp_xi[idx].z_eta_x);   // x_eta*z-z_eta*x
                        temp[5][iTemp] = 0.5 * (temp_xi[idx].y_eta_x - temp_xi[idx].x_eta_y);   // y_eta*x-x_eta*y
                    }
                    coef_eta[0] += NodeDifferece6th(temp[0]);
                    coef_eta[1] += NodeDifferece6th(temp[1]);
                    coef_eta[2] += NodeDifferece6th(temp[2]);
                    coef_zeta[0] += NodeDifferece6th(temp[3]);
                    coef_zeta[1] += NodeDifferece6th(temp[4]);
                    coef_zeta[2] += NodeDifferece6th(temp[5]);
                }
            }
        }
        // j direction
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 3; j < nj - 3; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    auto coef_xi = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, k));
                    auto coef_eta = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, k));
                    auto coef_zeta = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        int idx = Idx(i, j - 3 + iTemp, k);
                        temp[0][iTemp] = 0.5 * (temp_eta[idx].z_zeta_y - temp_eta[idx].y_zeta_z); // z_zeta*y-y_zeta*z
                        temp[1][iTemp] = 0.5 * (temp_eta[idx].x_zeta_z - temp_eta[idx].z_zeta_x); // x_zeta*z-z_zeta*x
                        temp[2][iTemp] = 0.5 * (temp_eta[idx].y_zeta_x - temp_eta[idx].x_zeta_y); //  y_zeta*x-x_zeta*y
                        temp[3][iTemp] = 0.5 * (temp_eta[idx].y_xi_z - temp_eta[idx].z_xi_y);     // y_xi*z-z_xi*y
                        temp[4][iTemp] = 0.5 * (temp_eta[idx].z_xi_x - temp_eta[idx].x_xi_z);     // z_xi*x-x_xi*z
                        temp[5][iTemp] = 0.5 * (temp_eta[idx].x_xi_y - temp_eta[idx].y_xi_x);     // x_xi*y-y_xi*x
                    }
                    coef_xi[0] += NodeDifferece6th(temp[0]);
                    coef_xi[1] += NodeDifferece6th(temp[1]);
                    coef_xi[2] += NodeDifferece6th(temp[2]);
                    coef_zeta[0] += NodeDifferece6th(temp[3]);
                    coef_zeta[1] += NodeDifferece6th(temp[4]);
                    coef_zeta[2] += NodeDifferece6th(temp[5]);
                }
            }
        }
        // k direction
        for (int k = 3; k < nk - 3; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    auto coef_xi = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, k));
                    auto coef_eta = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, k));
                    auto coef_zeta = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, k));
                    for (int iTemp = 0; iTemp < 6; iTemp++)
                    {
                        int idx = Idx(i, j, k - 3 + iTemp);
                        temp[0][iTemp] = 0.5 * (temp_zeta[idx].y_eta_z - temp_zeta[idx].z_eta_y); // y_eta*z-z_eta*y
                        temp[1][iTemp] = 0.5 * (temp_zeta[idx].z_eta_x - temp_zeta[idx].x_eta_z); // z_eta*x-x_eta*z
                        temp[2][iTemp] = 0.5 * (temp_zeta[idx].x_eta_y - temp_zeta[idx].y_eta_x); // x_eta*y-y_eta*x
                        temp[3][iTemp] = 0.5 * (temp_zeta[idx].z_xi_y - temp_zeta[idx].y_xi_z);   // z_xi*y-y_xi*z
                        temp[4][iTemp] = 0.5 * (temp_zeta[idx].x_xi_z - temp_zeta[idx].z_xi_x);   // x_xi*z-z_xi*x
                        temp[5][iTemp] = 0.5 * (temp_zeta[idx].y_xi_x - temp_zeta[idx].x_xi_y);   //  y_xi*x-x_xi*y
                    }
                    coef_xi[0] += NodeDifferece6th(temp[0]);
                    coef_xi[1] += NodeDifferece6th(temp[1]);
                    coef_xi[2] += NodeDifferece6th(temp[2]);
                    coef_eta[0] += NodeDifferece6th(temp[3]);
                    coef_eta[1] += NodeDifferece6th(temp[4]);
                    coef_eta[2] += NodeDifferece6th(temp[5]);
                }
            }
        }
        // 第六步：计算边界点度量系数
        // i=1,2;ni-3,ni-2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                auto coef_xi1 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(1, j, k));
                auto coef_eta1 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(1, j, k));
                auto coef_zeta1 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(1, j, k));
                auto coef_xi2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(2, j, k));
                auto coef_eta2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(2, j, k));
                auto coef_zeta2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(2, j, k));
                for (int iTemp = 0; iTemp < 5; iTemp++)
                {
                    int idx = Idx(0 + iTemp, j, k);
                    temp[0][iTemp] = 0.5 * (temp_xi[idx].y_zeta_z - temp_xi[idx].z_zeta_y); // y_zeta*z-z_zeta*y
                    temp[1][iTemp] = 0.5 * (temp_xi[idx].z_zeta_x - temp_xi[idx].x_zeta_z); // z_zeta*x-x_zeta*z
                    temp[2][iTemp] = 0.5 * (temp_xi[idx].x_zeta_y - temp_xi[idx].y_zeta_x); // x_zeta*y-y_zeta*x
                    temp[3][iTemp] = 0.5 * (temp_xi[idx].z_eta_y - temp_xi[idx].y_eta_z);   // z_eta*y-y_eta*z
                    temp[4][iTemp] = 0.5 * (temp_xi[idx].x_eta_z - temp_xi[idx].z_eta_x);   // x_eta*z-z_eta*x
                    temp[5][iTemp] = 0.5 * (temp_xi[idx].y_eta_x - temp_xi[idx].x_eta_y);   // y_eta*x-x_eta*y
                }
                coef_eta1[0] += NodeDifferece4thRight(temp[0]);
                coef_eta1[1] += NodeDifferece4thRight(temp[1]);
                coef_eta1[2] += NodeDifferece4thRight(temp[2]);
                coef_zeta1[0] += NodeDifferece4thRight(temp[3]);
                coef_zeta1[1] += NodeDifferece4thRight(temp[4]);
                coef_zeta1[2] += NodeDifferece4thRight(temp[5]);
                coef_eta2[0] += NodeDifferece4th(temp[0]);
                coef_eta2[1] += NodeDifferece4th(temp[1]);
                coef_eta2[2] += NodeDifferece4th(temp[2]);
                coef_zeta2[0] += NodeDifferece4th(temp[3]);
                coef_zeta2[1] += NodeDifferece4th(temp[4]);
                coef_zeta2[2] = NodeDifferece4th(temp[5]);
                auto coef_xi_n2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(ni - 2, j, k));
                auto coef_eta_n2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(ni - 2, j, k));
                auto coef_zeta_n2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(ni - 2, j, k));
                auto coef_xi_n3 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(ni - 3, j, k));
                auto coef_eta_n3 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(ni - 3, j, k));
                auto coef_zeta_n3 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(ni - 3, j, k));
                for (int iTemp = 0; iTemp < 5; iTemp++)
                {
                    int idx = Idx(ni - 2 - iTemp, j, k);
                    temp[0][iTemp] = 0.5 * (temp_xi[idx].y_zeta_z - temp_xi[idx].z_zeta_y); // y_zeta*z-z_zeta*y
                    temp[1][iTemp] = 0.5 * (temp_xi[idx].z_zeta_x - temp_xi[idx].x_zeta_z); // z_zeta*x-x_zeta*z
                    temp[2][iTemp] = 0.5 * (temp_xi[idx].x_zeta_y - temp_xi[idx].y_zeta_x); // x_zeta*y-y_zeta*x
                    temp[3][iTemp] = 0.5 * (temp_xi[idx].z_eta_y - temp_xi[idx].y_eta_z);   // z_eta*y-y_eta*z
                    temp[4][iTemp] = 0.5 * (temp_xi[idx].x_eta_z - temp_xi[idx].z_eta_x);   // x_eta*z-z_eta*x
                    temp[5][iTemp] = 0.5 * (temp_xi[idx].y_eta_x - temp_xi[idx].x_eta_y);   // y_eta*x-x_eta*y
                }
                coef_eta_n2[0] += NodeDifferece4thLeft(temp[0]);
                coef_eta_n2[1] += NodeDifferece4thLeft(temp[1]);
                coef_eta_n2[2] += NodeDifferece4thLeft(temp[2]);
                coef_zeta_n2[0] += NodeDifferece4thLeft(temp[3]);
                coef_zeta_n2[1] += NodeDifferece4thLeft(temp[4]);
                coef_zeta_n2[2] += NodeDifferece4thLeft(temp[5]);
                coef_eta_n3[0] += -NodeDifferece4th(temp[0]);
                coef_eta_n3[1] += -NodeDifferece4th(temp[1]);
                coef_eta_n3[2] += -NodeDifferece4th(temp[2]);
                coef_zeta_n3[0] += -NodeDifferece4th(temp[3]);
                coef_zeta_n3[1] += -NodeDifferece4th(temp[4]);
                coef_zeta_n3[2] += -NodeDifferece4th(temp[5]);
            }
        }
        // j=1,2;nj-3,nj-2
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                auto coef_xi1 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, 1, k));
                auto coef_eta1 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, 1, k));
                auto coef_zeta1 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, 1, k));
                auto coef_xi2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, 2, k));
                auto coef_eta2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, 2, k));
                auto coef_zeta2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, 2, k));
                for (int iTemp = 0; iTemp < 5; iTemp++)
                {
                    int idx = Idx(i, 0 + iTemp, k);
                    temp[0][iTemp] = 0.5 * (temp_eta[idx].z_zeta_y - temp_eta[idx].y_zeta_z); // z_zeta*y-y_zeta*z
                    temp[1][iTemp] = 0.5 * (temp_eta[idx].x_zeta_z - temp_eta[idx].z_zeta_x); // x_zeta*z-z_zeta*x
                    temp[2][iTemp] = 0.5 * (temp_eta[idx].y_zeta_x - temp_eta[idx].x_zeta_y); //  y_zeta*x-x_zeta*y
                    temp[3][iTemp] = 0.5 * (temp_eta[idx].y_xi_z - temp_eta[idx].z_xi_y);     // y_xi*z-z_xi*y
                    temp[4][iTemp] = 0.5 * (temp_eta[idx].z_xi_x - temp_eta[idx].x_xi_z);     // z_xi*x-x_xi*z
                    temp[5][iTemp] = 0.5 * (temp_eta[idx].x_xi_y - temp_eta[idx].y_xi_x);     // x_xi*y-y_xi*x
                }
                coef_xi1[0] += NodeDifferece4thRight(temp[0]);
                coef_xi1[1] += NodeDifferece4thRight(temp[1]);
                coef_xi1[2] += NodeDifferece4thRight(temp[2]);
                coef_zeta1[0] += NodeDifferece4thRight(temp[3]);
                coef_zeta1[1] += NodeDifferece4thRight(temp[4]);
                coef_zeta1[2] += NodeDifferece4thRight(temp[5]);
                coef_xi2[0] += NodeDifferece4th(temp[0]);
                coef_xi2[1] += NodeDifferece4th(temp[1]);
                coef_xi2[2] += NodeDifferece4th(temp[2]);
                coef_zeta2[0] += NodeDifferece4th(temp[3]);
                coef_zeta2[1] += NodeDifferece4th(temp[4]);
                coef_zeta2[2] += NodeDifferece4th(temp[5]);
                auto coef_xi_n2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, nj - 2, k));
                auto coef_eta_n2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, nj - 2, k));
                auto coef_zeta_n2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, nj - 2, k));
                auto coef_xi_n3 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, nj - 3, k));
                auto coef_eta_n3 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, nj - 3, k));
                auto coef_zeta_n3 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, nj - 3, k));
                for (int iTemp = 0; iTemp < 5; iTemp++)
                {
                    int idx = Idx(i, nj - 2 - iTemp, k);
                    temp[0][iTemp] = 0.5 * (temp_eta[idx].z_zeta_y - temp_eta[idx].y_zeta_z); // z_zeta*y-y_zeta*z
                    temp[1][iTemp] = 0.5 * (temp_eta[idx].x_zeta_z - temp_eta[idx].z_zeta_x); // x_zeta*z-z_zeta*x
                    temp[2][iTemp] = 0.5 * (temp_eta[idx].y_zeta_x - temp_eta[idx].x_zeta_y); //  y_zeta*x-x_zeta*y
                    temp[3][iTemp] = 0.5 * (temp_eta[idx].y_xi_z - temp_eta[idx].z_xi_y);     // y_xi*z-z_xi*y
                    temp[4][iTemp] = 0.5 * (temp_eta[idx].z_xi_x - temp_eta[idx].x_xi_z);     // z_xi*x-x_xi*z
                    temp[5][iTemp] = 0.5 * (temp_eta[idx].x_xi_y - temp_eta[idx].y_xi_x);     // x_xi*y-y_xi*x
                }
                coef_xi_n2[0] += NodeDifferece4thLeft(temp[0]);
                coef_xi_n2[1] += NodeDifferece4thLeft(temp[1]);
                coef_xi_n2[2] += NodeDifferece4thLeft(temp[2]);
                coef_zeta_n2[0] += NodeDifferece4thLeft(temp[3]);
                coef_zeta_n2[1] += NodeDifferece4thLeft(temp[4]);
                coef_zeta_n2[2] += NodeDifferece4thLeft(temp[5]);
                coef_xi_n3[0] += -NodeDifferece4th(temp[0]);
                coef_xi_n3[1] += -NodeDifferece4th(temp[1]);
                coef_xi_n3[2] += -NodeDifferece4th(temp[2]);
                coef_zeta_n3[0] += -NodeDifferece4th(temp[3]);
                coef_zeta_n3[1] += -NodeDifferece4th(temp[4]);
                coef_zeta_n3[2] += -NodeDifferece4th(temp[5]);
            }
        }
        // k=1,2;nk-3,nk-2
        for (int j = 1; j < nj - 1; ++j)
        {
            for (int i = 1; i < ni - 1; ++i)
            {
                auto coef_xi1 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, 1));
                auto coef_eta1 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, 1));
                auto coef_zeta1 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, 1));
                auto coef_xi2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, 2));
                auto coef_eta2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, 2));
                auto coef_zeta2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, 2));
                for (int iTemp = 0; iTemp < 5; iTemp++)
                {
                    int idx = Idx(i, j, 0 + iTemp);
                    temp[0][iTemp] = 0.5 * (temp_zeta[idx].y_eta_z - temp_zeta[idx].z_eta_y); // y_eta*z-z_eta*y
                    temp[1][iTemp] = 0.5 * (temp_zeta[idx].z_eta_x - temp_zeta[idx].x_eta_z); // z_eta*x-x_eta*z
                    temp[2][iTemp] = 0.5 * (temp_zeta[idx].x_eta_y - temp_zeta[idx].y_eta_x); // x_eta*y-y_eta*x
                    temp[3][iTemp] = 0.5 * (temp_zeta[idx].z_xi_y - temp_zeta[idx].y_xi_z);   // z_xi*y-y_xi*z
                    temp[4][iTemp] = 0.5 * (temp_zeta[idx].x_xi_z - temp_zeta[idx].z_xi_x);   // x_xi*z-z_xi*x
                    temp[5][iTemp] = 0.5 * (temp_zeta[idx].y_xi_x - temp_zeta[idx].x_xi_y);   //  y_xi*x-x_xi*y
                }
                coef_xi1[0] += NodeDifferece4thRight(temp[0]);
                coef_xi1[1] += NodeDifferece4thRight(temp[1]);
                coef_xi1[2] += NodeDifferece4thRight(temp[2]);
                coef_eta1[0] += NodeDifferece4thRight(temp[3]);
                coef_eta1[1] += NodeDifferece4thRight(temp[4]);
                coef_eta1[2] += NodeDifferece4thRight(temp[5]);
                coef_xi2[0] += NodeDifferece4th(temp[0]);
                coef_xi2[1] += NodeDifferece4th(temp[1]);
                coef_xi2[2] += NodeDifferece4th(temp[2]);
                coef_eta2[0] += NodeDifferece4th(temp[3]);
                coef_eta2[1] += NodeDifferece4th(temp[4]);
                coef_eta2[2] += NodeDifferece4th(temp[5]);
                auto coef_xi_n2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, nk - 2));
                auto coef_eta_n2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, nk - 2));
                auto coef_zeta_n2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, nk - 2));
                auto coef_xi_n3 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, nk - 3));
                auto coef_eta_n3 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, nk - 3));
                auto coef_zeta_n3 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, nk - 3));
                for (int iTemp = 0; iTemp < 5; iTemp++)
                {
                    int idx = Idx(i, j, nk - 2 - iTemp);
                    temp[0][iTemp] = 0.5 * (temp_zeta[idx].y_eta_z - temp_zeta[idx].z_eta_y); // y_eta*z-z_eta*y
                    temp[1][iTemp] = 0.5 * (temp_zeta[idx].z_eta_x - temp_zeta[idx].x_eta_z); // z_eta*x-x_eta*z
                    temp[2][iTemp] = 0.5 * (temp_zeta[idx].x_eta_y - temp_zeta[idx].y_eta_x); // x_eta*y-y_eta*x
                    temp[3][iTemp] = 0.5 * (temp_zeta[idx].z_xi_y - temp_zeta[idx].y_xi_z);   // z_xi*y-y_xi*z
                    temp[4][iTemp] = 0.5 * (temp_zeta[idx].x_xi_z - temp_zeta[idx].z_xi_x);   // x_xi*z-z_xi*x
                    temp[5][iTemp] = 0.5 * (temp_zeta[idx].y_xi_x - temp_zeta[idx].x_xi_y);   //  y_xi*x-x_xi*y
                }
                coef_xi_n2[0] += NodeDifferece4thLeft(temp[0]);
                coef_xi_n2[1] += NodeDifferece4thLeft(temp[1]);
                coef_xi_n2[2] += NodeDifferece4thLeft(temp[2]);
                coef_eta_n2[0] += NodeDifferece4thLeft(temp[3]);
                coef_eta_n2[1] += NodeDifferece4thLeft(temp[4]);
                coef_eta_n2[2] += NodeDifferece4thLeft(temp[5]);
                coef_xi_n3[0] += -NodeDifferece4th(temp[0]);
                coef_xi_n3[1] += -NodeDifferece4th(temp[1]);
                coef_xi_n3[2] += -NodeDifferece4th(temp[2]);
                coef_eta_n3[0] += -NodeDifferece4th(temp[3]);
                coef_eta_n3[1] += -NodeDifferece4th(temp[4]);
                coef_eta_n3[2] += -NodeDifferece4th(temp[5]);
            }
        }
    }
    void NSSolverStruct::CalcMetricsS4()
    {
        // TODO
    }
    void NSSolverStruct::Preprocess()
    {
        NSSolver::Preprocess();
        auto para = GetPara();
        if (para->GetInterSchme() == InterpolationScheme::Grad)
        {
            CalcPrimGradWLS();
            CalcLimiterVK();
        }
    }
    void NSSolverStruct::CalcPrimGradWLS()
    {
        auto grid = GetGrid();
        auto node = grid->GetNode();
        auto idx_proxy = GetIdxProxy();
        auto data_manager = GetDataManager();
        int ni, nj, nk;
        ni = idx_proxy->GetNi();
        nj = idx_proxy->GetNj();
        nk = idx_proxy->GetNk();
        Matrix3d A, A_inv;
        DVector3D b, grad;
        int temp_i[6], temp_j[6], temp_k[6];
        double delta_x, delta_y, delta_z, delta_val;
        double weight;
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    A.setZero();
                    int idx = idx_proxy->GetIdx(i, j, k);
                    auto coord = node->GetCoord(i, j, k);
                    temp_i[0] = i - 1, temp_j[0] = j, temp_k[0] = k;
                    temp_i[1] = i + 1, temp_j[1] = j, temp_k[1] = k;
                    temp_i[2] = i, temp_j[2] = j - 1, temp_k[2] = k;
                    temp_i[3] = i, temp_j[3] = j + 1, temp_k[3] = k;
                    temp_i[4] = i, temp_j[4] = j, temp_k[4] = k - 1;
                    temp_i[5] = i, temp_j[5] = j, temp_k[5] = k + 1;
                    for (int iNeigh = 0; iNeigh < 6; ++iNeigh)
                    {
                        auto coord_neigh = node->GetCoord(temp_i[iNeigh], temp_j[iNeigh], temp_k[iNeigh]);
                        delta_x = coord_neigh[0] - coord[0];
                        delta_y = coord_neigh[1] - coord[1];
                        delta_z = coord_neigh[2] - coord[2];
                        weight = 1.0 / sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
                        A(0, 0) += weight * delta_x * delta_x;
                        A(0, 1) += weight * delta_x * delta_y;
                        A(0, 2) += weight * delta_x * delta_z;
                        A(1, 0) += weight * delta_y * delta_x;
                        A(1, 1) += weight * delta_y * delta_y;
                        A(1, 2) += weight * delta_y * delta_z;
                        A(2, 0) += weight * delta_z * delta_x;
                        A(2, 1) += weight * delta_z * delta_y;
                        A(2, 2) += weight * delta_z * delta_z;
                    }
                    A_inv = A.inverse();
                    for (int idx_eq = 0; idx_eq < 5; ++idx_eq)
                    {
                        b.setZero();
                        for (int iNeigh = 0; iNeigh < 6; ++iNeigh)
                        {
                            delta_val = data_manager->GetPrim(idx_eq, idx_proxy->GetIdx(temp_i[iNeigh], temp_j[iNeigh], temp_k[iNeigh])) -
                                        data_manager->GetPrim(idx_eq, idx);
                            auto coord_neigh = node->GetCoord(temp_i[iNeigh], temp_j[iNeigh], temp_k[iNeigh]);
                            delta_x = coord_neigh[0] - coord[0];
                            delta_y = coord_neigh[1] - coord[1];
                            delta_z = coord_neigh[2] - coord[2];
                            weight = 1.0 / sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
                            b[0] += weight * delta_x * delta_val;
                            b[1] += weight * delta_y * delta_val;
                            b[2] += weight * delta_z * delta_val;
                        }

                        grad = A_inv * b;
                        for (size_t iDim = 0; iDim < 3; ++iDim)
                        {
                            data_manager->SetPrimitiveGrad(idx_eq, iDim, idx, grad(iDim));
                        }
                    }
                }
            }
        }
    }
    void NSSolverStruct::CalcLimiterVK()
    {
        auto grid = GetGrid();
        auto node = grid->GetNode();
        auto data_manager = GetDataManager();
        auto idx_proxy = GetIdxProxy();
        double eps = 1e-6;
        double vk_coef = 1.0e-5;
        int ni, nj, nk;
        ni = idx_proxy->GetNi();
        nj = idx_proxy->GetNj();
        nk = idx_proxy->GetNk();
        double max_val, min_val;
        int temp_i[6], temp_j[6], temp_k[6];
        for (int idx_eq = 0; idx_eq < 5; idx_eq++)
        {

            for (int k = 1; k < nk - 1; ++k)
            {
                for (int j = 1; j < nj - 1; ++j)
                {
                    for (int i = 1; i < ni - 1; ++i)
                    {
                        int idx = idx_proxy->GetIdx(i, j, k);
                        auto coord = node->GetCoord(i, j, k);
                        temp_i[0] = i - 1, temp_j[0] = j, temp_k[0] = k;
                        temp_i[1] = i + 1, temp_j[1] = j, temp_k[1] = k;
                        temp_i[2] = i, temp_j[2] = j - 1, temp_k[2] = k;
                        temp_i[3] = i, temp_j[3] = j + 1, temp_k[3] = k;
                        temp_i[4] = i, temp_j[4] = j, temp_k[4] = k - 1;
                        temp_i[5] = i, temp_j[5] = j, temp_k[5] = k + 1;
                        max_val = min_val = data_manager->GetPrim(idx_eq, idx);
                        for (int iNeigh = 0; iNeigh < 6; ++iNeigh)
                        {
                            int idx_neigh = idx_proxy->GetIdx(temp_i[iNeigh], temp_j[iNeigh], temp_k[iNeigh]);
                            max_val = Max(max_val, data_manager->GetPrim(idx_eq, idx_neigh));
                            min_val = Min(min_val, data_manager->GetPrim(idx_eq, idx_neigh));
                        }
                        eps = vk_coef * (max_val - min_val);
                        eps = eps * eps + SMALL_NUMBER;
                        double delta_max = max_val - data_manager->GetPrim(idx_eq, idx);
                        double delta_min = min_val - data_manager->GetPrim(idx_eq, idx);
                        double temp_coef = LARGE_NUMBER;
                        for (int iNeigh = 0; iNeigh < 6; ++iNeigh)
                        {
                            int idx_neigh = idx_proxy->GetIdx(temp_i[iNeigh], temp_j[iNeigh], temp_k[iNeigh]);
                            auto coord_neigh = node->GetCoord(temp_i[iNeigh], temp_j[iNeigh], temp_k[iNeigh]);
                            double delta2 = 0.0;
                            for (int iDim = 0; iDim < 3; ++iDim)
                            {
                                delta2 += (coord_neigh[iDim] - coord[iDim]) * data_manager->GetPrimGrad(idx_eq, iDim, idx);
                            }
                            delta2 *= 0.5;
                            if (delta2 > 0.0)
                            {
                                temp_coef = Min(temp_coef, LimiterVK(max_val - data_manager->GetPrim(idx_eq, idx), delta2, eps));
                            }
                            else if (delta2 < 0.0)
                            {
                                temp_coef = Min(temp_coef, LimiterVK(min_val - data_manager->GetPrim(idx_eq, idx), delta2, eps));
                            }
                            else
                            {
                                temp_coef = Min(temp_coef, 1.0);
                            }
                        }
                        data_manager->SetLimiter(idx_eq, idx, temp_coef);
                    }
                }
            }
        }
    }
    void NSSolverStruct::CalcTimeStepLocal()
    {
        auto grid = GetGrid();
        auto data_manager = GetDataManager();
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        auto para = GetPara();
        double cfl = para->GetCflNumber();
        double gamma = GetGas()->GetGamma();
        double min_dt = LARGE_NUMBER;
        for (int k = ks; k <= ke; ++k)
        {
            for (int j = js; j <= je; ++j)
            {
                for (int i = is; i <= ie; ++i)
                {
                    int idx = m_idx_proxy->GetIdx(i, j, k);
                    auto xi = m_node_metrics->GetXi(idx);
                    auto eta = m_node_metrics->GetEta(idx);
                    auto zeta = m_node_metrics->GetZeta(idx);
                    auto jacobi = m_node_metrics->GetJacobian(idx);
                    double c = sqrt(gamma * data_manager->GetPressure(idx) / data_manager->GetDensity(idx));
                    double norm_xi = sqrt(xi[0] * xi[0] + xi[1] * xi[1] + xi[2] * xi[2]);
                    double norm_eta = sqrt(eta[0] * eta[0] + eta[1] * eta[1] + eta[2] * eta[2]);
                    double norm_zeta = sqrt(zeta[0] * zeta[0] + zeta[1] * zeta[1] + zeta[2] * zeta[2]);
                    double u_xi = data_manager->GetVelocity(0, idx) * xi[0] + data_manager->GetVelocity(1, idx) * xi[1] +
                                  data_manager->GetVelocity(2, idx) * xi[2];
                    double u_eta = data_manager->GetVelocity(0, idx) * eta[0] + data_manager->GetVelocity(1, idx) * eta[1] +
                                   data_manager->GetVelocity(2, idx) * eta[2];
                    double u_zeta = data_manager->GetVelocity(0, idx) * zeta[0] +
                                    data_manager->GetVelocity(1, idx) * zeta[1] +
                                    data_manager->GetVelocity(2, idx) * zeta[2];
                    double lamda = abs(u_xi) + abs(u_eta) + c * (norm_xi + norm_eta);
                    if (grid->GetDim() == 3)
                    {
                        lamda += abs(u_zeta) + c * norm_zeta;
                    }
                    lamda = lamda * jacobi;
                    data_manager->SetTimeStep(idx, cfl / lamda);
                    if (data_manager->GetTimeStep(idx) < min_dt)
                    {
                        min_dt = data_manager->GetTimeStep(idx);
                    }
                }
            }
        }
        GlobalData::Update("dt", min_dt);
    }
    void NSSolverStruct::ReduceTimeStep(double &dt)
    {
        auto grid = GetGrid();
        auto data_manager = GetDataManager();
        auto ni = m_idx_proxy->GetNi();
        auto nj = m_idx_proxy->GetNj();
        auto nk = m_idx_proxy->GetNk();
        for (int k = 0; k < nk; ++k)
        {
            for (int j = 0; j < nj; ++j)
            {
                for (int i = 0; i < ni; ++i)
                {
                    int idx = m_idx_proxy->GetIdx(i, j, k);
                    data_manager->SetTimeStep(idx, dt);
                }
            }
        }
    }
    void NSSolverStruct::RungeKutta()
    {
        auto grid = GetGrid();
        auto para = GetPara();
        auto data_manager = GetDataManager();
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        const DArray &rk_coef = para->GetRKCoef();
        int rkStage = rk_coef.size();
        double dt, jacobi;
        for (int iStage = 0; iStage < rkStage; ++iStage)
        {
            CalcResidual();
            for (int i = is; i <= ie; ++i)
            {
                for (int j = js; j <= je; ++j)
                {
                    for (int k = ks; k <= ke; ++k)
                    {
                        int idx = m_idx_proxy->GetIdx(i, j, k);
                        dt = data_manager->GetTimeStep(idx);
                        jacobi = m_node_metrics->GetJacobian(idx);
                        for (int idx_eq = 0; idx_eq < 5; ++idx_eq)
                        {
                            data_manager->SetCons(idx_eq, idx,
                                                  data_manager->GetCons(idx_eq, idx) +
                                                      dt * rk_coef[iStage] * data_manager->GetResidual(idx_eq, idx) * jacobi);
                        }
                    }
                }
            }
        }
    }
    void NSSolverStruct::Prim2Cons()
    {
        auto gas = GetGas();
        auto grid = GetGrid();
        auto data_manager = GetDataManager();
        auto ni = m_idx_proxy->GetNi();
        auto nj = m_idx_proxy->GetNj();
        auto nk = m_idx_proxy->GetNk();
        double prim[5], cons[5];
        for (int k = 0; k < nk; ++k)
        {
            for (int j = 0; j < nj; ++j)
            {
                for (int i = 0; i < ni; ++i)
                {
                    int idx = m_idx_proxy->GetIdx(i, j, k);
                    for (int idx_eq = 0; idx_eq < GetPara()->GetEqNum(); ++idx_eq)
                    {
                        prim[idx_eq] = data_manager->GetPrim(idx_eq, idx);
                    }
                    gas->Prim2Cons(prim, cons);
                    for (int iEqu = 0; iEqu < GetPara()->GetEqNum(); ++iEqu)
                    {
                        data_manager->SetCons(iEqu, idx, cons[iEqu]);
                    }
                }
            }
        }
    }
    void NSSolverStruct::Cons2Prim()
    {
        auto gas = GetGas();
        auto grid = GetGrid();
        auto data_manager = GetDataManager();
        auto ni = m_idx_proxy->GetNi();
        auto nj = m_idx_proxy->GetNj();
        auto nk = m_idx_proxy->GetNk();
        double prim[5], cons[5];
        for (int k = 0; k < nk; ++k)
        {
            for (int j = 0; j < nj; ++j)
            {
                for (int i = 0; i < ni; ++i)
                {
                    int idx = m_idx_proxy->GetIdx(i, j, k);
                    for (int idx_eq = 0; idx_eq < GetPara()->GetEqNum(); ++idx_eq)
                    {
                        cons[idx_eq] = data_manager->GetCons(idx_eq, idx);
                    }
                    gas->Cons2Prim(cons, prim);
                    for (int iEqu = 0; iEqu < GetPara()->GetEqNum(); ++iEqu)
                    {
                        data_manager->SetPrim(iEqu, idx, prim[iEqu]);
                    }
                }
            }
        }
    }
    void NSSolverStruct::ZeroResidual()
    {
        auto grid = GetGrid();
        auto data_manager = GetDataManager();
        auto ni = m_idx_proxy->GetNi();
        auto nj = m_idx_proxy->GetNj();
        auto nk = m_idx_proxy->GetNk();
        for (int k = 0; k < nk; ++k)
        {
            for (int j = 0; j < nj; ++j)
            {
                for (int i = 0; i < ni; ++i)
                {
                    int idx = m_idx_proxy->GetIdx(i, j, k);
                    for (int idx_eq = 0; idx_eq < GetPara()->GetEqNum(); ++idx_eq)
                    {
                        data_manager->SetResidual(idx_eq, idx, 0.0);
                    }
                }
            }
        }
    }

    void NSSolverStruct::CalcMidNode1st()
    {
        auto grid = GetGrid();
        auto data_manager = GetDataManager();
        auto ni = m_idx_proxy->GetNi();
        auto nj = m_idx_proxy->GetNj();
        auto nk = m_idx_proxy->GetNk();
        double value[5];
        double value_left[5], value_right[5];
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    int idx = m_idx_proxy->GetIdx(i, j, k);
                    int idx_left = m_idx_proxy->GetIdx(i - 1, j, k);
                    int idx_right = m_idx_proxy->GetIdx(i + 1, j, k);
                    for (int idx_eq = 0; idx_eq < 5; ++idx_eq)
                    {
                        data_manager->SetMidNodePrimLeft(idx_eq, 0, idx, data_manager->GetPrim(idx_eq, idx_left));
                        data_manager->SetMidNodePrimRight(idx_eq, 0, idx, data_manager->GetPrim(idx_eq, idx_right));
                    }
                    idx_left = m_idx_proxy->GetIdx(i, j - 1, k);
                    idx_right = m_idx_proxy->GetIdx(i, j + 1, k);
                    for (int idx_eq = 0; idx_eq < 5; ++idx_eq)
                    {
                        data_manager->SetMidNodePrimLeft(idx_eq, 1, idx, data_manager->GetPrim(idx_eq, idx_left));
                        data_manager->SetMidNodePrimRight(idx_eq, 1, idx, data_manager->GetPrim(idx_eq, idx_right));
                    }
                    if (grid->GetDim() == 3)
                    {
                        idx_left = m_idx_proxy->GetIdx(i, j, k - 1);
                        idx_right = m_idx_proxy->GetIdx(i, j, k + 1);
                        for (int idx_eq = 0; idx_eq < 5; ++idx_eq)
                        {
                            data_manager->SetMidNodePrimLeft(idx_eq, 2, idx, data_manager->GetPrim(idx_eq, idx_left));
                            data_manager->SetMidNodePrimRight(idx_eq, 2, idx, data_manager->GetPrim(idx_eq, idx_right));
                        }
                    }
                }
            }
        }
    }
    void NSSolverStruct::CalcMidNodeGrad()
    {
        auto grid = GetGrid();
        auto node = grid->GetNode();
        auto data_manager = GetDataManager();
        auto ni = m_idx_proxy->GetNi();
        auto nj = m_idx_proxy->GetNj();
        auto nk = m_idx_proxy->GetNk();
        double value[5];
        double coord_vec[3];
        double value_left[5], value_right[5];
        for (int k = 1; k < nk - 1; ++k)
        {
            for (int j = 1; j < nj - 1; ++j)
            {
                for (int i = 1; i < ni - 1; ++i)
                {
                    int idx = m_idx_proxy->GetIdx(i, j, k);
                    int idx_left = m_idx_proxy->GetIdx(i - 1, j, k);
                    int idx_right = m_idx_proxy->GetIdx(i + 1, j, k);
                    coord_vec[0] = node->GetCoord(i + 1, j, k)[0] - node->GetCoord(i, j, k)[0];
                    coord_vec[1] = node->GetCoord(i + 1, j, k)[1] - node->GetCoord(i, j, k)[1];
                    coord_vec[2] = node->GetCoord(i + 1, j, k)[2] - node->GetCoord(i, j, k)[2];
                    MidNodeGrad(idx, idx_right, coord_vec, value_left, value_right);
                    for (int idx_eq = 0; idx_eq < 5; ++idx_eq)
                    {
                        data_manager->SetMidNodePrim(idx_eq, 0, idx, value_left[idx_eq], value_right[idx_eq]);
                    }
                    coord_vec[0] = node->GetCoord(i, j, k)[0] - node->GetCoord(i - 1, j, k)[0];
                    coord_vec[1] = node->GetCoord(i, j, k)[1] - node->GetCoord(i - 1, j, k)[1];
                    coord_vec[2] = node->GetCoord(i, j, k)[2] - node->GetCoord(i - 1, j, k)[2];
                    MidNodeGrad(idx_left, idx, coord_vec, value_left, value_right);
                    for (int idx_eq = 0; idx_eq < 5; ++idx_eq)
                    {
                        data_manager->SetMidNodePrim(idx_eq, 0, idx, value_left[idx_eq], value_right[idx_eq]);
                    }
                    idx_left = m_idx_proxy->GetIdx(i, j - 1, k);
                    idx_right = m_idx_proxy->GetIdx(i, j + 1, k);
                    coord_vec[0] = node->GetCoord(i, j + 1, k)[0] - node->GetCoord(i, j, k)[0];
                    coord_vec[1] = node->GetCoord(i, j + 1, k)[1] - node->GetCoord(i, j, k)[1];
                    coord_vec[2] = node->GetCoord(i, j + 1, k)[2] - node->GetCoord(i, j, k)[2];
                    MidNodeGrad(idx, idx_right, coord_vec, value_left, value_right);
                    for (int idx_eq = 0; idx_eq < 5; ++idx_eq)
                    {
                        data_manager->SetMidNodePrim(idx_eq, 1, idx, value_left[idx_eq], value_right[idx_eq]);
                    }
                    coord_vec[0] = node->GetCoord(i, j, k)[0] - node->GetCoord(i, j - 1, k)[0];
                    coord_vec[1] = node->GetCoord(i, j, k)[1] - node->GetCoord(i, j - 1, k)[1];
                    coord_vec[2] = node->GetCoord(i, j, k)[2] - node->GetCoord(i, j - 1, k)[2];
                    MidNodeGrad(idx_left, idx, coord_vec, value_left, value_right);
                    for (int idx_eq = 0; idx_eq < 5; ++idx_eq)
                    {
                        data_manager->SetMidNodePrim(idx_eq, 1, idx, value_left[idx_eq], value_right[idx_eq]);
                    }
                    if (grid->GetDim() == 3)
                    {
                        idx_left = m_idx_proxy->GetIdx(i, j, k - 1);
                        idx_right = m_idx_proxy->GetIdx(i, j, k + 1);
                        coord_vec[0] = node->GetCoord(i, j, k + 1)[0] - node->GetCoord(i, j, k)[0];
                        coord_vec[1] = node->GetCoord(i, j, k + 1)[1] - node->GetCoord(i, j, k)[1];
                        coord_vec[2] = node->GetCoord(i, j, k + 1)[2] - node->GetCoord(i, j, k)[2];
                        MidNodeGrad(idx, idx_right, coord_vec, value_left, value_right);
                        for (int idx_eq = 0; idx_eq < 5; ++idx_eq)
                        {
                            data_manager->SetMidNodePrim(idx_eq, 2, idx, value_left[idx_eq], value_right[idx_eq]);
                        }
                        coord_vec[0] = node->GetCoord(i, j, k)[0] - node->GetCoord(i, j, k - 1)[0];
                        coord_vec[1] = node->GetCoord(i, j, k)[1] - node->GetCoord(i, j, k - 1)[1];
                        coord_vec[2] = node->GetCoord(i, j, k)[2] - node->GetCoord(i, j, k - 1)[2];
                        MidNodeGrad(idx_left, idx, coord_vec, value_left, value_right);
                        for (int idx_eq = 0; idx_eq < 5; ++idx_eq)
                        {
                            data_manager->SetMidNodePrim(idx_eq, 2, idx, value_left[idx_eq], value_right[idx_eq]);
                        }
                    }
                }
            }
        }
    }
    void NSSolverStruct::CalcMidNodePrimMUSCL()
    {
        auto grid = GetGrid();
        auto data_manager = GetDataManager();
        int equ_num = GetPara()->GetEqNum();
        auto ni = m_idx_proxy->GetNi();
        auto nj = m_idx_proxy->GetNj();
        auto nk = m_idx_proxy->GetNk();
        // MUSCL整点变量(i-2,i-1,i,i+1,i+2)的值
        double value_temp[5];
        int idx_temp[5];
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        double left_value, right_value;
        ///@note
        /// 从ks-1开始，到ke结束，因为对第一个计算点ks进行通量差分时，需要使用ks-1/2处的值
        for (int k = ks; k <= ke; ++k)
        {
            for (int j = js; j <= je; ++j)
            {
                for (int i = is; i <= ie; ++i)
                {
                    idx_temp[0] = m_idx_proxy->GetIdx(i - 2, j, k);
                    idx_temp[1] = m_idx_proxy->GetIdx(i - 1, j, k);
                    idx_temp[2] = m_idx_proxy->GetIdx(i, j, k);
                    idx_temp[3] = m_idx_proxy->GetIdx(i + 1, j, k);
                    idx_temp[4] = m_idx_proxy->GetIdx(i + 2, j, k);
                    for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
                    {
                        value_temp[0] = data_manager->GetPrim(idx_eq, idx_temp[0]);
                        value_temp[1] = data_manager->GetPrim(idx_eq, idx_temp[1]);
                        value_temp[2] = data_manager->GetPrim(idx_eq, idx_temp[2]);
                        value_temp[3] = data_manager->GetPrim(idx_eq, idx_temp[3]);
                        value_temp[4] = data_manager->GetPrim(idx_eq, idx_temp[4]);
                        MidNodeMUSCL(value_temp, left_value, right_value);
                        data_manager->SetMidNodePrim(idx_eq, 0, idx_temp[2], left_value, right_value);
                    }
                    idx_temp[0] = m_idx_proxy->GetIdx(i, j - 2, k);
                    idx_temp[1] = m_idx_proxy->GetIdx(i, j - 1, k);
                    idx_temp[3] = m_idx_proxy->GetIdx(i, j + 1, k);
                    idx_temp[4] = m_idx_proxy->GetIdx(i, j + 2, k);
                    for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
                    {
                        value_temp[0] = data_manager->GetPrim(idx_eq, idx_temp[0]);
                        value_temp[1] = data_manager->GetPrim(idx_eq, idx_temp[1]);
                        value_temp[2] = data_manager->GetPrim(idx_eq, idx_temp[2]);
                        value_temp[3] = data_manager->GetPrim(idx_eq, idx_temp[3]);
                        value_temp[4] = data_manager->GetPrim(idx_eq, idx_temp[4]);
                        MidNodeMUSCL(value_temp, left_value, right_value);
                        data_manager->SetMidNodePrim(idx_eq, 1, idx_temp[2], left_value, right_value);
                    }
                    if (grid->GetDim() == 3)
                    {
                        idx_temp[0] = m_idx_proxy->GetIdx(i, j, k - 2);
                        idx_temp[1] = m_idx_proxy->GetIdx(i, j, k - 1);
                        idx_temp[3] = m_idx_proxy->GetIdx(i, j, k + 1);
                        idx_temp[4] = m_idx_proxy->GetIdx(i, j, k + 2);
                        for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
                        {
                            value_temp[0] = data_manager->GetPrim(idx_eq, idx_temp[0]);
                            value_temp[1] = data_manager->GetPrim(idx_eq, idx_temp[1]);
                            value_temp[2] = data_manager->GetPrim(idx_eq, idx_temp[2]);
                            value_temp[3] = data_manager->GetPrim(idx_eq, idx_temp[3]);
                            value_temp[4] = data_manager->GetPrim(idx_eq, idx_temp[4]);
                            MidNodeMUSCL(value_temp, left_value, right_value);
                            data_manager->SetMidNodePrim(idx_eq, 2, idx_temp[2], left_value, right_value);
                        }
                    }
                }
            }
        }
    }

    void NSSolverStruct::CalcMidGhostNodePrimMUSCL()
    {
        auto grid = GetGrid();
        auto data_manager = GetDataManager();
        int equ_num = GetPara()->GetEqNum();
        auto ni = m_idx_proxy->GetNi();
        auto nj = m_idx_proxy->GetNj();
        auto nk = m_idx_proxy->GetNk();
        double value[5];
        double value_left[5], value_right[5];
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        int idx_temp[5];
        double left_value, right_value;
        int i, j, k;
        for (k = ks; k <= ke; ++k)
        {
            for (j = js; j <= je; ++j)
            {
                // i=is-1处的值
                // 不存在i-2处的值且不参与计算,用i-1处的值代替，不影响计算结果
                idx_temp[0] = m_idx_proxy->GetIdx(is - 2, j, k);
                idx_temp[1] = m_idx_proxy->GetIdx(is - 2, j, k);
                idx_temp[2] = m_idx_proxy->GetIdx(is - 1, j, k);
                idx_temp[3] = m_idx_proxy->GetIdx(is, j, k);
                idx_temp[4] = m_idx_proxy->GetIdx(is + 1, j, k);
                for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
                {
                    value[0] = data_manager->GetPrim(idx_eq, idx_temp[0]);
                    value[1] = data_manager->GetPrim(idx_eq, idx_temp[1]);
                    value[2] = data_manager->GetPrim(idx_eq, idx_temp[2]);
                    value[3] = data_manager->GetPrim(idx_eq, idx_temp[3]);
                    value[4] = data_manager->GetPrim(idx_eq, idx_temp[4]);
                    MidNodeMUSCL(value, left_value, right_value);
                    data_manager->SetMidNodePrim(idx_eq, 0, idx_temp[2], left_value, right_value);
                }
            }
        }

        for (k = ks; k <= ke; ++k)
        {
            for (i = is; i <= ie; ++i)
            {
                // j=js-1处的值，(js-1/2)右值会用到，但(js-1/2)左值在计算中不会使用
                // 不存在j-2处的值,用j-1处的值代替，不影响计算结果
                idx_temp[0] = m_idx_proxy->GetIdx(i, js - 2, k);
                idx_temp[1] = m_idx_proxy->GetIdx(i, js - 2, k);
                idx_temp[2] = m_idx_proxy->GetIdx(i, js - 1, k);
                idx_temp[3] = m_idx_proxy->GetIdx(i, js, k);
                idx_temp[4] = m_idx_proxy->GetIdx(i, js + 1, k);
                for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
                {
                    value[0] = data_manager->GetPrim(idx_eq, idx_temp[0]);
                    value[1] = data_manager->GetPrim(idx_eq, idx_temp[1]);
                    value[2] = data_manager->GetPrim(idx_eq, idx_temp[2]);
                    value[3] = data_manager->GetPrim(idx_eq, idx_temp[3]);
                    value[4] = data_manager->GetPrim(idx_eq, idx_temp[4]);
                    MidNodeMUSCL(value, left_value, right_value);
                    data_manager->SetMidNodePrim(idx_eq, 1, idx_temp[2], left_value, right_value);
                }
            }
        }

        if (grid->GetDim() == 3)
        {
            for (j = js; j <= je; ++j)
            {
                for (i = is; i <= ie; ++i)
                {
                    // k=ks-1处的值，(ks-1/2)右值会用到，但(ks-1/2)左值在计算中不会使用
                    // 不存在k-2处的值,用k-1处的值代替，不影响计算结果
                    idx_temp[0] = m_idx_proxy->GetIdx(i, j, ks - 2);
                    idx_temp[1] = m_idx_proxy->GetIdx(i, j, ks - 2);
                    idx_temp[2] = m_idx_proxy->GetIdx(i, j, ks - 1);
                    idx_temp[3] = m_idx_proxy->GetIdx(i, j, ks);
                    idx_temp[4] = m_idx_proxy->GetIdx(i, j, ks + 1);
                    for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
                    {
                        value[0] = data_manager->GetPrim(idx_eq, idx_temp[0]);
                        value[1] = data_manager->GetPrim(idx_eq, idx_temp[1]);
                        value[2] = data_manager->GetPrim(idx_eq, idx_temp[2]);
                        value[3] = data_manager->GetPrim(idx_eq, idx_temp[3]);
                        value[4] = data_manager->GetPrim(idx_eq, idx_temp[4]);
                        MidNodeMUSCL(value, left_value, right_value);
                        data_manager->SetMidNodePrim(idx_eq, 2, idx_temp[2], left_value, right_value);
                    }
                }
            }
        }
    }

    void NSSolverStruct::CalcMidNodePrimWCNS5()
    {
        auto grid = GetGrid();
        auto data_manager = GetDataManager();
        int equ_num = GetPara()->GetEqNum();
        auto ni = m_idx_proxy->GetNi();
        auto nj = m_idx_proxy->GetNj();
        auto nk = m_idx_proxy->GetNk();
        // WCNS整点变量(i-2,i-1,i,i+1,i+2)的值
        double value_temp[5];
        int idx_temp[5];
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        double left_value, right_value;
        ///@note
        /// 从ks-1开始，到ke结束，因为对第一个计算点ks进行通量差分时，需要使用ks-1/2处的值
        for (int k = ks - 1; k <= ke + 1; ++k)
        {
            for (int j = js - 1; j <= je + 1; ++j)
            {
                for (int i = is - 1; i <= ie + 1; ++i)
                {
                    idx_temp[0] = m_idx_proxy->GetIdx(i - 2, j, k);
                    idx_temp[1] = m_idx_proxy->GetIdx(i - 1, j, k);
                    idx_temp[2] = m_idx_proxy->GetIdx(i, j, k);
                    idx_temp[3] = m_idx_proxy->GetIdx(i + 1, j, k);
                    idx_temp[4] = m_idx_proxy->GetIdx(i + 2, j, k);
                    for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
                    {
                        value_temp[0] = data_manager->GetPrim(idx_eq, idx_temp[0]);
                        value_temp[1] = data_manager->GetPrim(idx_eq, idx_temp[1]);
                        value_temp[2] = data_manager->GetPrim(idx_eq, idx_temp[2]);
                        value_temp[3] = data_manager->GetPrim(idx_eq, idx_temp[3]);
                        value_temp[4] = data_manager->GetPrim(idx_eq, idx_temp[4]);
                        MidNodeWCNS5(value_temp, left_value, right_value);
                        data_manager->SetMidNodePrimLeft(idx_eq, 0, idx_temp[2], left_value);
                        data_manager->SetMidNodePrimRight(idx_eq, 0, idx_temp[1], right_value);
                    }
                    idx_temp[0] = m_idx_proxy->GetIdx(i, j - 2, k);
                    idx_temp[1] = m_idx_proxy->GetIdx(i, j - 1, k);
                    idx_temp[3] = m_idx_proxy->GetIdx(i, j + 1, k);
                    idx_temp[4] = m_idx_proxy->GetIdx(i, j + 2, k);
                    for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
                    {
                        value_temp[0] = data_manager->GetPrim(idx_eq, idx_temp[0]);
                        value_temp[1] = data_manager->GetPrim(idx_eq, idx_temp[1]);
                        value_temp[2] = data_manager->GetPrim(idx_eq, idx_temp[2]);
                        value_temp[3] = data_manager->GetPrim(idx_eq, idx_temp[3]);
                        value_temp[4] = data_manager->GetPrim(idx_eq, idx_temp[4]);
                        MidNodeWCNS5(value_temp, left_value, right_value);
                        data_manager->SetMidNodePrimLeft(idx_eq, 1, idx_temp[2], left_value);
                        data_manager->SetMidNodePrimRight(idx_eq, 1, idx_temp[1], right_value);
                    }
                    if (grid->GetDim() == 3)
                    {
                        idx_temp[0] = m_idx_proxy->GetIdx(i, j, k - 2);
                        idx_temp[1] = m_idx_proxy->GetIdx(i, j, k - 1);
                        idx_temp[3] = m_idx_proxy->GetIdx(i, j, k + 1);
                        idx_temp[4] = m_idx_proxy->GetIdx(i, j, k + 2);
                        for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
                        {
                            value_temp[0] = data_manager->GetPrim(idx_eq, idx_temp[0]);
                            value_temp[1] = data_manager->GetPrim(idx_eq, idx_temp[1]);
                            value_temp[2] = data_manager->GetPrim(idx_eq, idx_temp[2]);
                            value_temp[3] = data_manager->GetPrim(idx_eq, idx_temp[3]);
                            value_temp[4] = data_manager->GetPrim(idx_eq, idx_temp[4]);
                            MidNodeWCNS5(value_temp, left_value, right_value);
                            data_manager->SetMidNodePrimLeft(idx_eq, 2, idx_temp[2], left_value);
                            data_manager->SetMidNodePrimRight(idx_eq, 2, idx_temp[1], right_value);
                        }
                    }
                }
            }
        }
    }

    void NSSolverStruct::CalcMidGhostNodePrimWCNS5()
    {
        auto grid = GetGrid();
        auto data_manager = GetDataManager();
        int equ_num = GetPara()->GetEqNum();
        auto ni = m_idx_proxy->GetNi();
        auto nj = m_idx_proxy->GetNj();
        auto nk = m_idx_proxy->GetNk();
        double value_left[4], value_right[4];
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        // 左侧的模板点(0,1,2,3),右侧的模板点(N-4,N-3,N-2,N-1)
        int idx_temp_left[4], idx_temp_right[4];
        double left_value, right_value;
        int i, j, k;
        for (k = ks; k <= ke; ++k)
        {
            for (j = js; j <= je; ++j)
            {
                idx_temp_left[0] = m_idx_proxy->GetIdx(0, j, k);
                idx_temp_left[1] = m_idx_proxy->GetIdx(1, j, k);
                idx_temp_left[2] = m_idx_proxy->GetIdx(2, j, k);
                idx_temp_left[3] = m_idx_proxy->GetIdx(3, j, k);
                idx_temp_right[0] = m_idx_proxy->GetIdx(ni - 4, j, k);
                idx_temp_right[1] = m_idx_proxy->GetIdx(ni - 3, j, k);
                idx_temp_right[2] = m_idx_proxy->GetIdx(ni - 2, j, k);
                idx_temp_right[3] = m_idx_proxy->GetIdx(ni - 1, j, k);
                for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
                {
                    value_left[0] = data_manager->GetPrim(idx_eq, idx_temp_left[0]);
                    value_left[1] = data_manager->GetPrim(idx_eq, idx_temp_left[1]);
                    value_left[2] = data_manager->GetPrim(idx_eq, idx_temp_left[2]);
                    value_left[3] = data_manager->GetPrim(idx_eq, idx_temp_left[3]);
                    value_right[0] = data_manager->GetPrim(idx_eq, idx_temp_right[0]);
                    value_right[1] = data_manager->GetPrim(idx_eq, idx_temp_right[1]);
                    value_right[2] = data_manager->GetPrim(idx_eq, idx_temp_right[2]);
                    value_right[3] = data_manager->GetPrim(idx_eq, idx_temp_right[3]);
                    left_value = right_value =
                        (5 * value_left[0] + 15 * value_left[1] - 5 * value_left[2] + value_left[3]) / 16;
                    data_manager->SetMidNodePrim(idx_eq, 0, idx_temp_left[0], left_value, right_value);
                    left_value = right_value =
                        (-value_left[0] + 9 * value_left[1] + 9 * value_left[2] - value_left[3]) / 16;
                    data_manager->SetMidNodePrim(idx_eq, 0, idx_temp_left[1], left_value, right_value);
                    left_value = right_value =
                        (5 * value_right[0] + 15 * value_right[1] - 5 * value_right[2] + value_right[3]) / 16;
                    data_manager->SetMidNodePrim(idx_eq, 0, idx_temp_right[1], left_value, right_value);
                    left_value = right_value =
                        (-value_right[0] + 9 * value_right[1] + 9 * value_right[2] - value_right[3]) / 16;
                    data_manager->SetMidNodePrim(idx_eq, 0, idx_temp_right[2], left_value, right_value);
                }
            }
        }

        for (k = ks; k <= ke; ++k)
        {
            for (i = is; i <= ie; ++i)
            {
                // j=js-1处的值，(js-1/2)右值会用到，但(js-1/2)左值在计算中不会使用
                // 不存在j-2处的值,用j-1处的值代替，不影响计算结果
                idx_temp_left[0] = m_idx_proxy->GetIdx(i, 0, k);
                idx_temp_left[1] = m_idx_proxy->GetIdx(i, 1, k);
                idx_temp_left[2] = m_idx_proxy->GetIdx(i, 2, k);
                idx_temp_left[3] = m_idx_proxy->GetIdx(i, 3, k);
                idx_temp_right[0] = m_idx_proxy->GetIdx(i, nj - 4, k);
                idx_temp_right[1] = m_idx_proxy->GetIdx(i, nj - 3, k);
                idx_temp_right[2] = m_idx_proxy->GetIdx(i, nj - 2, k);
                idx_temp_right[3] = m_idx_proxy->GetIdx(i, nj - 1, k);
                for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
                {
                    value_left[0] = data_manager->GetPrim(idx_eq, idx_temp_left[0]);
                    value_left[1] = data_manager->GetPrim(idx_eq, idx_temp_left[1]);
                    value_left[2] = data_manager->GetPrim(idx_eq, idx_temp_left[2]);
                    value_left[3] = data_manager->GetPrim(idx_eq, idx_temp_left[3]);
                    value_right[0] = data_manager->GetPrim(idx_eq, idx_temp_right[0]);
                    value_right[1] = data_manager->GetPrim(idx_eq, idx_temp_right[1]);
                    value_right[2] = data_manager->GetPrim(idx_eq, idx_temp_right[2]);
                    value_right[3] = data_manager->GetPrim(idx_eq, idx_temp_right[3]);
                    left_value = right_value =
                        (5 * value_left[0] + 15 * value_left[1] - 5 * value_left[2] + value_left[3]) / 16;
                    data_manager->SetMidNodePrim(idx_eq, 1, idx_temp_left[0], left_value, right_value);
                    left_value = right_value =
                        (-value_left[0] + 9 * value_left[1] + 9 * value_left[2] - value_left[3]) / 16;
                    data_manager->SetMidNodePrim(idx_eq, 1, idx_temp_left[1], left_value, right_value);
                    left_value = right_value =
                        (5 * value_right[0] + 15 * value_right[1] - 5 * value_right[2] + value_right[3]) / 16;
                    data_manager->SetMidNodePrim(idx_eq, 1, idx_temp_right[1], left_value, right_value);
                    left_value = right_value =
                        (-value_right[0] + 9 * value_right[1] + 9 * value_right[2] - value_right[3]) / 16;
                    data_manager->SetMidNodePrim(idx_eq, 1, idx_temp_right[2], left_value, right_value);
                }
            }
        }

        if (grid->GetDim() == 3)
        {
            for (j = js; j <= je; ++j)
            {
                for (i = is; i <= ie; ++i)
                {
                    // k=ks-1处的值，(ks-1/2)右值会用到，但(ks-1/2)左值在计算中不会使用
                    // 不存在k-2处的值,用k-1处的值代替，不影响计算结果
                    idx_temp_left[0] = m_idx_proxy->GetIdx(i, j, 0);
                    idx_temp_left[1] = m_idx_proxy->GetIdx(i, j, 1);
                    idx_temp_left[2] = m_idx_proxy->GetIdx(i, j, 2);
                    idx_temp_left[3] = m_idx_proxy->GetIdx(i, j, 3);
                    idx_temp_right[0] = m_idx_proxy->GetIdx(i, j, nk - 4);
                    idx_temp_right[1] = m_idx_proxy->GetIdx(i, j, nk - 3);
                    idx_temp_right[2] = m_idx_proxy->GetIdx(i, j, nk - 2);
                    idx_temp_right[3] = m_idx_proxy->GetIdx(i, j, nk - 1);
                    for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
                    {
                        value_left[0] = data_manager->GetPrim(idx_eq, idx_temp_left[0]);
                        value_left[1] = data_manager->GetPrim(idx_eq, idx_temp_left[1]);
                        value_left[2] = data_manager->GetPrim(idx_eq, idx_temp_left[2]);
                        value_left[3] = data_manager->GetPrim(idx_eq, idx_temp_left[3]);
                        value_right[0] = data_manager->GetPrim(idx_eq, idx_temp_right[0]);
                        value_right[1] = data_manager->GetPrim(idx_eq, idx_temp_right[1]);
                        value_right[2] = data_manager->GetPrim(idx_eq, idx_temp_right[2]);
                        value_right[3] = data_manager->GetPrim(idx_eq, idx_temp_right[3]);
                        left_value = right_value =
                            (5 * value_left[0] + 15 * value_left[1] - 5 * value_left[2] + value_left[3]) / 16;
                        data_manager->SetMidNodePrim(idx_eq, 2, idx_temp_left[0], left_value, right_value);
                        left_value = right_value =
                            (-value_left[0] + 9 * value_left[1] + 9 * value_left[2] - value_left[3]) / 16;
                        data_manager->SetMidNodePrim(idx_eq, 2, idx_temp_left[1], left_value, right_value);
                        left_value = right_value =
                            (5 * value_right[0] + 15 * value_right[1] - 5 * value_right[2] + value_right[3]) / 16;
                        data_manager->SetMidNodePrim(idx_eq, 2, idx_temp_right[1], left_value, right_value);
                        left_value = right_value =
                            (-value_right[0] + 9 * value_right[1] + 9 * value_right[2] - value_right[3]) / 16;
                        data_manager->SetMidNodePrim(idx_eq, 2, idx_temp_right[2], left_value, right_value);
                    }
                }
            }
        }
    }

    void NSSolverStruct::MidNodeGrad(int idx_left, int idx_right, double *coord_vec, double *value_left,
                                     double *value_right)
    {
        int equ_num = GetPara()->GetEqNum();
        auto data_manager = GetDataManager();
        for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
        {
            value_left[idx_eq] =
                data_manager->GetPrim(idx_eq, idx_left) + 0.5 * data_manager->GetLimiter(idx_eq, idx_left) *
                                                              (coord_vec[0] * data_manager->GetPrimGrad(idx_eq, 0, idx_left) +
                                                               coord_vec[1] * data_manager->GetPrimGrad(idx_eq, 1, idx_left) +
                                                               coord_vec[2] * data_manager->GetPrimGrad(idx_eq, 2, idx_left));
            value_right[idx_eq] =
                data_manager->GetPrim(idx_eq, idx_right) - 0.5 * data_manager->GetLimiter(idx_eq, idx_right) *
                                                               (coord_vec[0] * data_manager->GetPrimGrad(idx_eq, 0, idx_right) +
                                                                coord_vec[1] * data_manager->GetPrimGrad(idx_eq, 1, idx_right) +
                                                                coord_vec[2] * data_manager->GetPrimGrad(idx_eq, 2, idx_right));
        }
    }
    void NSSolverStruct::MidNode1st(int index_left, int index_right, double *value_rec_left, double *value_rec_right)
    {
        int equ_num = GetPara()->GetEqNum();
        auto grid = GetGrid();
        auto data_manager = GetDataManager();
        for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
        {
            value_rec_left[idx_eq] = data_manager->GetPrim(idx_eq, index_left);
            value_rec_right[idx_eq] = data_manager->GetPrim(idx_eq, index_right);
        }
    }

    void NSSolverStruct::MidNodeMUSCL(const double *value, double &value_left, double &value_right)
    {
        auto para = GetPara();
        auto data_manager = GetDataManager();
        double k = -1.0;
        double delta_plus = value[3] - value[2];
        double delta_minus = value[2] - value[1];
        double delta = LimiterVanLeer(delta_plus, delta_minus);
        value_left = value[2] + 0.125 * ((1.0 - k) * delta + (1.0 + k) * delta);
        delta_plus = value[4] - value[3];
        delta_minus = value[3] - value[2];
        delta = LimiterVanLeer(delta_plus, delta_minus);
        value_right = value[3] - 0.125 * ((1.0 - k) * delta + (1.0 + k) * delta);
        if (isnan(value_left) || isnan(value_right))
        {
            Log::error("MUSCL interpolation failed!");
            Log::error("value[0] = {}, value[1] = {}, value[2] = {}, value[3] = {}, "
                       "value[4] = {}",
                       value[0], value[1], value[2], value[3], value[4]);
            exit(0);
        }
    }

    // 参考《计算空气动力学》pp.84
    void NSSolverStruct::MidNodeWCNS5(const double *value, double &value_left, double &value_right)
    {
        double h = 0.025;
        double f_left, f_right;
        double f[3];
        double omega_left[3], omega_right[5];
        double s_left, s_right;
        double s[3];
        double beta_left[3], beta_right[3];
        double c_left[3], c_right[3];
        c_left[0] = 1.0 / 16.0;
        c_left[1] = 10.0 / 16.0;
        c_left[2] = 5.0 / 16.0;
        c_right[0] = 5.0 / 16.0;
        c_right[1] = 10.0 / 16.0;
        c_right[2] = 1.0 / 16.0;
        f[0] = 0.5 * (value[0] - 4.0 * value[1] + 3.0 * value[2]) / h;
        f[1] = 0.5 * (value[3] - value[1]) / h;
        f[2] = 0.5 * (-3.0 * value[2] + 4.0 * value[3] - value[4]) / h;
        s[0] = (value[0] - 2.0 * value[1] + value[2]) / h / h;
        s[1] = (value[1] - 2.0 * value[2] + value[3]) / h / h;
        s[2] = (value[2] - 2.0 * value[3] + value[4]) / h / h;
        double IS[3];
        for (int i = 0; i < 3; ++i)
        {
            IS[i] = pow(f[i] * h, 2) + pow(s[i] * h * h, 2);
        }
        double eps = 1.0e-6;
        for (int i = 0; i < 3; ++i)
        {
            beta_left[i] = c_left[i] / pow(IS[i] + eps, 2);
            beta_right[i] = c_right[i] / pow(IS[i] + eps, 2);
        }
        for (int i = 0; i < 3; ++i)
        {
            omega_left[i] = beta_left[i] / (beta_left[0] + beta_left[1] + beta_left[2]);
            omega_right[i] = beta_right[i] / (beta_right[0] + beta_right[1] + beta_right[2]);
        }
        f_left = f_right = s_left = s_right = 0.0;
        for (int i = 0; i < 3; ++i)
        {
            f_left += omega_left[i] * f[i];
            s_left += omega_left[i] * s[i];
            f_right += omega_right[i] * f[i];
            s_right += omega_right[i] * s[i];
        }
        value_left = value[2] + 0.5 * h * f_left + 0.125 * h * h * s_left;
        value_right = value[2] - 0.5 * h * f_right + 0.125 * h * h * s_right;
    }

    void NSSolverStruct::CalcInviscidResidual()
    {
        auto para = GetPara();
        auto inter_scheme = para->GetInterSchme();
        if (inter_scheme == InterpolationScheme::FirstOrder)
        {
            CalcInviscidResidual1st();
        }
        else if (inter_scheme == InterpolationScheme::MUSCL)
        {
            CalcInviscidResidualMUSCL();
        }
        else if (inter_scheme == InterpolationScheme::WCNS5)
        {
            CalcInviscidResidualWCNS5();
        }
        else if (inter_scheme == InterpolationScheme::Grad)
        {
            CalcInviscidResidualGrad();
        }
        else
        {
            Log::error("Interpolation scheme not found!");
            CalcInviscidResidual1st();
        }
    }
    void NSSolverStruct::CalcInviscidResidual1st()
    {
        CalcMidNode1st();
        CalcMidGhostNodePrimMUSCL();
        FluxDifference2nd();
    }
    void NSSolverStruct::CalcInviscidResidualGrad()
    {
        CalcMidNodeGrad();
        CalcMidGhostNodePrimMUSCL();
        FluxDifference2nd();
    }
    void NSSolverStruct::CalcInviscidResidualMUSCL()
    {
        CalcMidNodePrimMUSCL();
        CalcMidGhostNodePrimMUSCL();
        FluxDifference2nd();
    }
    void NSSolverStruct::CalcInviscidResidualWCNS5()
    {
        CalcMidNodePrimWCNS5();
        CalcMidGhostNodePrimWCNS5();
        FluxDifference6th();
    }
    void NSSolverStruct::CalcViscousResidual()
    {
        // todo
    }
    void NSSolverStruct::CalcViscousFlux()
    {
        // todo
    }
    void NSSolverStruct::CalcViscousFluxGrad()
    {
        // todo
    }
    void NSSolverStruct::SourceResidual()
    {
        // todo
    }
    void NSSolverStruct::BoundaryCondition()
    {
        auto grid = GetGrid();
        auto data_manager = GetDataManager();
        auto bound_map = grid->GetBoundMap();
        for (auto &boundary : bound_map->GetBoundaryMap())
        {
            auto &bound_name = boundary.first;
            auto &bound = boundary.second;
            if (bound_name == "hole")
                continue;
            if (bound_name == "riemann")
            {
#pragma omp parallel for
                for (int iBound = 0; iBound < bound.size(); ++iBound)
                {
                    RiemannBC(bound[iBound]);
                }
            }
            else if (bound_name == "inlet")
            {
#pragma omp parallel for
                for (int iBound = 0; iBound < bound.size(); ++iBound)
                {
                    InletBC(bound[iBound]);
                }
            }
            else if (bound_name == "outlet")
            {
#pragma omp parallel for
                for (int iBound = 0; iBound < bound.size(); ++iBound)
                {
                    OutletBC(bound[iBound]);
                }
            }
            else if (bound_name == "wall")
            {
#pragma omp parallel for
                for (int iBound = 0; iBound < bound.size(); ++iBound)
                {
                    WallBC(bound[iBound]);
                }
            }
            else
            {
                Log::error("Boundary condition not found");
                exit(0);
            }
        }
    }
    void NSSolverStruct::InletBC(BoundStruct &bound)
    {
        auto grid = GetGrid();
        auto data_manager = GetDataManager();
        int ghost_size = grid->GetGhostLevel();
        int i_bound, j_bound, k_bound;
        bound.GetIdxBound(i_bound, j_bound, k_bound);
        auto bound_direction = bound.GetDirection();
        double prim_far[5];
        prim_far[0] = GetPara()->GetInflowDensity();
        prim_far[1] = GetPara()->GetInflowVelocityX();
        prim_far[2] = GetPara()->GetInflowVelocityY();
        prim_far[3] = GetPara()->GetInflowVelocityZ();
        prim_far[4] = GetPara()->GetInflowPressure();
        prim_far[0] = 1.0;
        prim_far[1] = 0.0;
        prim_far[2] = 0.0;
        prim_far[3] = 0.0;
        prim_far[4] = 1.0;
        double cons_far[5];
        GetGas()->Prim2Cons(prim_far, cons_far);
        int i_ghost, j_ghost, k_ghost;
        int idx_ghost;
        for (int iGhost = 1; iGhost <= ghost_size; ++iGhost)
        {
            i_ghost = i_bound + iGhost * bound_direction[0];
            j_ghost = j_bound + iGhost * bound_direction[1];
            k_ghost = k_bound + iGhost * bound_direction[2];
            idx_ghost = m_idx_proxy->GetIdx(i_ghost, j_ghost, k_ghost);
            for (int idx_eq = 0; idx_eq < 5; ++idx_eq)
            {
                data_manager->SetPrim(idx_eq, idx_ghost, prim_far[idx_eq]);
                data_manager->SetCons(idx_eq, idx_ghost, cons_far[idx_eq]);
            }
        }
    }
    void NSSolverStruct::OutletBC(BoundStruct &bound)
    {
        auto grid = GetGrid();
        auto data_manager = GetDataManager();
        int ghost_size = grid->GetGhostLevel();
        int i_bound, j_bound, k_bound;
        bound.GetIdxBound(i_bound, j_bound, k_bound);
        int idx_bound = m_idx_proxy->GetIdx(i_bound, j_bound, k_bound);
        int i_ghost, j_ghost, k_ghost;
        for (int iGhost = 1; iGhost <= ghost_size; ++iGhost)
        {
            i_ghost = i_bound + iGhost * bound.GetDirection()[0];
            j_ghost = j_bound + iGhost * bound.GetDirection()[1];
            k_ghost = k_bound + iGhost * bound.GetDirection()[2];
            int idx_ghost = m_idx_proxy->GetIdx(i_ghost, j_ghost, k_ghost);
            for (int idx_eq = 0; idx_eq < 5; ++idx_eq)
            {
                data_manager->SetPrim(idx_eq, idx_ghost, data_manager->GetPrim(idx_eq, idx_bound));
                data_manager->SetCons(idx_eq, idx_ghost, data_manager->GetCons(idx_eq, idx_bound));
            }
        }
    }
    void NSSolverStruct::WallBC(BoundStruct &bound)
    {
        auto grid = GetGrid();
        auto data_manager = GetDataManager();
        int ghost_size = grid->GetGhostLevel();
        auto norm_bnd = bound.GetNormBound();
        int i_bound, j_bound, k_bound;
        bound.GetIdxBound(i_bound, j_bound, k_bound);
        int idx_bound = m_idx_proxy->GetIdx(i_bound, j_bound, k_bound);
        int i_ref, j_ref, k_ref;
        int i_ghost, j_ghost, k_ghost;
        for (int iGhost = 1; iGhost <= ghost_size; ++iGhost)
        {
            i_ghost = i_bound + iGhost * bound.GetDirection()[0];
            j_ghost = j_bound + iGhost * bound.GetDirection()[1];
            k_ghost = k_bound + iGhost * bound.GetDirection()[2];
            int idx_ghost = m_idx_proxy->GetIdx(i_ghost, j_ghost, k_ghost);
            i_ref = i_ghost + bound.GetDirection()[0];
            j_ref = j_ghost + bound.GetDirection()[1];
            k_ref = k_ghost + bound.GetDirection()[2];
            int idx_ref = m_idx_proxy->GetIdx(i_ref, j_ref, k_ref);
            double prim_ghost[5];
            for (int idx_eq = 0; idx_eq < data_manager->GetEqNum(); ++idx_eq)
            {
                prim_ghost[idx_eq] = data_manager->GetPrim(idx_eq, idx_ref);
            }
            double vel_ref[3] = {prim_ghost[1], prim_ghost[2], prim_ghost[3]};
            double vn_ref = vel_ref[0] * norm_bnd[0] + vel_ref[1] * norm_bnd[1] + vel_ref[2] * norm_bnd[2];
            prim_ghost[1] = prim_ghost[1] - 2.0 * vn_ref * norm_bnd[0];
            prim_ghost[2] = prim_ghost[2] - 2.0 * vn_ref * norm_bnd[1];
            prim_ghost[3] = prim_ghost[3] - 2.0 * vn_ref * norm_bnd[2];
            double cons_ghost[5];
            GetGas()->Prim2Cons(prim_ghost, cons_ghost);
            for (int idx_eq = 0; idx_eq < 5; ++idx_eq)
            {
                data_manager->SetPrim(idx_eq, idx_ghost, prim_ghost[idx_eq]);
                data_manager->SetCons(idx_eq, idx_bound, cons_ghost[idx_eq]);
            }
        }
    }
    void NSSolverStruct::RiemannBC(BoundStruct &bound)
    {
    }
    void NSSolverStruct::SymmetryBC(BoundStruct &bound)
    {
    }
    void NSSolverStruct::CheckPrimtive()
    {
    }
    void NSSolverStruct::CheckResidual()
    {
    }
    void NSSolverStruct::FixPrimtive()
    {
    }
    void NSSolverStruct::BackupField(std::string &back_folder)
    {
    }
    GridStruct *NSSolverStruct::GetGrid()
    {
        return static_cast<GridStruct *>(FlowSolver::GetGrid());
    }

    Metrics *NSSolverStruct::GetMidMetricsI()
    {
        return m_metrics_half_i;
    }

    Metrics *NSSolverStruct::GetMidMetricsJ()
    {
        return m_metrics_half_j;
    }

    Metrics *NSSolverStruct::GetMidMetricsK()
    {
        return m_metrics_half_k;
    }

    double NSSolverStruct::NodeDifferece2nd(double *mid_data)
    {
        return 0.5 * (mid_data[1] - mid_data[0]);
    }

    double NSSolverStruct::NodeDifferece4th(double *mid_data)
    {
        return 27.0 / 24.0 * (mid_data[2] - mid_data[1]) - 1.0 / 24.0 * (mid_data[3] - mid_data[0]);
    }

    double NSSolverStruct::NodeDifferece4thLeft(double *mid_data)
    {
        return 1.0 / 24.0 *
               (22.0 * mid_data[0] - 17.0 * mid_data[1] - 9.0 * mid_data[2] + 5.0 * mid_data[3] - 1.0 * mid_data[4]);
    }

    double NSSolverStruct::NodeDifferece4thRight(double *mid_data)
    {
        return -1.0 / 24.0 *
               (22.0 * mid_data[0] - 17.0 * mid_data[1] - 9.0 * mid_data[2] + 5.0 * mid_data[3] - 1.0 * mid_data[4]);
    }

    double NSSolverStruct::NodeDifferece6th(double *mid_data)
    {
        return 75.0 / 64.0 * (mid_data[3] - mid_data[2]) - 25.0 / 384.0 * (mid_data[4] - mid_data[1]) +
               3.0 / 640.0 * (mid_data[5] - mid_data[0]);
    }

    double NSSolverStruct::MidNodeInter2nd(double *node_data)
    {
        return 0.5 * (node_data[0] + node_data[1]);
    }

    double NSSolverStruct::MidNodeInter2ndLeft(double *node_data)
    {
        return 1.5 * node_data[0] - 0.5 * node_data[1];
    }

    double NSSolverStruct::MidNodeInter2ndRight(double *node_data)
    {
        return -0.5 * node_data[0] + 1.5 * node_data[1];
    }

    double NSSolverStruct::MidNodeInter4th(double *node_data)
    {
        return 1.0 / 16.0 * (-1.0 * node_data[0] + 9.0 * node_data[1] + 9.0 * node_data[2] - 1.0 * node_data[3]);
    }

    double NSSolverStruct::MidNodeInter4thLeft1(double *node_data)
    {
        return 1.0 / 16.0 * (-5.0 * node_data[0] + 21.0 * node_data[1] - 35.0 * node_data[2] + 35.0 * node_data[3]);
    }

    double NSSolverStruct::MidNodeInter4thLeft2(double *node_data)
    {
        return 1.0 / 16.0 * (1.0 * node_data[0] - 5.0 * node_data[1] + 15.0 * node_data[2] + 5.0 * node_data[3]);
    }

    double NSSolverStruct::MidNodeInter4thRight1(double *node_data)
    {
        return 1.0 / 16.0 * (35.0 * node_data[0] - 35.0 * node_data[1] + 21.0 * node_data[2] - 5.0 * node_data[3]);
    }

    double NSSolverStruct::MidNodeInter4thRight2(double *node_data)
    {
        return 1.0 / 16.0 * (5.0 * node_data[0] + 15.0 * node_data[1] - 5.0 * node_data[2] + 1.0 * node_data[3]);
    }

    double NSSolverStruct::MidNodeInter6th(double *node_data)
    {
        return 75.0 / 128.0 * (node_data[2] + node_data[3]) - 25.0 / 256.0 * (node_data[1] + node_data[4]) +
               3.0 / 256.0 * (node_data[0] + node_data[5]);
    }

} // namespace  zaran