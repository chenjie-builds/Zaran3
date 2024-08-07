#include "NSSolverStruct.h"
namespace  zaran
{
    NSSolverStruct::NSSolverStruct(int index, string name, FlowSolverPara* para, GridStruct* grid, DataManagerNSStruct* data_manager)
        :NSSolver(index, name, para, grid, data_manager)
    {
        int ni, nj, nk;
        m_idx_proxy = new StructIdxProxy(grid);
        ni = m_idx_proxy->GetNi();
        nj = m_idx_proxy->GetNj();
        nk = m_idx_proxy->GetNk();
        int node_num = ni * nj * nk;
        m_node_metrics = new Metrics(node_num);
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
    }
    DataManagerNSStruct* NSSolverStruct::GetDataManager()
    {
        return static_cast<DataManagerNSStruct*>(NSSolver::GetDataManager());
    }
    void NSSolverStruct::InitFieldFarFlow()
    {
        auto grid = GetGrid();
        auto data_manager = GetDataManager();
        auto ni = m_idx_proxy->GetNi();
        auto nj = m_idx_proxy->GetNj();
        auto nk = m_idx_proxy->GetNk();
        FlowSolverPara* para = GetPara();
        double prim_far[5];
        prim_far[0] = para->GetInflowDensity();
        prim_far[1] = para->GetInflowVelocityX();
        prim_far[2] = para->GetInflowVelocityY();
        prim_far[3] = para->GetInflowVelocityZ();
        prim_far[4] = para->GetInflowPressure();
        for (int iVal = 0; iVal < para->GetEquNum(); ++iVal)
        {
            for (int k = 0; k < nk; ++k)
            {
                for (int j = 0; j < nj; ++j)
                {
                    for (int i = 0; i < ni; ++i)
                    {
                        int idx = m_idx_proxy->GetIdx(i, j, k);
                        for (int iVal = 0; iVal < para->GetEquNum(); ++iVal)
                        {
                            data_manager->SetPrim(iVal, idx, prim_far[iVal]);
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
        FlowSolverPara* para = GetPara();
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
                    for (int iVal = 0; iVal < para->GetEquNum(); ++iVal)
                    {
                        data_manager->SetPrim(iVal, idx, prim_far[iVal]);
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
                    // data_manager->SetPrim(0, idx, pow(0.2 * i, 2) * pow(0.2 * j, 2.3));

                }
            }
        }

    }
    void NSSolverStruct::InitFieldBackup()
    {
    }
    void NSSolverStruct::CalcMetrics()
    {
        Log::info("Compute Struct Coordination Transformation Coefficients...");
        auto grid = GetGrid();
        auto data_manager = GetDataManager();
        auto ni = m_idx_proxy->GetNi();
        auto nj = m_idx_proxy->GetNj();
        auto nk = m_idx_proxy->GetNk();
        auto node = grid->GetNode();
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        const  double* xRight, * xLeft, * yRight, * yLeft, * zRight, * zLeft;
        xRight = xLeft = yRight = yLeft = zRight = zLeft = nullptr;
        double max_jacobian = -LARGE_NUMBER;
        double min_jacobian = LARGE_NUMBER;
        int max_jacobian_node = -1;
        int min_jacobian_node = -1;
        for (int i = is; i <= ie; ++i)
        {
            for (int j = js; j <= je; ++j)
            {
                for (int k = ks; k <= ke; ++k)
                {
                    int idx = m_idx_proxy->GetIdx(i, j, k);
                    xLeft = node->GetCoord(i - 1, j, k);
                    xRight = node->GetCoord(i + 1, j, k);
                    yLeft = node->GetCoord(i, j - 1, k);
                    yRight = node->GetCoord(i, j + 1, k);
                    if (grid->GetDim() == 3)
                    {
                        zLeft = node->GetCoord(i, j, k - 1);
                        zRight = node->GetCoord(i, j, k + 1);
                    }
                    m_node_metrics->CalcMetric(idx, xRight, xLeft, yRight, yLeft, zRight, zLeft);
                    if (m_node_metrics->GetJacobian(idx) > max_jacobian)
                    {
                        max_jacobian = m_node_metrics->GetJacobian(idx);
                        max_jacobian_node = idx;
                    }
                    if (m_node_metrics->GetJacobian(idx) < min_jacobian)
                    {
                        min_jacobian = m_node_metrics->GetJacobian(idx);
                        min_jacobian_node = idx;
                    }
                }
            }
        }
        Log::info("Compute Struct Coordination Transformation Coefficients Done.");
        int max_i, max_j, max_k;
        m_idx_proxy->GetIdxStruct(max_jacobian_node, max_i, max_j, max_k);
        Log::info("Max Jacobian = {0}, Node = ({}, {}, {})", max_jacobian, max_i, max_j, max_k);
        int min_i, min_j, min_k;
        m_idx_proxy->GetIdxStruct(min_jacobian_node, min_i, min_j, min_k);
        Log::info("Min Jacobian = {0}, Node = ({}, {}, {})", min_jacobian, min_i, min_j, min_k);
    }
    void NSSolverStruct::CalcMetricsCMM1()
    {
        Log::info("Compute Struct Coordination Transformation Coefficients...");
        auto grid = GetGrid();
        auto data_manager = GetDataManager();
        auto ni = m_idx_proxy->GetNi();
        auto nj = m_idx_proxy->GetNj();
        auto nk = m_idx_proxy->GetNk();
        //CMM1计算度量系数时需要的中间数据，逆变换度量系数与坐标乘积
        struct MetricsCMM1
        {
            double x_eta_y;
            double x_zeta_z;
            double y_eta_z;
            double y_zeta_z;
            double z_xi_x;
            double z_eta_y;
            double z_zeta_x;
            double x[4];//x_xi, x_eta, x_zeta, x_tau
            double y[4];//y_xi, y_eta, y_zeta, y_tau
        };
        //计算度量系数时的中间数据
        std::vector<std::vector<std::vector<MetricsCMM1>>> metrics_temp;
        auto node = grid->GetNode();
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        const  double* xRight, * xLeft, * yRight, * yLeft, * zRight, * zLeft;
        xRight = xLeft = yRight = yLeft = zRight = zLeft = nullptr;
        double max_jacobian = -LARGE_NUMBER;
        double min_jacobian = LARGE_NUMBER;
        int max_jacobian_node = -1;
        int min_jacobian_node = -1;

        for (int i = is; i <= ie; ++i)
        {
            for (int j = js; j <= je; ++j)
            {
                for (int k = ks; k <= ke; ++k)
                {


                }
            }
        }
        Log::info("Compute Struct Coordination Transformation Coefficients Done.");

    }
    void NSSolverStruct::CalcMetricsCMM2()
    {
    }
    void NSSolverStruct::CalcMetricsSCMM()
    {
    }
    void NSSolverStruct::CalcMetricsFSCMM()
    {
    }
    void NSSolverStruct::Preprocess()
    {
        NSSolver::Preprocess();
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
                    double u_xi = data_manager->GetVelocity(0, idx) * xi[0] + data_manager->GetVelocity(1, idx) * xi[1] + data_manager->GetVelocity(2, idx) * xi[2];
                    double u_eta = data_manager->GetVelocity(0, idx) * eta[0] + data_manager->GetVelocity(1, idx) * eta[1] + data_manager->GetVelocity(2, idx) * eta[2];
                    double u_zeta = data_manager->GetVelocity(0, idx) * zeta[0] + data_manager->GetVelocity(1, idx) * zeta[1] + data_manager->GetVelocity(2, idx) * zeta[2];
                    double lamda = abs(u_xi) + abs(u_eta) + c * (norm_xi + norm_eta );
                    if(grid->GetDim() == 3)
                    {
                        lamda += abs(u_zeta) + c * norm_zeta;
                    }
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
    void NSSolverStruct::ReduceTimeStep(double& dt)
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
        const DArray& rk_coef = para->GetRKCoef();
        int rkStage = rk_coef.size();
        double dt, jacobi;
        for (int iStage = 0; iStage < rkStage; ++iStage)
        {
            CalcResidual();
            for (int i = is;i <= ie;++i)
            {
                for (int j = js;j <= je;++j)
                {
                    for (int k = ks;k <= ke;++k)
                    {
                        int idx = m_idx_proxy->GetIdx(i, j, k);
                        dt = data_manager->GetTimeStep(idx);
                        jacobi = m_node_metrics->GetJacobian(idx);
                        for (int iVal = 0; iVal < 5; ++iVal)
                        {
                            data_manager->SetCons(iVal, idx, data_manager->GetCons(iVal, idx) + dt * rk_coef[iStage] * data_manager->GetResidual(iVal, idx) * jacobi);
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
                    for (int iVal = 0; iVal < GetPara()->GetEquNum(); ++iVal)
                    {
                        prim[iVal] = data_manager->GetPrim(iVal, idx);
                    }
                    gas->Prim2Cons(prim, cons);
                    for (int iEqu = 0; iEqu < GetPara()->GetEquNum(); ++iEqu)
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
                    for (int iVal = 0; iVal < GetPara()->GetEquNum(); ++iVal)
                    {
                        cons[iVal] = data_manager->GetCons(iVal, idx);
                    }
                    gas->Cons2Prim(cons, prim);
                    for (int iEqu = 0; iEqu < GetPara()->GetEquNum(); ++iEqu)
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
                    for (int iVal = 0; iVal < GetPara()->GetEquNum(); ++iVal)
                    {
                        data_manager->SetResidual(iVal, idx, 0.0);
                    }
                }
            }
        }
    }
    void NSSolverStruct::CalcMidNodePrim()
    {
        int interplot_type = 0;
        if (interplot_type == 0)
        {
            CalcMidNodePrimMUSCL();
        }
        else if (interplot_type == 1)
        {
            CalcMidNodePrimWCNS5();
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
        for (int k = 0; k < nk; ++k)
        {
            for (int j = 0; j < nj; ++j)
            {
                for (int i = 0; i < ni; ++i)
                {
                    int idx = m_idx_proxy->GetIdx(i, j, k);
                    int idx_left = m_idx_proxy->GetIdx(i - 1, j, k);
                    int idx_right = m_idx_proxy->GetIdx(i + 1, j, k);
                    for (int iVal = 0; iVal < 5; ++iVal)
                    {
                        data_manager->SetMidNodePrimLeft(iVal, 0, idx, data_manager->GetPrim(iVal, idx_left));
                        data_manager->SetMidNodePrimRight(iVal, 0, idx, data_manager->GetPrim(iVal, idx_right));
                    }
                    idx_left = m_idx_proxy->GetIdx(i, j - 1, k);
                    idx_right = m_idx_proxy->GetIdx(i, j + 1, k);
                    for (int iVal = 0; iVal < 5; ++iVal)
                    {
                        data_manager->SetMidNodePrimLeft(iVal, 1, idx, data_manager->GetPrim(iVal, idx_left));
                        data_manager->SetMidNodePrimRight(iVal, 1, idx, data_manager->GetPrim(iVal, idx_right));
                    }
                    if (grid->GetDim() == 3)
                    {
                        idx_left = m_idx_proxy->GetIdx(i, j, k - 1);
                        idx_right = m_idx_proxy->GetIdx(i, j, k + 1);
                        for (int iVal = 0; iVal < 5; ++iVal)
                        {
                            data_manager->SetMidNodePrimLeft(iVal, 2, idx, data_manager->GetPrim(iVal, idx_left));
                            data_manager->SetMidNodePrimRight(iVal, 2, idx, data_manager->GetPrim(iVal, idx_right));
                        }
                    }
                }
            }
        }
    }
    void NSSolverStruct::CalcMidNodePrimGrad()
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
        for (int k = 0; k < nk; ++k)
        {
            for (int j = 0; j < nj; ++j)
            {
                for (int i = 0; i < ni; ++i)
                {
                    int idx = m_idx_proxy->GetIdx(i, j, k);
                    int idx_left = m_idx_proxy->GetIdx(i - 1, j, k);
                    int idx_right = m_idx_proxy->GetIdx(i + 1, j, k);
                    coord_vec[0] = node->GetCoord(i + 1, j, k)[0] - node->GetCoord(i, j, k)[0];
                    coord_vec[1] = node->GetCoord(i + 1, j, k)[1] - node->GetCoord(i, j, k)[1];
                    coord_vec[2] = node->GetCoord(i + 1, j, k)[2] - node->GetCoord(i, j, k)[2];
                    MidNodeGrad(idx, idx_right, coord_vec, value_left, value_right);
                    for (int iVal = 0; iVal < 5; ++iVal)
                    {
                        data_manager->SetMidNodePrim(iVal, 0, idx, value_left[iVal], value_right[iVal]);
                    }
                    coord_vec[0] = node->GetCoord(i, j, k)[0] - node->GetCoord(i - 1, j, k)[0];
                    coord_vec[1] = node->GetCoord(i, j, k)[1] - node->GetCoord(i - 1, j, k)[1];
                    coord_vec[2] = node->GetCoord(i, j, k)[2] - node->GetCoord(i - 1, j, k)[2];
                    MidNodeGrad(idx_left, idx, coord_vec, value_left, value_right);
                    for (int iVal = 0; iVal < 5; ++iVal)
                    {
                        data_manager->SetMidNodePrim(iVal, 0, idx, value_left[iVal], value_right[iVal]);
                    }
                    idx_left = m_idx_proxy->GetIdx(i, j - 1, k);
                    idx_right = m_idx_proxy->GetIdx(i, j + 1, k);
                    coord_vec[0] = node->GetCoord(i, j + 1, k)[0] - node->GetCoord(i, j, k)[0];
                    coord_vec[1] = node->GetCoord(i, j + 1, k)[1] - node->GetCoord(i, j, k)[1];
                    coord_vec[2] = node->GetCoord(i, j + 1, k)[2] - node->GetCoord(i, j, k)[2];
                    MidNodeGrad(idx, idx_right, coord_vec, value_left, value_right);
                    for (int iVal = 0; iVal < 5; ++iVal)
                    {
                        data_manager->SetMidNodePrim(iVal, 1, idx, value_left[iVal], value_right[iVal]);
                    }
                    coord_vec[0] = node->GetCoord(i, j, k)[0] - node->GetCoord(i, j - 1, k)[0];
                    coord_vec[1] = node->GetCoord(i, j, k)[1] - node->GetCoord(i, j - 1, k)[1];
                    coord_vec[2] = node->GetCoord(i, j, k)[2] - node->GetCoord(i, j - 1, k)[2];
                    MidNodeGrad(idx_left, idx, coord_vec, value_left, value_right);
                    for (int iVal = 0; iVal < 5; ++iVal)
                    {
                        data_manager->SetMidNodePrim(iVal, 1, idx, value_left[iVal], value_right[iVal]);
                    }
                    if (grid->GetDim() == 3)
                    {
                        idx_left = m_idx_proxy->GetIdx(i, j, k - 1);
                        idx_right = m_idx_proxy->GetIdx(i, j, k + 1);
                        coord_vec[0] = node->GetCoord(i, j, k + 1)[0] - node->GetCoord(i, j, k)[0];
                        coord_vec[1] = node->GetCoord(i, j, k + 1)[1] - node->GetCoord(i, j, k)[1];
                        coord_vec[2] = node->GetCoord(i, j, k + 1)[2] - node->GetCoord(i, j, k)[2];
                        MidNodeGrad(idx, idx_right, coord_vec, value_left, value_right);
                        for (int iVal = 0; iVal < 5; ++iVal)
                        {
                            data_manager->SetMidNodePrim(iVal, 2, idx, value_left[iVal], value_right[iVal]);
                        }
                        coord_vec[0] = node->GetCoord(i, j, k)[0] - node->GetCoord(i, j, k - 1)[0];
                        coord_vec[1] = node->GetCoord(i, j, k)[1] - node->GetCoord(i, j, k - 1)[1];
                        coord_vec[2] = node->GetCoord(i, j, k)[2] - node->GetCoord(i, j, k - 1)[2];
                        MidNodeGrad(idx_left, idx, coord_vec, value_left, value_right);
                        for (int iVal = 0; iVal < 5; ++iVal)
                        {
                            data_manager->SetMidNodePrim(iVal, 2, idx, value_left[iVal], value_right[iVal]);
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
        int equ_num = GetPara()->GetEquNum();
        auto ni = m_idx_proxy->GetNi();
        auto nj = m_idx_proxy->GetNj();
        auto nk = m_idx_proxy->GetNk();
        //MUSCL整点变量(i-2,i-1,i,i+1,i+2)的值
        double value_temp[5];
        int idx_temp[5];
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        double left_value, right_value;
        ///@note 从ks-1开始，到ke结束，因为对第一个计算点ks进行通量差分时，需要使用ks-1/2处的值
        for (int k = ks;k <= ke;++k)
        {
            for (int j = js;j <= je;++j)
            {
                for (int i = is;i <= ie;++i)
                {
                    idx_temp[0] = m_idx_proxy->GetIdx(i - 2, j, k);
                    idx_temp[1] = m_idx_proxy->GetIdx(i - 1, j, k);
                    idx_temp[2] = m_idx_proxy->GetIdx(i, j, k);
                    idx_temp[3] = m_idx_proxy->GetIdx(i + 1, j, k);
                    idx_temp[4] = m_idx_proxy->GetIdx(i + 2, j, k);
                    for (int iVal = 0; iVal < equ_num; ++iVal)
                    {
                        value_temp[0] = data_manager->GetPrim(iVal, idx_temp[0]);
                        value_temp[1] = data_manager->GetPrim(iVal, idx_temp[1]);
                        value_temp[2] = data_manager->GetPrim(iVal, idx_temp[2]);
                        value_temp[3] = data_manager->GetPrim(iVal, idx_temp[3]);
                        value_temp[4] = data_manager->GetPrim(iVal, idx_temp[4]);
                        MidNodeMUSCL(value_temp, left_value, right_value);
                        data_manager->SetMidNodePrim(iVal, 0, idx_temp[2], left_value, right_value);
                    }
                    idx_temp[0] = m_idx_proxy->GetIdx(i, j - 2, k);
                    idx_temp[1] = m_idx_proxy->GetIdx(i, j - 1, k);
                    idx_temp[3] = m_idx_proxy->GetIdx(i, j + 1, k);
                    idx_temp[4] = m_idx_proxy->GetIdx(i, j + 2, k);
                    for (int iVal = 0; iVal < equ_num; ++iVal)
                    {
                        value_temp[0] = data_manager->GetPrim(iVal, idx_temp[0]);
                        value_temp[1] = data_manager->GetPrim(iVal, idx_temp[1]);
                        value_temp[2] = data_manager->GetPrim(iVal, idx_temp[2]);
                        value_temp[3] = data_manager->GetPrim(iVal, idx_temp[3]);
                        value_temp[4] = data_manager->GetPrim(iVal, idx_temp[4]);
                        MidNodeMUSCL(value_temp, left_value, right_value);
                        data_manager->SetMidNodePrim(iVal, 1, idx_temp[2], left_value, right_value);
                    }
                    if (grid->GetDim() == 3)
                    {
                        idx_temp[0] = m_idx_proxy->GetIdx(i, j, k - 2);
                        idx_temp[1] = m_idx_proxy->GetIdx(i, j, k - 1);
                        idx_temp[3] = m_idx_proxy->GetIdx(i, j, k + 1);
                        idx_temp[4] = m_idx_proxy->GetIdx(i, j, k + 2);
                        for (int iVal = 0; iVal < equ_num; ++iVal)
                        {
                            value_temp[0] = data_manager->GetPrim(iVal, idx_temp[0]);
                            value_temp[1] = data_manager->GetPrim(iVal, idx_temp[1]);
                            value_temp[2] = data_manager->GetPrim(iVal, idx_temp[2]);
                            value_temp[3] = data_manager->GetPrim(iVal, idx_temp[3]);
                            value_temp[4] = data_manager->GetPrim(iVal, idx_temp[4]);
                            MidNodeMUSCL(value_temp, left_value, right_value);
                            data_manager->SetMidNodePrim(iVal, 2, idx_temp[2], left_value, right_value);
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
        int equ_num = GetPara()->GetEquNum();
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
        for (k = ks;k <= ke;++k)
        {
            for (j = js;j <= je;++j)
            {
                //i=is-1处的值
                //不存在i-2处的值且不参与计算,用i-1处的值代替，不影响计算结果
                idx_temp[0] = m_idx_proxy->GetIdx(is - 2, j, k);
                idx_temp[1] = m_idx_proxy->GetIdx(is - 2, j, k);
                idx_temp[2] = m_idx_proxy->GetIdx(is - 1, j, k);
                idx_temp[3] = m_idx_proxy->GetIdx(is, j, k);
                idx_temp[4] = m_idx_proxy->GetIdx(is + 1, j, k);
                for (int iVal = 0; iVal < equ_num; ++iVal)
                {
                    value[0] = data_manager->GetPrim(iVal, idx_temp[0]);
                    value[1] = data_manager->GetPrim(iVal, idx_temp[1]);
                    value[2] = data_manager->GetPrim(iVal, idx_temp[2]);
                    value[3] = data_manager->GetPrim(iVal, idx_temp[3]);
                    value[4] = data_manager->GetPrim(iVal, idx_temp[4]);
                    MidNodeMUSCL(value, left_value, right_value);
                    data_manager->SetMidNodePrim(iVal, 0, idx_temp[2], left_value, right_value);
                }
            }
        }

        for (k = ks;k <= ke;++k)
        {
            for (i = is;i <= ie;++i)
            {
                //j=js-1处的值，(js-1/2)右值会用到，但(js-1/2)左值在计算中不会使用
                //不存在j-2处的值,用j-1处的值代替，不影响计算结果
                idx_temp[0] = m_idx_proxy->GetIdx(i, js - 2, k);
                idx_temp[1] = m_idx_proxy->GetIdx(i, js - 2, k);
                idx_temp[2] = m_idx_proxy->GetIdx(i, js - 1, k);
                idx_temp[3] = m_idx_proxy->GetIdx(i, js, k);
                idx_temp[4] = m_idx_proxy->GetIdx(i, js + 1, k);
                for (int iVal = 0; iVal < equ_num; ++iVal)
                {
                    value[0] = data_manager->GetPrim(iVal, idx_temp[0]);
                    value[1] = data_manager->GetPrim(iVal, idx_temp[1]);
                    value[2] = data_manager->GetPrim(iVal, idx_temp[2]);
                    value[3] = data_manager->GetPrim(iVal, idx_temp[3]);
                    value[4] = data_manager->GetPrim(iVal, idx_temp[4]);
                    MidNodeMUSCL(value, left_value, right_value);
                    data_manager->SetMidNodePrim(iVal, 1, idx_temp[2], left_value, right_value);
                }
            }
        }

        if (grid->GetDim() == 3)
        {
            for (j = js;j <= je;++j)
            {
                for (i = is;i <= ie;++i)
                {
                    //k=ks-1处的值，(ks-1/2)右值会用到，但(ks-1/2)左值在计算中不会使用
                    //不存在k-2处的值,用k-1处的值代替，不影响计算结果
                    idx_temp[0] = m_idx_proxy->GetIdx(i, j, ks - 2);
                    idx_temp[1] = m_idx_proxy->GetIdx(i, j, ks - 2);
                    idx_temp[2] = m_idx_proxy->GetIdx(i, j, ks - 1);
                    idx_temp[3] = m_idx_proxy->GetIdx(i, j, ks);
                    idx_temp[4] = m_idx_proxy->GetIdx(i, j, ks + 1);
                    for (int iVal = 0; iVal < equ_num; ++iVal)
                    {
                        value[0] = data_manager->GetPrim(iVal, idx_temp[0]);
                        value[1] = data_manager->GetPrim(iVal, idx_temp[1]);
                        value[2] = data_manager->GetPrim(iVal, idx_temp[2]);
                        value[3] = data_manager->GetPrim(iVal, idx_temp[3]);
                        value[4] = data_manager->GetPrim(iVal, idx_temp[4]);
                        MidNodeMUSCL(value, left_value, right_value);
                        data_manager->SetMidNodePrim(iVal, 2, idx_temp[2], left_value, right_value);
                    }
                }
            }
        }
    }

    void NSSolverStruct::CalcMidNodePrimWCNS5()
    {
        auto grid = GetGrid();
        auto data_manager = GetDataManager();
        int equ_num = GetPara()->GetEquNum();
        auto ni = m_idx_proxy->GetNi();
        auto nj = m_idx_proxy->GetNj();
        auto nk = m_idx_proxy->GetNk();
        //WCNS整点变量(i-2,i-1,i,i+1,i+2)的值
        double value_temp[5];
        int idx_temp[5];
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        double left_value, right_value;
        ///@note 从ks-1开始，到ke结束，因为对第一个计算点ks进行通量差分时，需要使用ks-1/2处的值
        for (int k = ks - 1;k <= ke + 1;++k)
        {
            for (int j = js - 1;j <= je + 1;++j)
            {
                for (int i = is - 1;i <= ie + 1;++i)
                {
                    idx_temp[0] = m_idx_proxy->GetIdx(i - 2, j, k);
                    idx_temp[1] = m_idx_proxy->GetIdx(i - 1, j, k);
                    idx_temp[2] = m_idx_proxy->GetIdx(i, j, k);
                    idx_temp[3] = m_idx_proxy->GetIdx(i + 1, j, k);
                    idx_temp[4] = m_idx_proxy->GetIdx(i + 2, j, k);
                    for (int iVal = 0; iVal < equ_num; ++iVal)
                    {
                        value_temp[0] = data_manager->GetPrim(iVal, idx_temp[0]);
                        value_temp[1] = data_manager->GetPrim(iVal, idx_temp[1]);
                        value_temp[2] = data_manager->GetPrim(iVal, idx_temp[2]);
                        value_temp[3] = data_manager->GetPrim(iVal, idx_temp[3]);
                        value_temp[4] = data_manager->GetPrim(iVal, idx_temp[4]);
                        MidNodeWCNS5(value_temp, left_value, right_value);
                        data_manager->SetMidNodePrimLeft(iVal, 0, idx_temp[2], left_value);
                        data_manager->SetMidNodePrimRight(iVal, 0, idx_temp[1], right_value);
                    }
                    idx_temp[0] = m_idx_proxy->GetIdx(i, j - 2, k);
                    idx_temp[1] = m_idx_proxy->GetIdx(i, j - 1, k);
                    idx_temp[3] = m_idx_proxy->GetIdx(i, j + 1, k);
                    idx_temp[4] = m_idx_proxy->GetIdx(i, j + 2, k);
                    for (int iVal = 0; iVal < equ_num; ++iVal)
                    {
                        value_temp[0] = data_manager->GetPrim(iVal, idx_temp[0]);
                        value_temp[1] = data_manager->GetPrim(iVal, idx_temp[1]);
                        value_temp[2] = data_manager->GetPrim(iVal, idx_temp[2]);
                        value_temp[3] = data_manager->GetPrim(iVal, idx_temp[3]);
                        value_temp[4] = data_manager->GetPrim(iVal, idx_temp[4]);
                        MidNodeWCNS5(value_temp, left_value, right_value);
                        data_manager->SetMidNodePrimLeft(iVal, 1, idx_temp[2], left_value);
                        data_manager->SetMidNodePrimRight(iVal, 1, idx_temp[1], right_value);
                    }
                    if (grid->GetDim() == 3)
                    {
                        idx_temp[0] = m_idx_proxy->GetIdx(i, j, k - 2);
                        idx_temp[1] = m_idx_proxy->GetIdx(i, j, k - 1);
                        idx_temp[3] = m_idx_proxy->GetIdx(i, j, k + 1);
                        idx_temp[4] = m_idx_proxy->GetIdx(i, j, k + 2);
                        for (int iVal = 0; iVal < equ_num; ++iVal)
                        {
                            value_temp[0] = data_manager->GetPrim(iVal, idx_temp[0]);
                            value_temp[1] = data_manager->GetPrim(iVal, idx_temp[1]);
                            value_temp[2] = data_manager->GetPrim(iVal, idx_temp[2]);
                            value_temp[3] = data_manager->GetPrim(iVal, idx_temp[3]);
                            value_temp[4] = data_manager->GetPrim(iVal, idx_temp[4]);
                            MidNodeWCNS5(value_temp, left_value, right_value);
                            data_manager->SetMidNodePrimLeft(iVal, 2, idx_temp[2], left_value);
                            data_manager->SetMidNodePrimRight(iVal, 2, idx_temp[1], right_value);
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
        int equ_num = GetPara()->GetEquNum();
        auto ni = m_idx_proxy->GetNi();
        auto nj = m_idx_proxy->GetNj();
        auto nk = m_idx_proxy->GetNk();
        double value_left[4], value_right[4];
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        //左侧的模板点(0,1,2,3),右侧的模板点(N-4,N-3,N-2,N-1)
        int idx_temp_left[4], idx_temp_right[4];
        double left_value, right_value;
        int i, j, k;
        for (k = ks;k <= ke;++k)
        {
            for (j = js;j <= je;++j)
            {
                idx_temp_left[0] = m_idx_proxy->GetIdx(0, j, k);
                idx_temp_left[1] = m_idx_proxy->GetIdx(1, j, k);
                idx_temp_left[2] = m_idx_proxy->GetIdx(2, j, k);
                idx_temp_left[3] = m_idx_proxy->GetIdx(3, j, k);
                idx_temp_right[0] = m_idx_proxy->GetIdx(ni - 4, j, k);
                idx_temp_right[1] = m_idx_proxy->GetIdx(ni - 3, j, k);
                idx_temp_right[2] = m_idx_proxy->GetIdx(ni - 2, j, k);
                idx_temp_right[3] = m_idx_proxy->GetIdx(ni - 1, j, k);
                for (int iVal = 0; iVal < equ_num; ++iVal)
                {
                    value_left[0] = data_manager->GetPrim(iVal, idx_temp_left[0]);
                    value_left[1] = data_manager->GetPrim(iVal, idx_temp_left[1]);
                    value_left[2] = data_manager->GetPrim(iVal, idx_temp_left[2]);
                    value_left[3] = data_manager->GetPrim(iVal, idx_temp_left[3]);
                    value_right[0] = data_manager->GetPrim(iVal, idx_temp_right[0]);
                    value_right[1] = data_manager->GetPrim(iVal, idx_temp_right[1]);
                    value_right[2] = data_manager->GetPrim(iVal, idx_temp_right[2]);
                    value_right[3] = data_manager->GetPrim(iVal, idx_temp_right[3]);
                    left_value = right_value = (5 * value_left[0] + 15 * value_left[1] - 5 * value_left[2] + value_left[3]) / 16;
                    data_manager->SetMidNodePrim(iVal, 0, idx_temp_left[0], left_value, right_value);
                    left_value = right_value = (-value_left[0] + 9 * value_left[1] + 9 * value_left[2] - value_left[3]) / 16;
                    data_manager->SetMidNodePrim(iVal, 0, idx_temp_left[1], left_value, right_value);
                    left_value = right_value = (5 * value_right[0] + 15 * value_right[1] - 5 * value_right[2] + value_right[3]) / 16;
                    data_manager->SetMidNodePrim(iVal, 0, idx_temp_right[1], left_value, right_value);
                    left_value = right_value = (-value_right[0] + 9 * value_right[1] + 9 * value_right[2] - value_right[3]) / 16;
                    data_manager->SetMidNodePrim(iVal, 0, idx_temp_right[2], left_value, right_value);
                }
            }
        }

        for (k = ks;k <= ke;++k)
        {
            for (i = is;i <= ie;++i)
            {
                //j=js-1处的值，(js-1/2)右值会用到，但(js-1/2)左值在计算中不会使用
                //不存在j-2处的值,用j-1处的值代替，不影响计算结果
                idx_temp_left[0] = m_idx_proxy->GetIdx(i, 0, k);
                idx_temp_left[1] = m_idx_proxy->GetIdx(i, 1, k);
                idx_temp_left[2] = m_idx_proxy->GetIdx(i, 2, k);
                idx_temp_left[3] = m_idx_proxy->GetIdx(i, 3, k);
                idx_temp_right[0] = m_idx_proxy->GetIdx(i, nj - 4, k);
                idx_temp_right[1] = m_idx_proxy->GetIdx(i, nj - 3, k);
                idx_temp_right[2] = m_idx_proxy->GetIdx(i, nj - 2, k);
                idx_temp_right[3] = m_idx_proxy->GetIdx(i, nj - 1, k);
                for (int iVal = 0; iVal < equ_num; ++iVal)
                {
                    value_left[0] = data_manager->GetPrim(iVal, idx_temp_left[0]);
                    value_left[1] = data_manager->GetPrim(iVal, idx_temp_left[1]);
                    value_left[2] = data_manager->GetPrim(iVal, idx_temp_left[2]);
                    value_left[3] = data_manager->GetPrim(iVal, idx_temp_left[3]);
                    value_right[0] = data_manager->GetPrim(iVal, idx_temp_right[0]);
                    value_right[1] = data_manager->GetPrim(iVal, idx_temp_right[1]);
                    value_right[2] = data_manager->GetPrim(iVal, idx_temp_right[2]);
                    value_right[3] = data_manager->GetPrim(iVal, idx_temp_right[3]);
                    left_value = right_value = (5 * value_left[0] + 15 * value_left[1] - 5 * value_left[2] + value_left[3]) / 16;
                    data_manager->SetMidNodePrim(iVal, 1, idx_temp_left[0], left_value, right_value);
                    left_value = right_value = (-value_left[0] + 9 * value_left[1] + 9 * value_left[2] - value_left[3]) / 16;
                    data_manager->SetMidNodePrim(iVal, 1, idx_temp_left[1], left_value, right_value);
                    left_value = right_value = (5 * value_right[0] + 15 * value_right[1] - 5 * value_right[2] + value_right[3]) / 16;
                    data_manager->SetMidNodePrim(iVal, 1, idx_temp_right[1], left_value, right_value);
                    left_value = right_value = (-value_right[0] + 9 * value_right[1] + 9 * value_right[2] - value_right[3]) / 16;
                    data_manager->SetMidNodePrim(iVal, 1, idx_temp_right[2], left_value, right_value);

                }
            }
        }

        if (grid->GetDim() == 3)
        {
            for (j = js;j <= je;++j)
            {
                for (i = is;i <= ie;++i)
                {
                    //k=ks-1处的值，(ks-1/2)右值会用到，但(ks-1/2)左值在计算中不会使用
                    //不存在k-2处的值,用k-1处的值代替，不影响计算结果
                    idx_temp_left[0] = m_idx_proxy->GetIdx(i, j, 0);
                    idx_temp_left[1] = m_idx_proxy->GetIdx(i, j, 1);
                    idx_temp_left[2] = m_idx_proxy->GetIdx(i, j, 2);
                    idx_temp_left[3] = m_idx_proxy->GetIdx(i, j, 3);
                    idx_temp_right[0] = m_idx_proxy->GetIdx(i, j, nk - 4);
                    idx_temp_right[1] = m_idx_proxy->GetIdx(i, j, nk - 3);
                    idx_temp_right[2] = m_idx_proxy->GetIdx(i, j, nk - 2);
                    idx_temp_right[3] = m_idx_proxy->GetIdx(i, j, nk - 1);
                    for (int iVal = 0; iVal < equ_num; ++iVal)
                    {
                        value_left[0] = data_manager->GetPrim(iVal, idx_temp_left[0]);
                        value_left[1] = data_manager->GetPrim(iVal, idx_temp_left[1]);
                        value_left[2] = data_manager->GetPrim(iVal, idx_temp_left[2]);
                        value_left[3] = data_manager->GetPrim(iVal, idx_temp_left[3]);
                        value_right[0] = data_manager->GetPrim(iVal, idx_temp_right[0]);
                        value_right[1] = data_manager->GetPrim(iVal, idx_temp_right[1]);
                        value_right[2] = data_manager->GetPrim(iVal, idx_temp_right[2]);
                        value_right[3] = data_manager->GetPrim(iVal, idx_temp_right[3]);
                        left_value = right_value = (5 * value_left[0] + 15 * value_left[1] - 5 * value_left[2] + value_left[3]) / 16;
                        data_manager->SetMidNodePrim(iVal, 2, idx_temp_left[0], left_value, right_value);
                        left_value = right_value = (-value_left[0] + 9 * value_left[1] + 9 * value_left[2] - value_left[3]) / 16;
                        data_manager->SetMidNodePrim(iVal, 2, idx_temp_left[1], left_value, right_value);
                        left_value = right_value = (5 * value_right[0] + 15 * value_right[1] - 5 * value_right[2] + value_right[3]) / 16;
                        data_manager->SetMidNodePrim(iVal, 2, idx_temp_right[1], left_value, right_value);
                        left_value = right_value = (-value_right[0] + 9 * value_right[1] + 9 * value_right[2] - value_right[3]) / 16;
                        data_manager->SetMidNodePrim(iVal, 2, idx_temp_right[2], left_value, right_value);
                    }
                }
            }
        }
    }

    void NSSolverStruct::MidNodeGrad(int idx_left, int idx_right, double* coord_vec, double* value_left, double* value_right)
    {
        int equ_num = GetPara()->GetEquNum();
        auto data_manager = GetDataManager();
        for (int iVal = 0; iVal < equ_num; ++iVal)
        {
            value_left[iVal] = data_manager->GetPrim(iVal, idx_left) + 0.5 * data_manager->GetLimiter(iVal, idx_left) *
                (coord_vec[0] * data_manager->GetPrimGrad(iVal, 0, idx_left) +
                    coord_vec[1] * data_manager->GetPrimGrad(iVal, 1, idx_left) +
                    coord_vec[2] * data_manager->GetPrimGrad(iVal, 2, idx_left));
            value_right[iVal] = data_manager->GetPrim(iVal, idx_right) - 0.5 * data_manager->GetLimiter(iVal, idx_right) *
                (coord_vec[0] * data_manager->GetPrimGrad(iVal, 0, idx_right) +
                    coord_vec[1] * data_manager->GetPrimGrad(iVal, 1, idx_right) +
                    coord_vec[2] * data_manager->GetPrimGrad(iVal, 2, idx_right));
        }
    }
    void NSSolverStruct::MidNode1st(int index_left, int index_right, double* value_rec_left, double* value_rec_right)
    {
        int equ_num = GetPara()->GetEquNum();
        auto  grid = GetGrid();
        auto data_manager = GetDataManager();
        for (int iVal = 0; iVal < equ_num; ++iVal)
        {
            value_rec_left[iVal] = data_manager->GetPrim(iVal, index_left);
            value_rec_right[iVal] = data_manager->GetPrim(iVal, index_right);
        }
    }

    void NSSolverStruct::MidNodeMUSCL(const double* value, double& value_left, double& value_right)
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
            Log::error("value[0] = {}, value[1] = {}, value[2] = {}, value[3] = {}, value[4] = {}", value[0], value[1], value[2], value[3], value[4]);
            exit(0);
        }
    }

    //参考《计算空气动力学》pp.84
    void NSSolverStruct::MidNodeWCNS5(const double* value, double& value_left, double& value_right)
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
        for (int i = 0;i < 3;++i)
        {
            IS[i] = pow(f[i] * h, 2) + pow(s[i] * h * h, 2);
        }
        double eps = 1.0e-6;
        for (int i = 0; i < 3; ++i)
        {
            beta_left[i] = c_left[i] / pow(IS[i] + eps, 2);
            beta_right[i] = c_right[i] / pow(IS[i] + eps, 2);
        }
        for (int i = 0;i < 3;++i)
        {
            omega_left[i] = beta_left[i] / (beta_left[0] + beta_left[1] + beta_left[2]);
            omega_right[i] = beta_right[i] / (beta_right[0] + beta_right[1] + beta_right[2]);
        }
        f_left = f_right = s_left = s_right = 0.0;
        for (int i = 0;i < 3;++i)
        {
            f_left += omega_left[i] * f[i];
            s_left += omega_left[i] * s[i];
            f_right += omega_right[i] * f[i];
            s_right += omega_right[i] * s[i];
        }
        value_left = value[2] + 0.5 * h * f_left + 0.125 * h * h * s_left;
        value_right = value[2] - 0.5 * h * f_right + 0.125 * h * h * s_right;
    }
    void NSSolverStruct::FluxDifferenceMUSCL()
    {
        auto grid = GetGrid();
        auto para = GetPara();
        auto data_manager = GetDataManager();
        auto gas = GetGas();
        auto equ_num = para->GetEquNum();
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        RiemannSolverPara riemann_para[6];
        for (int i = 0; i < 6; ++i)
        {
            riemann_para[i].gamma_left = riemann_para[i].gamma_right = gas->GetGamma();
        }
        //当前节点的编号
        int idx;
        //差分模板的编号
        int idx_temp[5];
        //差分模板的值
        double value[5];
        double res_tmp[5];
        for (int i = is;i <= ie;++i)
        {
            for (int j = js;j <= je;++j)
            {
                for (int k = ks;k <= ke;++k)
                {
                    for (int iVal = 0;iVal < equ_num;++iVal)
                    {
                        res_tmp[iVal] = data_manager->GetResidual(iVal, m_idx_proxy->GetIdx(i, j, k));
                    }
                    // i direction
                    idx = m_idx_proxy->GetIdx(i, j, k);
                    double jacobi = m_node_metrics->GetJacobian(idx);
                    riemann_para[0].norm(0) = m_node_metrics->GetXi(idx)[0];
                    riemann_para[0].norm(1) = m_node_metrics->GetXi(idx)[1];
                    riemann_para[0].norm(2) = m_node_metrics->GetXi(idx)[2];
                    riemann_para[0].nt = m_node_metrics->GetXi(idx)[3];
                    riemann_para[1].norm(0) = m_node_metrics->GetXi(idx)[0];
                    riemann_para[1].norm(1) = m_node_metrics->GetXi(idx)[1];
                    riemann_para[1].norm(2) = m_node_metrics->GetXi(idx)[2];
                    riemann_para[1].nt = m_node_metrics->GetXi(idx)[3];
                    for (int iVal = 0; iVal < equ_num; ++iVal)
                    {
                        riemann_para[0].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, 0, m_idx_proxy->GetIdx(i, j, k));
                        riemann_para[0].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, 0, m_idx_proxy->GetIdx(i, j, k));
                        riemann_para[1].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, 0, m_idx_proxy->GetIdx(i - 1, j, k));
                        riemann_para[1].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, 0, m_idx_proxy->GetIdx(i - 1, j, k));
                    }
                    m_riemann_solver->Solver(riemann_para[0]);
                    m_riemann_solver->Solver(riemann_para[1]);
                    for (int iVal = 0;iVal < equ_num;++iVal)
                    {
                        res_tmp[iVal] -= (riemann_para[0].flux[iVal] - riemann_para[1].flux[iVal]) / jacobi;
                    }
                    // j direction
                    riemann_para[2].norm(0) = m_node_metrics->GetEta(idx)[0];
                    riemann_para[2].norm(1) = m_node_metrics->GetEta(idx)[1];
                    riemann_para[2].norm(2) = m_node_metrics->GetEta(idx)[2];
                    riemann_para[2].nt = m_node_metrics->GetEta(idx)[3];
                    riemann_para[3].norm(0) = m_node_metrics->GetEta(idx)[0];
                    riemann_para[3].norm(1) = m_node_metrics->GetEta(idx)[1];
                    riemann_para[3].norm(2) = m_node_metrics->GetEta(idx)[2];
                    riemann_para[3].nt = m_node_metrics->GetEta(idx)[3];
                    for (int iVal = 0; iVal < equ_num; ++iVal)
                    {
                        riemann_para[2].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, 1, m_idx_proxy->GetIdx(i, j, k));
                        riemann_para[2].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, 1, m_idx_proxy->GetIdx(i, j, k));
                        riemann_para[3].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, 1, m_idx_proxy->GetIdx(i, j - 1, k));
                        riemann_para[3].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, 1, m_idx_proxy->GetIdx(i, j - 1, k));
                    }
                    m_riemann_solver->Solver(riemann_para[2]);
                    m_riemann_solver->Solver(riemann_para[3]);
                    for (int iVal = 0; iVal < equ_num; ++iVal)
                    {
                        res_tmp[iVal] -= (riemann_para[2].flux[iVal] - riemann_para[3].flux[iVal]) / jacobi;
                    }

                    if (grid->GetDim() == 3)
                    {
                        // k direction
                        riemann_para[4].norm(0) = m_node_metrics->GetZeta(idx)[0];
                        riemann_para[4].norm(1) = m_node_metrics->GetZeta(idx)[1];
                        riemann_para[4].norm(2) = m_node_metrics->GetZeta(idx)[2];
                        riemann_para[4].nt = m_node_metrics->GetZeta(idx)[3];
                        riemann_para[5].norm(0) = m_node_metrics->GetZeta(idx)[0];
                        riemann_para[5].norm(1) = m_node_metrics->GetZeta(idx)[1];
                        riemann_para[5].norm(2) = m_node_metrics->GetZeta(idx)[2];
                        riemann_para[5].nt = m_node_metrics->GetZeta(idx)[3];
                        idx_temp[0] = m_idx_proxy->GetIdx(i, j, k - 2);
                        idx_temp[1] = m_idx_proxy->GetIdx(i, j, k - 1);
                        idx_temp[2] = m_idx_proxy->GetIdx(i, j, k);
                        idx_temp[3] = m_idx_proxy->GetIdx(i, j, k + 1);
                        idx_temp[4] = m_idx_proxy->GetIdx(i, j, k + 2);
                        for (int iVal = 0; iVal < equ_num; ++iVal)
                        {
                            riemann_para[4].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, 2, m_idx_proxy->GetIdx(i, j, k));
                            riemann_para[4].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, 2, m_idx_proxy->GetIdx(i, j, k));
                            riemann_para[5].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, 2, m_idx_proxy->GetIdx(i, j, k - 1));
                            riemann_para[5].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, 2, m_idx_proxy->GetIdx(i, j, k - 1));
                        }
                        m_riemann_solver->Solver(riemann_para[4]);
                        m_riemann_solver->Solver(riemann_para[5]);
                        for (int iVal = 0; iVal < equ_num; ++iVal)
                        {
                            res_tmp[iVal] -= (riemann_para[4].flux[iVal] - riemann_para[5].flux[iVal]) / jacobi;
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
    void NSSolverStruct::FluxDifferenceWCNS5()
    {
        auto grid = GetGrid();
        auto para = GetPara();
        auto node = grid->GetNode();
        auto data_manager = GetDataManager();
        auto gas = GetGas();
        auto equ_num = para->GetEquNum();
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        RiemannSolverPara riemann_para[6];
        for (int i = 0; i < 6; ++i)
        {
            riemann_para[i].gamma_left = riemann_para[i].gamma_right = gas->GetGamma();
        }
        //当前节点的编号
        int idx;
        //差分模板的编号
        int idx_temp[5];
        //差分模板的值
        double value[5];
        double res_tmp[5];
        for (int i = is;i <= ie;++i)
        {
            for (int j = js;j <= je;++j)
            {
                for (int k = ks;k <= ke;++k)
                {
                    for (int iVal = 0; iVal < equ_num; ++iVal)
                    {
                        res_tmp[iVal] = data_manager->GetResidual(iVal, m_idx_proxy->GetIdx(i, j, k));
                    }
                    // i direction
                    idx = m_idx_proxy->GetIdx(i, j, k);
                    double jacobi = m_node_metrics->GetJacobian(idx);
                    for (int i = 0;i < 6;++i)
                    {
                        riemann_para[i].norm(0) = m_node_metrics->GetXi(idx)[0];
                        riemann_para[i].norm(1) = m_node_metrics->GetXi(idx)[1];
                        riemann_para[i].norm(2) = m_node_metrics->GetXi(idx)[2];
                        riemann_para[i].nt = m_node_metrics->GetXi(idx)[3];
                    }
                    // Log::info("(i,j,k) = ({},{},{}), Coord: ({},{},{})", i, j, k, node->GetCoord(i, j, k)[0], node->GetCoord(i, j, k)[1], node->GetCoord(i, j, k)[2]);
                    for (int iVal = 0; iVal < equ_num; ++iVal)
                    {
                        riemann_para[0].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, 0, m_idx_proxy->GetIdx(i - 3, j, k));
                        riemann_para[0].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, 0, m_idx_proxy->GetIdx(i - 3, j, k));
                        riemann_para[1].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, 0, m_idx_proxy->GetIdx(i - 2, j, k));
                        riemann_para[1].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, 0, m_idx_proxy->GetIdx(i - 2, j, k));
                        riemann_para[2].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, 0, m_idx_proxy->GetIdx(i - 1, j, k));
                        riemann_para[2].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, 0, m_idx_proxy->GetIdx(i - 1, j, k));
                        riemann_para[3].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, 0, m_idx_proxy->GetIdx(i, j, k));
                        riemann_para[3].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, 0, m_idx_proxy->GetIdx(i, j, k));
                        riemann_para[4].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, 0, m_idx_proxy->GetIdx(i + 1, j, k));
                        riemann_para[4].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, 0, m_idx_proxy->GetIdx(i + 1, j, k));
                        riemann_para[5].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, 0, m_idx_proxy->GetIdx(i + 2, j, k));
                        riemann_para[5].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, 0, m_idx_proxy->GetIdx(i + 2, j, k));
                    }
                    for (int i = 0; i < 6; ++i)
                    {
                        m_riemann_solver->Solver(riemann_para[i]);
                    }
                    for (int iVar = 0; iVar < equ_num; ++iVar)
                    {
                        res_tmp[iVar] -= 75.0 / 64.0 * (riemann_para[3].flux[iVar] - riemann_para[2].flux[iVar]) - 25.0 / 384.0 * (riemann_para[4].flux[iVar] - riemann_para[1].flux[iVar]) + 3.0 / 640.0 * (riemann_para[5].flux[iVar] - riemann_para[0].flux[iVar]);
                    }

                    // j direction
                    for (int i = 0;i < 6;++i)
                    {
                        riemann_para[i].norm(0) = m_node_metrics->GetEta(idx)[0];
                        riemann_para[i].norm(1) = m_node_metrics->GetEta(idx)[1];
                        riemann_para[i].norm(2) = m_node_metrics->GetEta(idx)[2];
                        riemann_para[i].nt = m_node_metrics->GetEta(idx)[3];
                    }
                    for (int iVal = 0; iVal < equ_num; ++iVal)
                    {
                        riemann_para[0].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, 1, m_idx_proxy->GetIdx(i, j - 3, k));
                        riemann_para[0].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, 1, m_idx_proxy->GetIdx(i, j - 3, k));
                        riemann_para[1].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, 1, m_idx_proxy->GetIdx(i, j - 2, k));
                        riemann_para[1].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, 1, m_idx_proxy->GetIdx(i, j - 2, k));
                        riemann_para[2].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, 1, m_idx_proxy->GetIdx(i, j - 1, k));
                        riemann_para[2].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, 1, m_idx_proxy->GetIdx(i, j - 1, k));
                        riemann_para[3].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, 1, m_idx_proxy->GetIdx(i, j, k));
                        riemann_para[3].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, 1, m_idx_proxy->GetIdx(i, j, k));
                        riemann_para[4].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, 1, m_idx_proxy->GetIdx(i, j + 1, k));
                        riemann_para[4].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, 1, m_idx_proxy->GetIdx(i, j + 1, k));
                        riemann_para[5].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, 1, m_idx_proxy->GetIdx(i, j + 2, k));
                        riemann_para[5].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, 1, m_idx_proxy->GetIdx(i, j + 2, k));
                    }
                    for (int i = 0; i < 6; ++i)
                    {
                        m_riemann_solver->Solver(riemann_para[i]);
                    }
                    for (int iVar = 0; iVar < equ_num; ++iVar)
                    {
                        res_tmp[iVar] -= 75.0 / 64.0 * (riemann_para[3].flux[iVar] - riemann_para[2].flux[iVar]) - 25.0 / 384.0 * (riemann_para[4].flux[iVar] - riemann_para[1].flux[iVar]) + 3.0 / 640.0 * (riemann_para[5].flux[iVar] - riemann_para[0].flux[iVar]);
                    }
                    // k direction
                    if (grid->GetDim() == 3)
                    {
                        for (int i = 0;i < 6;++i)
                        {
                            riemann_para[i].norm(0) = m_node_metrics->GetZeta(idx)[0];
                            riemann_para[i].norm(1) = m_node_metrics->GetZeta(idx)[1];
                            riemann_para[i].norm(2) = m_node_metrics->GetZeta(idx)[2];
                            riemann_para[i].nt = m_node_metrics->GetZeta(idx)[3];
                        }
                        for (int iVal = 0; iVal < equ_num; ++iVal)
                        {
                            riemann_para[0].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, 2, m_idx_proxy->GetIdx(i, j, k - 3));
                            riemann_para[0].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, 2, m_idx_proxy->GetIdx(i, j, k - 3));
                            riemann_para[1].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, 2, m_idx_proxy->GetIdx(i, j, k - 2));
                            riemann_para[1].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, 2, m_idx_proxy->GetIdx(i, j, k - 2));
                            riemann_para[2].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, 2, m_idx_proxy->GetIdx(i, j, k - 1));
                            riemann_para[2].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, 2, m_idx_proxy->GetIdx(i, j, k - 1));
                            riemann_para[3].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, 2, m_idx_proxy->GetIdx(i, j, k));
                            riemann_para[3].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, 2, m_idx_proxy->GetIdx(i, j, k));
                            riemann_para[4].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, 2, m_idx_proxy->GetIdx(i, j, k + 1));
                            riemann_para[4].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, 2, m_idx_proxy->GetIdx(i, j, k + 1));
                            riemann_para[5].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, 2, m_idx_proxy->GetIdx(i, j, k + 2));
                            riemann_para[5].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, 2, m_idx_proxy->GetIdx(i, j, k + 2));
                        }
                        for (int i = 0; i < 6; ++i)
                        {
                            m_riemann_solver->Solver(riemann_para[i]);
                        }
                        for (int iVar = 0; iVar < equ_num; ++iVar)
                        {
                            res_tmp[iVar] -= 75.0 / 64.0 * (riemann_para[3].flux[iVar] - riemann_para[2].flux[iVar]) - 25.0 / 384.0 * (riemann_para[4].flux[iVar] - riemann_para[1].flux[iVar]) + 3.0 / 640.0 * (riemann_para[5].flux[iVar] - riemann_para[0].flux[iVar]);
                        }
                    }
                    for (int iVar = 0; iVar < equ_num; ++iVar)
                    {
                        data_manager->SetResidual(iVar, idx, res_tmp[iVar] / jacobi);
                    }
                }
            }
        }
    }
    void NSSolverStruct::ConvectiveResidual()
    {
        // ConvectiveResidualMUSCL();
        ConvectiveResidualWCNS5();
        return;
        auto grid = GetGrid();
        auto para = GetPara();
        auto data_manager = GetDataManager();
        auto gas = GetGas();
        auto equ_num = para->GetEquNum();
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        RiemannSolverPara riemann_para[6];
        for (int i = 0; i < 6; ++i)
        {
            riemann_para[i].gamma_left = riemann_para[i].gamma_right = gas->GetGamma();
        }
        int idx;
        for (int i = is;i <= ie;++i)
        {
            for (int j = js;j <= je;++j)
            {
                for (int k = ks;k <= ke;++k)
                {
                    idx = m_idx_proxy->GetIdx(i, j, k);
                    double jacobi = m_node_metrics->GetJacobian(idx);
                    // i direction
                    riemann_para[0].norm(0) = m_node_metrics->GetXi(idx)[0];
                    riemann_para[0].norm(1) = m_node_metrics->GetXi(idx)[1];
                    riemann_para[0].norm(2) = m_node_metrics->GetXi(idx)[2];
                    riemann_para[0].nt = m_node_metrics->GetXi(idx)[3];
                    riemann_para[1].norm(0) = m_node_metrics->GetXi(idx)[0];
                    riemann_para[1].norm(1) = m_node_metrics->GetXi(idx)[1];
                    riemann_para[1].norm(2) = m_node_metrics->GetXi(idx)[2];
                    riemann_para[1].nt = m_node_metrics->GetXi(idx)[3];
                    MidNode1st(idx, m_idx_proxy->GetIdx(i + 1, j, k), &riemann_para[0].prim_left(0), &riemann_para[0].prim_right(0));
                    MidNode1st(m_idx_proxy->GetIdx(i - 1, j, k), idx, &riemann_para[1].prim_left(0), &riemann_para[1].prim_right(0));
                    // j direction
                    riemann_para[2].norm(0) = m_node_metrics->GetEta(idx)[0];
                    riemann_para[2].norm(1) = m_node_metrics->GetEta(idx)[1];
                    riemann_para[2].norm(2) = m_node_metrics->GetEta(idx)[2];
                    riemann_para[2].nt = m_node_metrics->GetEta(idx)[3];
                    riemann_para[3].norm(0) = m_node_metrics->GetEta(idx)[0];
                    riemann_para[3].norm(1) = m_node_metrics->GetEta(idx)[1];
                    riemann_para[3].norm(2) = m_node_metrics->GetEta(idx)[2];
                    riemann_para[3].nt = m_node_metrics->GetEta(idx)[3];
                    MidNode1st(idx, m_idx_proxy->GetIdx(i, j + 1, k), &riemann_para[2].prim_left(0), &riemann_para[2].prim_right(0));
                    MidNode1st(m_idx_proxy->GetIdx(i, j - 1, k), idx, &riemann_para[3].prim_left(0), &riemann_para[3].prim_right(0));
                    // k direction
                    riemann_para[4].norm(0) = m_node_metrics->GetZeta(idx)[0];
                    riemann_para[4].norm(1) = m_node_metrics->GetZeta(idx)[1];
                    riemann_para[4].norm(2) = m_node_metrics->GetZeta(idx)[2];
                    riemann_para[4].nt = m_node_metrics->GetZeta(idx)[3];
                    riemann_para[5].norm(0) = m_node_metrics->GetZeta(idx)[0];
                    riemann_para[5].norm(1) = m_node_metrics->GetZeta(idx)[1];
                    riemann_para[5].norm(2) = m_node_metrics->GetZeta(idx)[2];
                    riemann_para[5].nt = m_node_metrics->GetZeta(idx)[3];
                    MidNode1st(idx, m_idx_proxy->GetIdx(i, j, k + 1), &riemann_para[4].prim_left(0), &riemann_para[4].prim_right(0));
                    MidNode1st(m_idx_proxy->GetIdx(i, j, k - 1), idx, &riemann_para[5].prim_left(0), &riemann_para[5].prim_right(0));
                    for (int i = 0; i < 6; ++i)
                    {
                        m_riemann_solver->Solver(riemann_para[i]);
                    }
                    for (int iVar = 0; iVar < equ_num; ++iVar)
                    {
                        double flux = (riemann_para[0].flux[iVar] - riemann_para[1].flux[iVar] + riemann_para[2].flux[iVar] - riemann_para[3].flux[iVar] + riemann_para[4].flux[iVar] - riemann_para[5].flux[iVar]) / jacobi;
                        data_manager->SetResidual(iVar, idx, data_manager->GetResidual(iVar, idx) - flux);
                    }
                }
            }
        }
    }
    void NSSolverStruct::ConvectiveResidual1st()
    {
        CalcMidNode1st();
    }
    void NSSolverStruct::ConvectiveResidualGrad()
    {

    }
    void NSSolverStruct::ConvectiveResidualMUSCL()
    {
        CalcMidNodePrimMUSCL();
        CalcMidGhostNodePrimMUSCL();
        FluxDifferenceMUSCL();
    }
    void NSSolverStruct::ConvectiveResidualWCNS5()
    {
        CalcMidNodePrimWCNS5();
        CalcMidGhostNodePrimWCNS5();
        FluxDifferenceWCNS5();
    }
    void NSSolverStruct::ViscousResidual()
    {
        //todo
    }
    void NSSolverStruct::CalcViscousFlux()
    {
        //todo
    }
    void NSSolverStruct::CalcViscousFluxGrad()
    {
        //todo
    }
    void NSSolverStruct::SourceResidual()
    {
        //todo
    }
    void NSSolverStruct::BoundaryCondition()
    {
        auto grid = GetGrid();
        auto data_manager = GetDataManager();
        auto bound_map = grid->GetBoundMap();
        for (auto& boundary : bound_map->GetBoundaryMap())
        {
            auto& bound_name = boundary.first;
            auto& bound = boundary.second;
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
    void NSSolverStruct::InletBC(BoundStruct& bound)
    {
        auto grid = GetGrid();
        auto data_manager = GetDataManager();
        int ghost_size = grid->GetGhostSize();
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
        for (int iGhost = 1;iGhost <= ghost_size;++iGhost)
        {
            i_ghost = i_bound + iGhost * bound_direction[0];
            j_ghost = j_bound + iGhost * bound_direction[1];
            k_ghost = k_bound + iGhost * bound_direction[2];
            idx_ghost = m_idx_proxy->GetIdx(i_ghost, j_ghost, k_ghost);
            for (int iVal = 0; iVal < 5; ++iVal)
            {
                data_manager->SetPrim(iVal, idx_ghost, prim_far[iVal]);
                data_manager->SetCons(iVal, idx_ghost, cons_far[iVal]);
            }
        }
    }
    void NSSolverStruct::OutletBC(BoundStruct& bound)
    {
        auto grid = GetGrid();
        auto data_manager = GetDataManager();
        int ghost_size = grid->GetGhostSize();
        int i_bound, j_bound, k_bound;
        bound.GetIdxBound(i_bound, j_bound, k_bound);
        int idx_bound = m_idx_proxy->GetIdx(i_bound, j_bound, k_bound);
        int i_ghost, j_ghost, k_ghost;
        for (int iGhost = 1;iGhost <= ghost_size;++iGhost)
        {
            i_ghost = i_bound + iGhost * bound.GetDirection()[0];
            j_ghost = j_bound + iGhost * bound.GetDirection()[1];
            k_ghost = k_bound + iGhost * bound.GetDirection()[2];
            int idx_ghost = m_idx_proxy->GetIdx(i_ghost, j_ghost, k_ghost);
            for (int iVal = 0; iVal < 5; ++iVal)
            {
                data_manager->SetPrim(iVal, idx_ghost, data_manager->GetPrim(iVal, idx_bound));
                data_manager->SetCons(iVal, idx_ghost, data_manager->GetCons(iVal, idx_bound));
            }
        }
    }
    void NSSolverStruct::WallBC(BoundStruct& bound)
    {
        auto grid = GetGrid();
        auto data_manager = GetDataManager();
        int ghost_size = grid->GetGhostSize();
        int i_bound, j_bound, k_bound;
        bound.GetIdxBound(i_bound, j_bound, k_bound);
        int idx_bound = m_idx_proxy->GetIdx(i_bound, j_bound, k_bound);
        int i_ref, j_ref, k_ref;
        int i_ghost, j_ghost, k_ghost;
        for (int iGhost = 1;iGhost <= ghost_size;++iGhost)
        {
            i_ghost = i_bound + iGhost * bound.GetDirection()[0];
            j_ghost = j_bound + iGhost * bound.GetDirection()[1];
            k_ghost = k_bound + iGhost * bound.GetDirection()[2];
            int idx_ghost = m_idx_proxy->GetIdx(i_ghost, j_ghost, k_ghost);
            i_ref = i_ghost + bound.GetDirection()[0];
            j_ref = j_ghost + bound.GetDirection()[1];
            k_ref = k_ghost + bound.GetDirection()[2];
            int idx_ref = m_idx_proxy->GetIdx(i_ref, j_ref, k_ref);
            for (int iVal = 0; iVal < 5; ++iVal)
            {
                data_manager->SetPrim(iVal, idx_bound, data_manager->GetPrim(iVal, idx_ghost));
                data_manager->SetCons(iVal, idx_bound, data_manager->GetCons(iVal, idx_ghost));
            }
        }
    }
    void NSSolverStruct::RiemannBC(BoundStruct& bound)
    {

    }
    void NSSolverStruct::SymmetryBC(BoundStruct& bound)
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
    void NSSolverStruct::BackupField(std::string& back_folder)
    {

    }
    GridStruct* NSSolverStruct::GetGrid()
    {
        return static_cast<GridStruct*>(FlowSolver::GetGrid());
    }
} // namespace  zaran