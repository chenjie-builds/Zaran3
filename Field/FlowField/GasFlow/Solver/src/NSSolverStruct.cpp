#include "NSSolverStruct.h"
namespace  zaran
{
    NSSolverStruct::NSSolverStruct(int index, string name, FlowSolverPara* para, GridStruct* grid, FieldData* fieldData, DataManagerNS* data_manager)
        :NSSolver(index, name, para, grid, fieldData), m_data_manager(data_manager)
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
    void NSSolverStruct::InitFieldFarFlow()
    {
        auto grid = GetGrid();
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
                            m_data_manager->SetPrim(iVal, idx, prim_far[iVal]);
                        }
                    }
                }
            }
        }
    }
    void NSSolverStruct::InitFieldFarFieldNoVelocity()
    {
        auto grid = GetGrid();
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
                        m_data_manager->SetPrim(iVal, idx, prim_far[iVal]);
                    }
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
    void NSSolverStruct::Preprocess()
    {
        NSSolver::Preprocess();
    }
    void NSSolverStruct::CalcTimeStepLocal()
    {
        auto grid = GetGrid();
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
                for (int i = is; i <=ie; ++i)
                {
                    int idx = m_idx_proxy->GetIdx(i, j, k);
                    auto xi = m_node_metrics->GetXi(idx);
                    auto eta = m_node_metrics->GetEta(idx);
                    auto zeta = m_node_metrics->GetZeta(idx);
                    auto jacobi = m_node_metrics->GetJacobian(idx);
                    double c = sqrt(gamma * m_data_manager->GetPressure(idx) / m_data_manager->GetDensity(idx));
                    double norm_xi = sqrt(xi[0] * xi[0] + xi[1] * xi[1] + xi[2] * xi[2]);
                    double norm_eta = sqrt(eta[0] * eta[0] + eta[1] * eta[1] + eta[2] * eta[2]);
                    double norm_zeta = sqrt(zeta[0] * zeta[0] + zeta[1] * zeta[1] + zeta[2] * zeta[2]);
                    double u_xi = m_data_manager->GetVelocity(0, idx) * xi[0] + m_data_manager->GetVelocity(1, idx) * xi[1] + m_data_manager->GetVelocity(2, idx) * xi[2];
                    double u_eta = m_data_manager->GetVelocity(0, idx) * eta[0] + m_data_manager->GetVelocity(1, idx) * eta[1] + m_data_manager->GetVelocity(2, idx) * eta[2];
                    double u_zeta = m_data_manager->GetVelocity(0, idx) * zeta[0] + m_data_manager->GetVelocity(1, idx) * zeta[1] + m_data_manager->GetVelocity(2, idx) * zeta[2];
                    double lamda = abs(u_xi) + abs(u_eta) + abs(u_zeta) + c * (norm_xi + norm_eta + norm_zeta);
                    m_data_manager->SetTimeStep(idx, cfl / lamda);
                    if (m_data_manager->GetTimeStep(idx) < min_dt)
                    {
                        min_dt = m_data_manager->GetTimeStep(idx);
                    }
                }
            }
        }
        GlobalData::Update("dt", min_dt);
    }
    void NSSolverStruct::ReduceTimeStep(double& dt)
    {
        auto grid = GetGrid();
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
                    m_data_manager->SetTimeStep(idx, dt);
                }
            }
        }
    }
    void NSSolverStruct::RungeKutta()
    {
        auto grid = GetGrid();
        auto para = GetPara();
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
                        dt = m_data_manager->GetTimeStep(idx);
                        jacobi = m_node_metrics->GetJacobian(idx);
                        for (int iVal = 0; iVal < 5; ++iVal)
                        {
                            m_data_manager->SetCons(iVal, idx, m_data_manager->GetCons(iVal, idx) + dt * rk_coef[iStage] * m_data_manager->GetResidual(iVal, idx) * jacobi);
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
                        prim[iVal] = m_data_manager->GetPrim(iVal, idx);
                    }
                    gas->Prim2Cons(prim, cons);
                    for (int iEqu = 0; iEqu < GetPara()->GetEquNum(); ++iEqu)
                    {
                        m_data_manager->SetCons(iEqu, idx, cons[iEqu]);
                    }
                }
            }
        }

    }
    void NSSolverStruct::Cons2Prim()
    {
        auto gas = GetGas();
        auto grid = GetGrid();
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
                        cons[iVal] = m_data_manager->GetCons(iVal, idx);
                    }
                    gas->Cons2Prim(cons, prim);
                    for (int iEqu = 0; iEqu < GetPara()->GetEquNum(); ++iEqu)
                    {
                        m_data_manager->SetPrim(iEqu, idx, prim[iEqu]);
                    }
                }
            }
        }
    }
    void NSSolverStruct::ZeroResidual()
    {
        auto grid = GetGrid();
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
                        m_data_manager->SetResidual(iVal, idx, 0.0);
                    }
                }
            }
        }
    }
    void NSSolverStruct::MidNodeRec2nd(int idx_left, int idx_right, double* value_rec_left, double* value_rec_right)
    {
        int equ_num = GetPara()->GetEquNum();
        auto grid = GetGrid();
        auto node = grid->GetNode();
        int il, jl, kl, ir, jr, kr;
        m_idx_proxy->GetIdxStruct(idx_left, il, jl, kl);
        m_idx_proxy->GetIdxStruct(idx_right, ir, jr, kr);
        double vec_node2neighbor[3];
        vec_node2neighbor[0] = node->GetCoord(il, jl, kl)[0] - node->GetCoord(ir, jr, kr)[0]; // (x_left - x_right)
        vec_node2neighbor[1] = node->GetCoord(il, jl, kl)[1] - node->GetCoord(ir, jr, kr)[1]; // (y_left - y_right)
        vec_node2neighbor[2] = node->GetCoord(il, jl, kl)[2] - node->GetCoord(ir, jr, kr)[2]; // (z_left - z_right)
        for (int iVal = 0; iVal < equ_num; ++iVal)
        {
            value_rec_left[iVal] = m_data_manager->GetPrim(iVal, idx_left) + 0.5 * m_data_manager->GetLimiter(iVal, idx_left) *
                (vec_node2neighbor[0] * m_data_manager->GetPrimGrad(iVal, 0, idx_left) +
                    vec_node2neighbor[1] * m_data_manager->GetPrimGrad(iVal, 1, idx_left) +
                    vec_node2neighbor[2] * m_data_manager->GetPrimGrad(iVal, 2, idx_left));
            value_rec_right[iVal] = m_data_manager->GetPrim(iVal, idx_right) - 0.5 * m_data_manager->GetLimiter(iVal, idx_right) *
                (vec_node2neighbor[0] * m_data_manager->GetPrimGrad(iVal, 0, idx_right) +
                    vec_node2neighbor[1] * m_data_manager->GetPrimGrad(iVal, 1, idx_right) +
                    vec_node2neighbor[2] * m_data_manager->GetPrimGrad(iVal, 2, idx_right));
        }
    }
    void NSSolverStruct::MidNodeRec1st(int index_left, int index_right, double* value_rec_left, double* value_rec_right)
    {
        int equ_num = GetPara()->GetEquNum();
        auto  grid = GetGrid();
        for (int iVal = 0; iVal < equ_num; ++iVal)
        {
            value_rec_left[iVal] = m_data_manager->GetPrim(iVal, index_left);
            value_rec_right[iVal] = m_data_manager->GetPrim(iVal, index_right);
        }
    }
    void NSSolverStruct::ConvectiveResidual()
    {
        auto grid = GetGrid();
        auto para = GetPara();
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
                    MidNodeRec2nd(idx, m_idx_proxy->GetIdx(i + 1, j, k), &riemann_para[0].prim_left(0), &riemann_para[0].prim_right(0));
                    MidNodeRec2nd(m_idx_proxy->GetIdx(i - 1, j, k), idx, &riemann_para[1].prim_left(0), &riemann_para[1].prim_right(0));
                    // j direction
                    riemann_para[2].norm(0) = m_node_metrics->GetEta(idx)[0];
                    riemann_para[2].norm(1) = m_node_metrics->GetEta(idx)[1];
                    riemann_para[2].norm(2) = m_node_metrics->GetEta(idx)[2];
                    riemann_para[2].nt = m_node_metrics->GetEta(idx)[3];
                    riemann_para[3].norm(0) = m_node_metrics->GetEta(idx)[0];
                    riemann_para[3].norm(1) = m_node_metrics->GetEta(idx)[1];
                    riemann_para[3].norm(2) = m_node_metrics->GetEta(idx)[2];
                    riemann_para[3].nt = m_node_metrics->GetEta(idx)[3];
                    MidNodeRec2nd(idx, m_idx_proxy->GetIdx(i, j + 1, k), &riemann_para[2].prim_left(0), &riemann_para[2].prim_right(0));
                    MidNodeRec2nd(m_idx_proxy->GetIdx(i, j - 1, k), idx, &riemann_para[3].prim_left(0), &riemann_para[3].prim_right(0));
                    // k direction
                    riemann_para[4].norm(0) = m_node_metrics->GetZeta(idx)[0];
                    riemann_para[4].norm(1) = m_node_metrics->GetZeta(idx)[1];
                    riemann_para[4].norm(2) = m_node_metrics->GetZeta(idx)[2];
                    riemann_para[4].nt = m_node_metrics->GetZeta(idx)[3];
                    riemann_para[5].norm(0) = m_node_metrics->GetZeta(idx)[0];
                    riemann_para[5].norm(1) = m_node_metrics->GetZeta(idx)[1];
                    riemann_para[5].norm(2) = m_node_metrics->GetZeta(idx)[2];
                    riemann_para[5].nt = m_node_metrics->GetZeta(idx)[3];
                    MidNodeRec2nd(idx, m_idx_proxy->GetIdx(i, j, k + 1), &riemann_para[4].prim_left(0), &riemann_para[4].prim_right(0));
                    MidNodeRec2nd(m_idx_proxy->GetIdx(i, j, k - 1), idx, &riemann_para[5].prim_left(0), &riemann_para[5].prim_right(0));
                    for (int i = 0; i < 6; ++i)
                    {
                        m_riemann_solver->Solver(riemann_para[i]);
                    }
                    for (int iVar = 0; iVar < equ_num; ++iVar)
                    {
                        double flux = (riemann_para[0].flux[iVar] - riemann_para[1].flux[iVar] + riemann_para[2].flux[iVar] - riemann_para[3].flux[iVar] + riemann_para[4].flux[iVar] - riemann_para[5].flux[iVar]) / jacobi;
                        m_data_manager->SetResidual(iVar, idx, m_data_manager->GetResidual(iVar, idx) - flux);
                    }
                }
            }
        }
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
        int i_ghost, j_ghost, k_ghost;
        bound.GetIdxGhost(i_ghost, j_ghost, k_ghost);
        int idx_ghost = m_idx_proxy->GetIdx(i_ghost, j_ghost, k_ghost);
        double prim_far[5];
        prim_far[0] = GetPara()->GetInflowDensity();
        prim_far[1] = GetPara()->GetInflowVelocityX();
        prim_far[2] = GetPara()->GetInflowVelocityY();
        prim_far[3] = GetPara()->GetInflowVelocityZ();
        prim_far[4] = GetPara()->GetInflowPressure();
        double cons_far[5];
        GetGas()->Prim2Cons(prim_far, cons_far);
        for (int iVal = 0; iVal < 5; ++iVal)
        {
            m_data_manager->SetPrim(iVal, idx_ghost, prim_far[iVal]);
            m_data_manager->SetCons(iVal, idx_ghost, cons_far[iVal]);
        }
    }
    void NSSolverStruct::OutletBC(BoundStruct& bound)
    {
        int i_bound, j_bound, k_bound;
        bound.GetIdxBound(i_bound, j_bound, k_bound);
        int idx_bound = m_idx_proxy->GetIdx(i_bound, j_bound, k_bound);
        int i_ghost, j_ghost, k_ghost;
        bound.GetIdxGhost(i_ghost, j_ghost, k_ghost);
        int idx_ghost = m_idx_proxy->GetIdx(i_ghost, j_ghost, k_ghost);
        for (int iVal = 0; iVal < 5; ++iVal)
        {
            m_data_manager->SetPrim(iVal, idx_ghost, m_data_manager->GetPrim(iVal, idx_bound));
            m_data_manager->SetCons(iVal, idx_ghost, m_data_manager->GetCons(iVal, idx_bound));
        }
    }
    void NSSolverStruct::WallBC(BoundStruct& bound)
    {
        int i_bound, j_bound, k_bound;
        bound.GetIdxBound(i_bound, j_bound, k_bound);
        int idx_bound = m_idx_proxy->GetIdx(i_bound, j_bound, k_bound);
        int i_ref, j_ref, k_ref;
        bound.GetIdxRef(i_ref, j_ref, k_ref);
        int idx_ref = m_idx_proxy->GetIdx(i_ref, j_ref, k_ref);
        int i_ghost, j_ghost, k_ghost;
        bound.GetIdxGhost(i_ghost, j_ghost, k_ghost);
        int idx_ghost = m_idx_proxy->GetIdx(i_ghost, j_ghost, k_ghost);
        for (int iVal = 0; iVal < 5; ++iVal)
        {
            m_data_manager->SetPrim(iVal, idx_bound, m_data_manager->GetPrim(iVal, idx_ghost));
            m_data_manager->SetCons(iVal, idx_bound, m_data_manager->GetCons(iVal, idx_ghost));
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