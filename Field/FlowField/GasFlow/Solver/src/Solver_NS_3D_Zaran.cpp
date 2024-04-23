#include "Solver_NS_3D_Zaran.h"
#include "Vanleer.h"
namespace zaran
{

    void Solver_NS_3D_Zaran::CreateData()
    {
        auto grid = GetGrid();
        int ni, nj, nk;
        grid->GetNodeNum(ni, nj, nk);
        int cell_num = (ni - 1) * (nj - 1) * (nk - 1);
        auto& dataPtr = GetFieldData();
        auto& data = *dataPtr;
        FieldDataType type = FieldDataType::real;
        //添加单元数据
        data.AddData("density", type, cell_num);
        data.AddData("velocity_x", type, cell_num);
        data.AddData("velocity_y", type, cell_num);
        data.AddData("velocity_w", type, cell_num);
        data.AddData("pressure", type, cell_num);
        data.AddData("conservative", type, cell_num * GetNumberOfEquations());
        data.AddData("residual", type, cell_num * GetNumberOfEquations());
        data.AddData("dt", type, cell_num);
        data.AddData("limiter", type, cell_num * GetNumberOfEquations());
        data.AddData("prim_grad", type, cell_num * GetNumberOfEquations() * 3);
        data.AddData("consRK0", type, cell_num);
        data.AddData("consRK1", type, cell_num);
        data.AddData("consRK2", type, cell_num);
        data.AddData("consRK3", type, cell_num);
        data.AddData("consRK4", type, cell_num);
        data.AddData("metrics", type, cell_num * 17);
        //添加边界节点数据
        auto& bound_patch = grid->GetBoundPatch();
        int bound_node_num = bound_patch.GetPatchNum();
        data.AddData("density_patch", type, bound_node_num);
        data.AddData("velocity_x_patch", type, bound_node_num);
        data.AddData("velocity_y_patch", type, bound_node_num);
        data.AddData("velocity_z_patch", type, bound_node_num);
        data.AddData("pressure_patch", type, bound_node_num);
        data.AddData("prim_grad_patch", type, bound_node_num * GetNumberOfEquations() * 3);
        data.AddData("limiter_patch", type, bound_node_num * GetNumberOfEquations());

    }

    void Solver_NS_3D_Zaran::RegisterFieldData()
    {
        NSSolver::RegisterFieldData();
        auto& data = *GetFieldData();
        m_prim_bound = new double* [5];
        auto grid = GetGrid();
        auto& bound_patch = grid->GetBoundPatch();
        int bound_node_num = bound_patch.GetPatchNum();
        data.GetData("density_patch", m_prim_bound[0]);
        data.GetData("velocity_x_patch", m_prim_bound[1]);
        data.GetData("velocity_y_patch", m_prim_bound[2]);
        data.GetData("velocity_y_patch", m_prim_bound[3]);
        data.GetData("pressure_patch", m_prim_bound[4]);
        data.GetData("prim_grad_patch", m_prim_grad_bound);
        data.GetData("limiter_patch", m_limiter_bound);
    }
    void Solver_NS_3D_Zaran::InitData()
    {
        NSSolver::InitData();
        FlowSolverPara* para = GetPara();
        double prim_init[5];
        prim_init[0] = para->GetInflowDensity();
        prim_init[1] = para->GetInflowVelocityX();
        prim_init[2] = para->GetInflowVelocityY();
        prim_init[3] = para->GetInflowVelocityZ();
        prim_init[4] = para->GetInflowPressure();
        auto grid = GetGrid();
        int n_patch = grid->GetBoundPatch().GetPatchNum();
        for (int iPatch = 0;iPatch < n_patch;++iPatch)
        {
            m_prim_bound[0][iPatch] = prim_init[0];
            m_prim_bound[1][iPatch] = prim_init[1];
            m_prim_bound[2][iPatch] = prim_init[2];
            m_prim_bound[3][iPatch] = prim_init[3];
            m_prim_bound[4][iPatch] = prim_init[4];
        }
    }
    void Solver_NS_3D_Zaran::CalcMetric()
    {
        CalcMetricStruct();
        CalcMetricMid();
    }
    void Solver_NS_3D_Zaran::CalcGradWLS()
    {
        CalcGradStructWLS();
        CalcGradMidWLS();
    }

    void Solver_NS_3D_Zaran::CalcTimeStepLocal()
    {
        CalcTimeStepLocalStruct();
        CalcTimeStepLocalMid();
    }

    void Solver_NS_3D_Zaran::ConvectiveResidual()
    {
        InviscidFluxStruct();
        InviscidFluxMid();
    }

    void Solver_NS_3D_Zaran::ViscousResidual()
    {
        ViscousFluxStruct();
        ViscousFluxMid();
    }

    void Solver_NS_3D_Zaran::SourceTermResidual()
    {
        //TODO
    }

    void Solver_NS_3D_Zaran::CalcLimiterVK()
    {
        //TODO
    }

    void Solver_NS_3D_Zaran::CalcLimiterBJ()
    {
        //TODO
    }


    void Solver_NS_3D_Zaran::CalcLimiterNone()
    {
        //TODO
    }

    void Solver_NS_3D_Zaran::CalcLimiterFirstOrder()
    {
        //TODO
    }

    void Solver_NS_3D_Zaran::CalcLimiterBound()
    {
        //TODO
    }
    void Solver_NS_3D_Zaran::RungeKutta()
    {
        auto& grid = GetGrid();
        FlowSolverPara* para = GetPara();
        const DArray& rkCoef = para->GetRKCoef();
        int rkStage = rkCoef.size();
        int ni, nj, nk;
        grid->GetNodeNum(ni, nj, nk);
        int cell_num = (ni - 1) * (nj - 1) * (nk - 1);
        for (int iStage = 0; iStage < rkStage; ++iStage)
        {
            CalcResidual();
#pragma omp parallel for
            for (int iNode = 0; iNode < cell_num; ++iNode)
            {
                for (int iVal = 0; iVal < 5; ++iVal)
                {
                    GetCons(iNode, iVal) -= rkCoef[iStage] * m_dt[iNode] * GetResidual(iNode, iVal);
                }
            }
        }
    }
    void Solver_NS_3D_Zaran::CalcMetricStruct()
    {
        auto grid = GetGrid();
        // 起始点和终止点的编号,s: start, e: end
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        int iCell;
        // grid->GetNodeIndex(i, j, k)的lamda表达式
        auto CellIndex = [&](int i, int j, int k) {return grid->GetCellIndex(i, j, k); };
        auto& cell_topo = grid->GetCellTopo();
        auto& cell_center_coord = cell_topo->GetCenterCoord();
        CoordTrans coordTrans;
        for (int k = ks; k < ke - 1; k++)
        {
            for (int j = js; j < je - 1; j++)
            {
                for (int i = is; i < ie - 1; i++)
                {
                    iCell = CellIndex(i, j, k);
                    coordTrans.CalcCoordTrans(3, cell_center_coord[CellIndex(i + 1, j, k)], cell_center_coord[CellIndex(i - 1, j, k)], cell_center_coord[CellIndex(i, j + 1, k)], cell_center_coord[CellIndex(i, j - 1, k)], cell_center_coord[CellIndex(i, j, k + 1)], cell_center_coord[CellIndex(i, j, k - 1)]);
                    GetMetricXi(iCell)[0] = coordTrans.GetXi()[0];
                    GetMetricXi(iCell)[1] = coordTrans.GetXi()[1];
                    GetMetricXi(iCell)[2] = coordTrans.GetXi()[2];
                    GetMetricXi(iCell)[3] = coordTrans.GetXi()[3];
                    GetMetricEta(iCell)[0] = coordTrans.GetEta()[0];
                    GetMetricEta(iCell)[1] = coordTrans.GetEta()[1];
                    GetMetricEta(iCell)[2] = coordTrans.GetEta()[2];
                    GetMetricEta(iCell)[3] = coordTrans.GetEta()[3];
                    GetMetricZeta(iCell)[0] = coordTrans.GetZeta()[0];
                    GetMetricZeta(iCell)[1] = coordTrans.GetZeta()[1];
                    GetMetricZeta(iCell)[2] = coordTrans.GetZeta()[2];
                    GetMetricZeta(iCell)[3] = coordTrans.GetZeta()[3];
                    GetMetricTau(iCell)[0] = coordTrans.GetTau()[0];
                    GetMetricTau(iCell)[1] = coordTrans.GetTau()[1];
                    GetMetricTau(iCell)[2] = coordTrans.GetTau()[2];
                    GetMetricTau(iCell)[3] = coordTrans.GetTau()[3];
                    GetMetricJacob(iCell) = coordTrans.J();
                }
            }
        }
    }

    void Solver_NS_3D_Zaran::CalcMetricMid()
    {
        auto grid = GetGrid();
        auto& bound_patch = grid->GetBoundPatch();
        auto& bound_node_index = bound_patch.GetIndex();
        auto& bound_node_coord = bound_patch.GetCoordinate();
        // 起始点和终止点的编号,s: start, e: end
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        int iCell;
        // grid->GetNodeIndex(i, j, k)的lamda表达式
        auto CellIndex = [&](int i, int j, int k) {return grid->GetCellIndex(i, j, k); };
        auto& cell_topo = grid->GetCellTopo();
        auto& cell_center_coord = cell_topo->GetCenterCoord();
        auto& cell_type = cell_topo->GetType();
        CoordTrans coordTrans;
        int i, j, k;
        Array<DVector3D*> neighbor_node_coord(6);
        IArray neighbor_node_index(6);
        for (int iPatch = 0;iPatch < bound_patch.GetPatchNum();iPatch++)
        {
            i = bound_node_index[iPatch][0];
            j = bound_node_index[iPatch][1];
            k = bound_node_index[iPatch][2];
            iCell = CellIndex(i, j, k);
            neighbor_node_index[0] = CellIndex(i + 1, j, k);
            neighbor_node_index[1] = CellIndex(i - 1, j, k);
            neighbor_node_index[2] = CellIndex(i, j + 1, k);
            neighbor_node_index[3] = CellIndex(i, j - 1, k);
            neighbor_node_index[4] = CellIndex(i, j, k + 1);
            neighbor_node_index[5] = CellIndex(i, j, k - 1);
            for (int iNeighbor = 0;iNeighbor < 6;iNeighbor++)
            {
                if (cell_type[neighbor_node_index[iNeighbor]] != CellType::Solid)
                    neighbor_node_coord[iNeighbor] = &cell_center_coord[neighbor_node_index[iNeighbor]];
                else
                    neighbor_node_coord[iNeighbor] = &bound_node_coord[iPatch];
            }
            coordTrans.CalcCoordTrans(3, *neighbor_node_coord[0], *neighbor_node_coord[1], *neighbor_node_coord[2], *neighbor_node_coord[3], *neighbor_node_coord[4], *neighbor_node_coord[5]);
            GetMetricXi(iCell)[0] = coordTrans.GetXi()[0];
            GetMetricXi(iCell)[1] = coordTrans.GetXi()[1];
            GetMetricXi(iCell)[2] = coordTrans.GetXi()[2];
            GetMetricXi(iCell)[3] = coordTrans.GetXi()[3];
            GetMetricEta(iCell)[0] = coordTrans.GetEta()[0];
            GetMetricEta(iCell)[1] = coordTrans.GetEta()[1];
            GetMetricEta(iCell)[2] = coordTrans.GetEta()[2];
            GetMetricEta(iCell)[3] = coordTrans.GetEta()[3];
            GetMetricZeta(iCell)[0] = coordTrans.GetZeta()[0];
            GetMetricZeta(iCell)[1] = coordTrans.GetZeta()[1];
            GetMetricZeta(iCell)[2] = coordTrans.GetZeta()[2];
            GetMetricZeta(iCell)[3] = coordTrans.GetZeta()[3];
            GetMetricTau(iCell)[0] = coordTrans.GetTau()[0];
            GetMetricTau(iCell)[1] = coordTrans.GetTau()[1];
            GetMetricTau(iCell)[2] = coordTrans.GetTau()[2];
            GetMetricTau(iCell)[3] = coordTrans.GetTau()[3];
            GetMetricJacob(iCell) = coordTrans.J();
            if (isinf(coordTrans.J()))
                Log::error("i:{},j:{},k:{},J:{}", i, j, k, coordTrans.J());
        }
    }

    void Solver_NS_3D_Zaran::CalcGradStructWLS()
    {
        //TODO
    }

    void Solver_NS_3D_Zaran::CalcGradMidWLS()
    {
        //TODO
    }

    void Solver_NS_3D_Zaran::CalcTimeStepLocalStruct()
    {
        auto grid = GetGrid();
        // 起始点和终止点的编号,s: start, e: end
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        int iCell;
        // grid->GetNodeIndex(i, j, k)的lamda表达式
        auto CellIndex = [&](int i, int j, int k) {return grid->GetCellIndex(i, j, k); };
        auto& cell_topo = grid->GetCellTopo();
        auto& cell_type = cell_topo->GetType();
        FlowSolverPara* para = GetPara();
        double cfl = para->GetCflNumber();
        double min_dt = LARGE_NUMBER;
        double gamma = 1.4;
#pragma omp parallel for private(iCell)
        for (int k = ks; k < ke - 1; k++)
        {
            for (int j = js; j < je - 1; j++)
            {
                for (int i = is; i < ie - 1; i++)
                {
                    iCell = CellIndex(i, j, k);
                    if (cell_type[iCell] != CellType::Fluid && cell_type[iCell] != CellType::FluidSolid)
                        continue;
                    double c = sqrt(gamma * m_prim[4][iCell] / m_prim[0][iCell]);
                    double normXi = sqrt(GetMetricXi(iCell)[0] * GetMetricXi(iCell)[0] + GetMetricXi(iCell)[1] * GetMetricXi(iCell)[1] + GetMetricXi(iCell)[2] * GetMetricXi(iCell)[2]);
                    double normEta = sqrt(GetMetricEta(iCell)[0] * GetMetricEta(iCell)[0] + GetMetricEta(iCell)[1] * GetMetricEta(iCell)[1] + GetMetricEta(iCell)[2] * GetMetricEta(iCell)[2]);
                    double normZeta = sqrt(GetMetricZeta(iCell)[0] * GetMetricZeta(iCell)[0] + GetMetricZeta(iCell)[1] * GetMetricZeta(iCell)[1] + GetMetricZeta(iCell)[2] * GetMetricZeta(iCell)[2]);
                    double uXi = m_prim[1][iCell] * GetMetricXi(iCell)[0] + m_prim[2][iCell] * GetMetricXi(iCell)[1] + m_prim[3][iCell] * GetMetricXi(iCell)[2];
                    double uEta = m_prim[1][iCell] * GetMetricEta(iCell)[0] + m_prim[2][iCell] * GetMetricEta(iCell)[1] + m_prim[3][iCell] * GetMetricEta(iCell)[2];
                    double uZeta = m_prim[1][iCell] * GetMetricZeta(iCell)[0] + m_prim[2][iCell] * GetMetricZeta(iCell)[1] + m_prim[3][iCell] * GetMetricZeta(iCell)[2];
                    double lambda = abs(uXi) + abs(uEta) + abs(uZeta) + c * (normXi + normEta + normZeta);
                    m_dt[iCell] = cfl / lambda;
                    min_dt = Min(min_dt, m_dt[iCell]);
                }
            }
        }
        GlobalData::Update("dt", min_dt);
    }

    void Solver_NS_3D_Zaran::CalcTimeStepLocalMid()
    {
        //TODO
    }
    void Solver_NS_3D_Zaran::InviscidFluxStruct()
    {
        auto& grid = GetGrid();
        auto& cell_topo = grid->GetCellTopo();
        auto& cell_type = cell_topo->GetType();
        auto& cell_center_coord = cell_topo->GetCenterCoord();
        // 起始点和终止点的编号,s: start, e: end
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        int iCell;
        DVector3D r, grad;
        int left_index, right_index;
        RiemannSolverPara riemann_para;
        riemann_para.gamma_left = riemann_para.gamma_right = 1.4;
#pragma omp parallel for private( iCell, r, grad, riemann_para, left_index, right_index)
        for (int k = ks; k < ke - 1; k++)
        {
            for (int j = js; j < je - 1; j++)
            {
                for (int i = is; i < ie - 1; i++)
                {
                    iCell = grid->GetCellIndex(i, j, k);
                    if (cell_type[iCell] != CellType::Fluid)
                        continue;
                    auto& jacobi = GetMetricJacob(iCell);
                    // i direction
                    riemann_para.norm(0) = GetMetricXi(iCell)[0];
                    riemann_para.norm(1) = GetMetricXi(iCell)[1];
                    riemann_para.norm(2) = GetMetricXi(iCell)[2];
                    riemann_para.nt = GetMetricXi(iCell)[3];
                    left_index = grid->GetCellIndex(i, j, k);
                    right_index = grid->GetCellIndex(i + 1, j, k);
                    r[0] = cell_center_coord[right_index].x() - cell_center_coord[left_index].x();
                    r[1] = cell_center_coord[right_index].y() - cell_center_coord[left_index].y();
                    r[2] = cell_center_coord[right_index].z() - cell_center_coord[left_index].z();
                    for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                    {
                        grad(0) = GetPrimGrad(left_index, iVal, 0);
                        grad(1) = GetPrimGrad(left_index, iVal, 1);
                        grad(2) = GetPrimGrad(left_index, iVal, 2);
                        riemann_para.prim_left(iVal) = m_prim[iVal][left_index] + 0.5 * GetLimiter(left_index, iVal) * grad.dot(r);
                        grad(0) = GetPrimGrad(right_index, iVal, 0);
                        grad(1) = GetPrimGrad(right_index, iVal, 1);
                        grad(2) = GetPrimGrad(right_index, iVal, 2);
                        riemann_para.prim_right(iVal) = m_prim[iVal][right_index] - 0.5 * GetLimiter(right_index, iVal) * grad.dot(r);
                    }
                    riemannSolver_->Solver(riemann_para);

                    for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                    {
                        GetResidual(iCell, iVal) += riemann_para.flux(iVal) / jacobi;
                    }

                    left_index = grid->GetCellIndex(i - 1, j, k);
                    right_index = grid->GetCellIndex(i, j, k);

                    r[0] = cell_center_coord[right_index].x() - cell_center_coord[left_index].x();
                    r[1] = cell_center_coord[right_index].y() - cell_center_coord[left_index].y();
                    r[2] = cell_center_coord[right_index].z() - cell_center_coord[left_index].z();
                    for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                    {

                        grad(0) = GetPrimGrad(left_index, iVal, 0);
                        grad(1) = GetPrimGrad(left_index, iVal, 1);
                        grad(2) = GetPrimGrad(left_index, iVal, 2);
                        riemann_para.prim_left(iVal) = m_prim[iVal][left_index] + 0.5 * GetLimiter(left_index, iVal) * grad.dot(r);
                        grad(0) = GetPrimGrad(right_index, iVal, 0);
                        grad(1) = GetPrimGrad(right_index, iVal, 1);
                        grad(2) = GetPrimGrad(right_index, iVal, 2);
                        riemann_para.prim_right(iVal) = m_prim[iVal][right_index] - 0.5 * GetLimiter(right_index, iVal) * grad.dot(r);
                    }
                    riemannSolver_->Solver(riemann_para);
                    for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                    {
                        GetResidual(iCell, iVal) -= riemann_para.flux(iVal) / jacobi;
                    }
                    // j direction
                    riemann_para.norm(0) = GetMetricEta(iCell)[0];
                    riemann_para.norm(1) = GetMetricEta(iCell)[1];
                    riemann_para.norm(2) = GetMetricEta(iCell)[2];
                    riemann_para.nt = GetMetricEta(iCell)[3];
                    left_index = grid->GetCellIndex(i, j, k);
                    right_index = grid->GetCellIndex(i, j + 1, k);
                    r[0] = cell_center_coord[right_index].x() - cell_center_coord[left_index].x();
                    r[1] = cell_center_coord[right_index].y() - cell_center_coord[left_index].y();
                    r[2] = cell_center_coord[right_index].z() - cell_center_coord[left_index].z();
                    for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                    {
                        grad(0) = GetPrimGrad(left_index, iVal, 0);
                        grad(1) = GetPrimGrad(left_index, iVal, 1);
                        grad(2) = GetPrimGrad(left_index, iVal, 2);
                        riemann_para.prim_left(iVal) = m_prim[iVal][left_index] + 0.5 * GetLimiter(left_index, iVal) * grad.dot(r);
                        grad(0) = GetPrimGrad(right_index, iVal, 0);
                        grad(1) = GetPrimGrad(right_index, iVal, 1);
                        grad(2) = GetPrimGrad(right_index, iVal, 2);
                        riemann_para.prim_right(iVal) = m_prim[iVal][right_index] - 0.5 * GetLimiter(right_index, iVal) * grad.dot(r);
                    }
                    riemannSolver_->Solver(riemann_para);
                    for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                    {
                        GetResidual(iCell, iVal) += riemann_para.flux(iVal) / jacobi;
                    }
                    left_index = grid->GetCellIndex(i, j - 1, k);
                    right_index = grid->GetCellIndex(i, j, k);

                    r[0] = cell_center_coord[right_index].x() - cell_center_coord[left_index].x();
                    r[1] = cell_center_coord[right_index].y() - cell_center_coord[left_index].y();
                    r[2] = cell_center_coord[right_index].z() - cell_center_coord[left_index].z();
                    for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                    {
                        grad(0) = GetPrimGrad(left_index, iVal, 0);
                        grad(1) = GetPrimGrad(left_index, iVal, 1);
                        grad(2) = GetPrimGrad(left_index, iVal, 2);
                        riemann_para.prim_left(iVal) = m_prim[iVal][left_index] + 0.5 * GetLimiter(left_index, iVal) * grad.dot(r);
                        grad(0) = GetPrimGrad(right_index, iVal, 0);
                        grad(1) = GetPrimGrad(right_index, iVal, 1);
                        grad(2) = GetPrimGrad(right_index, iVal, 2);
                        riemann_para.prim_right(iVal) = m_prim[iVal][right_index] - 0.5 * GetLimiter(right_index, iVal) * grad.dot(r);
                    }
                    riemannSolver_->Solver(riemann_para);
                    for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                    {
                        GetResidual(iCell, iVal) -= riemann_para.flux(iVal) / jacobi;
                    }
                    // k direction
                    riemann_para.norm(0) = GetMetricZeta(iCell)[0];
                    riemann_para.norm(1) = GetMetricZeta(iCell)[1];
                    riemann_para.norm(2) = GetMetricZeta(iCell)[2];
                    riemann_para.nt = GetMetricZeta(iCell)[3];
                    left_index = grid->GetCellIndex(i, j, k);
                    right_index = grid->GetCellIndex(i, j, k + 1);
                    r[0] = cell_center_coord[right_index].x() - cell_center_coord[left_index].x();
                    r[1] = cell_center_coord[right_index].y() - cell_center_coord[left_index].y();
                    r[2] = cell_center_coord[right_index].z() - cell_center_coord[left_index].z();
                    for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                    {
                        grad(0) = GetPrimGrad(left_index, iVal, 0);
                        grad(1) = GetPrimGrad(left_index, iVal, 1);
                        grad(2) = GetPrimGrad(left_index, iVal, 2);
                        riemann_para.prim_left(iVal) = m_prim[iVal][left_index] + 0.5 * GetLimiter(left_index, iVal) * grad.dot(r);
                        grad(0) = GetPrimGrad(right_index, iVal, 0);
                        grad(1) = GetPrimGrad(right_index, iVal, 1);
                        grad(2) = GetPrimGrad(right_index, iVal, 2);
                        riemann_para.prim_right(iVal) = m_prim[iVal][right_index] - 0.5 * GetLimiter(right_index, iVal) * grad.dot(r);
                    }
                    riemannSolver_->Solver(riemann_para);
                    for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                    {
                        GetResidual(iCell, iVal) += riemann_para.flux(iVal) / jacobi;
                    }
                    left_index = grid->GetCellIndex(i, j, k - 1);
                    right_index = grid->GetCellIndex(i, j, k);

                    r[0] = cell_center_coord[right_index].x() - cell_center_coord[left_index].x();
                    r[1] = cell_center_coord[right_index].y() - cell_center_coord[left_index].y();
                    r[2] = cell_center_coord[right_index].z() - cell_center_coord[left_index].z();
                    for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                    {
                        grad(0) = GetPrimGrad(left_index, iVal, 0);
                        grad(1) = GetPrimGrad(left_index, iVal, 1);
                        grad(2) = GetPrimGrad(left_index, iVal, 2);
                        riemann_para.prim_left(iVal) = m_prim[iVal][left_index] + 0.5 * GetLimiter(left_index, iVal) * grad.dot(r);
                        grad(0) = GetPrimGrad(right_index, iVal, 0);
                        grad(1) = GetPrimGrad(right_index, iVal, 1);
                        grad(2) = GetPrimGrad(right_index, iVal, 2);
                        riemann_para.prim_right(iVal) = m_prim[iVal][right_index] - 0.5 * GetLimiter(right_index, iVal) * grad.dot(r);
                    }
                    riemannSolver_->Solver(riemann_para);
                    for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                    {
                        GetResidual(iCell, iVal) -= riemann_para.flux(iVal) / jacobi;
                    }

                    if (isnan(GetResidual(iCell, 0)) || isnan(GetResidual(iCell, 1)) || isnan(GetResidual(iCell, 2)) || isnan(GetResidual(iCell, 3)) || isnan(GetResidual(iCell, 4)))
                    {
#ifdef USE_OMP

                        Log::error("i:{},j:{},k:{},cpu_index:{}", i, j, k, omp_get_thread_num());
#endif // DEBUG
                        Log::error("primL:{},{},{},{},{}", riemann_para.prim_left(0), riemann_para.prim_left(1), riemann_para.prim_left(2), riemann_para.prim_left(3), riemann_para.prim_left(4));
                        Log::error("primR:{},{},{},{},{}", riemann_para.prim_right(0), riemann_para.prim_right(1), riemann_para.prim_right(2), riemann_para.prim_right(3), riemann_para.prim_right(4));
                        Log::error("flux:{},{},{},{},{}", riemann_para.flux(0), riemann_para.flux(1), riemann_para.flux(2), riemann_para.flux(3), riemann_para.flux(4));
                        Log::error("norm:{},{},{}", riemann_para.norm(0), riemann_para.norm(1), riemann_para.norm(2));
                        exit(0);
                    }
                }
            }

        }
    }


    void Solver_NS_3D_Zaran::InviscidFluxMid()
    {
        auto& grid = GetGrid();
        auto& cell_topo = grid->GetCellTopo();
        auto& cell_type = cell_topo->GetType();
        auto& cell_center_coord = cell_topo->GetCenterCoord();
        auto& bound_patch = grid->GetBoundPatch();
        auto& bound_node_index = bound_patch.GetIndex();
        auto& bound_node_coord = bound_patch.GetCoordinate();

        // 起始点和终止点的编号,s: start, e: end
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        auto CellIndex = [&](int i, int j, int k) {return grid->GetCellIndex(i, j, k); };
        int iCell, i, j, k;
        DArray prim_left(5), prim_right(5);
        DArray prim_left_gradX(5), prim_right_gradX(5);
        DArray prim_left_gradY(5), prim_right_gradY(5);
        DArray prim_left_gradZ(5), prim_right_gradZ(5);
        DArray limiterCoef_left(5), limiterCoef_right(5);
        int left_index, right_index;
        DVector3D coord_left, coord_right, r, grad;
        RiemannSolverPara riemann_para;
        riemann_para.gamma_left = riemann_para.gamma_right = 1.4;
        for (int iPatch = 0;iPatch < bound_patch.GetPatchNum();iPatch++)
        {
            i = bound_node_index[iPatch][0];
            j = bound_node_index[iPatch][1];
            k = bound_node_index[iPatch][2];
            iCell = CellIndex(i, j, k);
            auto& jacobi = GetMetricJacob(iCell);
            // i direction
            riemann_para.norm(0) = GetMetricXi(iCell)[0];
            riemann_para.norm(1) = GetMetricXi(iCell)[1];
            riemann_para.norm(2) = GetMetricXi(iCell)[2];
            riemann_para.nt = GetMetricXi(iCell)[3];
            left_index = CellIndex(i, j, k);
            right_index = CellIndex(i + 1, j, k);
            coord_left = cell_center_coord[left_index];
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                prim_left[iVal] = m_prim[iVal][left_index];
                prim_left_gradX[iVal] = GetPrimGrad(left_index, iVal, 0);
                prim_left_gradY[iVal] = GetPrimGrad(left_index, iVal, 1);
                prim_left_gradZ[iVal] = GetPrimGrad(left_index, iVal, 2);
                limiterCoef_left[iVal] = GetLimiter(left_index, iVal);
            }
            if (cell_type[right_index] == CellType::Solid)
            {
                coord_right = bound_node_coord[iPatch];
                for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                {
                    prim_right[iVal] = m_prim_bound[iVal][iPatch];
                    prim_right_gradX[iVal] = GetPrimGradBound(iPatch, iVal, 0);
                    prim_right_gradY[iVal] = GetPrimGradBound(iPatch, iVal, 1);
                    prim_right_gradZ[iVal] = GetPrimGradBound(iPatch, iVal, 2);
                    limiterCoef_right[iVal] = GetLimiterBound(iPatch, iVal);
                }
            }
            else
            {
                coord_right = cell_center_coord[right_index];
                for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                {
                    prim_right[iVal] = m_prim[iVal][right_index];
                    prim_right_gradX[iVal] = GetPrimGrad(right_index, iVal, 0);
                    prim_right_gradY[iVal] = GetPrimGrad(right_index, iVal, 1);
                    prim_right_gradZ[iVal] = GetPrimGrad(right_index, iVal, 2);
                    limiterCoef_right[iVal] = GetLimiter(right_index, iVal);
                }
            }
            r[0] = coord_right.x() - coord_left.x();
            r[1] = coord_right.y() - coord_left.y();
            r[2] = coord_right.z() - coord_left.z();
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                grad(0) = prim_left_gradX[iVal];
                grad(1) = prim_left_gradY[iVal];
                grad(2) = prim_left_gradZ[iVal];
                riemann_para.prim_left(iVal) = prim_left[iVal] + 0.5 * limiterCoef_left[iVal] * grad.dot(r);
                grad(0) = prim_right_gradX[iVal];
                grad(1) = prim_right_gradY[iVal];
                grad(2) = prim_right_gradZ[iVal];
                riemann_para.prim_right(iVal) = prim_right[iVal] - 0.5 * limiterCoef_right[iVal] * grad.dot(r);
            }
            riemannSolver_->Solver(riemann_para);
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                GetResidual(iCell, iVal) += riemann_para.flux(iVal) / jacobi;
            }

            left_index = CellIndex(i - 1, j, k);
            right_index = CellIndex(i, j, k);
            coord_right = cell_center_coord[right_index];
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                prim_right[iVal] = m_prim[iVal][right_index];
                prim_right_gradX[iVal] = GetPrimGrad(right_index, iVal, 0);
                prim_right_gradY[iVal] = GetPrimGrad(right_index, iVal, 1);
                prim_right_gradZ[iVal] = GetPrimGrad(right_index, iVal, 2);
                limiterCoef_right[iVal] = GetLimiter(right_index, iVal);

            }
            if (cell_type[left_index] == CellType::Solid)
            {
                coord_left = bound_node_coord[iPatch];
                for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                {
                    prim_left[iVal] = m_prim_bound[iVal][iPatch];
                    prim_left_gradX[iVal] = GetPrimGradBound(iPatch, iVal, 0);
                    prim_left_gradY[iVal] = GetPrimGradBound(iPatch, iVal, 1);
                    prim_left_gradZ[iVal] = GetPrimGradBound(iPatch, iVal, 2);
                    limiterCoef_left[iVal] = GetLimiterBound(iPatch, iVal);
                }
            }
            else
            {
                coord_left = cell_center_coord[left_index];
                for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                {
                    prim_left[iVal] = m_prim[iVal][left_index];
                    prim_left_gradX[iVal] = GetPrimGrad(left_index, iVal, 0);
                    prim_left_gradY[iVal] = GetPrimGrad(left_index, iVal, 1);
                    prim_left_gradZ[iVal] = GetPrimGrad(left_index, iVal, 2);
                    limiterCoef_left[iVal] = GetLimiter(left_index, iVal);
                }
            }
            r[0] = coord_right.x() - coord_left.x();
            r[1] = coord_right.y() - coord_left.y();
            r[2] = coord_right.z() - coord_left.z();
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                grad(0) = prim_left_gradX[iVal];
                grad(1) = prim_left_gradY[iVal];
                grad(2) = prim_left_gradZ[iVal];
                riemann_para.prim_left(iVal) = prim_left[iVal] + 0.5 * limiterCoef_left[iVal] * grad.dot(r);
                grad(0) = prim_right_gradX[iVal];
                grad(1) = prim_right_gradY[iVal];
                grad(2) = prim_right_gradZ[iVal];
                riemann_para.prim_right(iVal) = prim_right[iVal] - 0.5 * limiterCoef_right[iVal] * grad.dot(r);
            }
            riemannSolver_->Solver(riemann_para);
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                GetResidual(iCell, iVal) -= riemann_para.flux(iVal) / jacobi;
            }
            // j direction
            riemann_para.norm(0) = GetMetricEta(iCell)[0];
            riemann_para.norm(1) = GetMetricEta(iCell)[1];
            riemann_para.norm(2) = GetMetricEta(iCell)[2];
            riemann_para.nt = GetMetricEta(iCell)[3];
            left_index = CellIndex(i, j, k);
            right_index = CellIndex(i, j + 1, k);
            coord_left = cell_center_coord[left_index];
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                prim_left[iVal] = m_prim[iVal][left_index];
                prim_left_gradX[iVal] = GetPrimGrad(left_index, iVal, 0);
                prim_left_gradY[iVal] = GetPrimGrad(left_index, iVal, 1);
                prim_left_gradZ[iVal] = GetPrimGrad(left_index, iVal, 2);
                limiterCoef_left[iVal] = GetLimiter(left_index, iVal);

            }
            if (cell_type[right_index] == CellType::Solid)
            {
                coord_right = bound_node_coord[iPatch];
                for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                {
                    prim_right[iVal] = m_prim_bound[iVal][iPatch];
                    prim_right_gradX[iVal] = GetPrimGradBound(iPatch, iVal, 0);
                    prim_right_gradY[iVal] = GetPrimGradBound(iPatch, iVal, 1);
                    prim_right_gradZ[iVal] = GetPrimGradBound(iPatch, iVal, 2);
                    limiterCoef_right[iVal] = GetLimiterBound(iPatch, iVal);
                }
            }
            else
            {
                coord_right = cell_center_coord[right_index];
                for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                {
                    prim_right[iVal] = m_prim[iVal][right_index];
                    prim_right_gradX[iVal] = GetPrimGrad(right_index, iVal, 0);
                    prim_right_gradY[iVal] = GetPrimGrad(right_index, iVal, 1);
                    prim_right_gradZ[iVal] = GetPrimGrad(right_index, iVal, 2);
                    limiterCoef_right[iVal] = GetLimiter(right_index, iVal);
                }
            }
            r[0] = coord_right.x() - coord_left.x();
            r[1] = coord_right.y() - coord_left.y();
            r[2] = coord_right.z() - coord_left.z();
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                grad(0) = prim_left_gradX[iVal];
                grad(1) = prim_left_gradY[iVal];
                grad(2) = prim_left_gradZ[iVal];
                riemann_para.prim_left(iVal) = prim_left[iVal] + 0.5 * limiterCoef_left[iVal] * grad.dot(r);
                grad(0) = prim_right_gradX[iVal];
                grad(1) = prim_right_gradY[iVal];
                grad(2) = prim_right_gradZ[iVal];
                riemann_para.prim_right(iVal) = prim_right[iVal] - 0.5 * limiterCoef_right[iVal] * grad.dot(r);
            }
            riemannSolver_->Solver(riemann_para);
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                GetResidual(iCell, iVal) += riemann_para.flux(iVal) / jacobi;
            }

            left_index = CellIndex(i, j - 1, k);
            right_index = CellIndex(i, j, k);
            coord_right = cell_center_coord[right_index];
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                prim_right[iVal] = m_prim[iVal][right_index];
                prim_right_gradX[iVal] = GetPrimGrad(right_index, iVal, 0);
                prim_right_gradY[iVal] = GetPrimGrad(right_index, iVal, 1);
                prim_right_gradZ[iVal] = GetPrimGrad(right_index, iVal, 2);
                limiterCoef_right[iVal] = GetLimiter(right_index, iVal);

            }
            if (cell_type[left_index] == CellType::Solid)
            {
                coord_left = bound_node_coord[iPatch];
                for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                {
                    prim_left[iVal] = m_prim_bound[iVal][iPatch];
                    prim_left_gradX[iVal] = GetPrimGradBound(iPatch, iVal, 0);
                    prim_left_gradY[iVal] = GetPrimGradBound(iPatch, iVal, 1);
                    prim_left_gradZ[iVal] = GetPrimGradBound(iPatch, iVal, 2);
                    limiterCoef_left[iVal] = GetLimiterBound(iPatch, iVal);
                }
            }
            else
            {
                coord_left = cell_center_coord[left_index];
                for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                {
                    prim_left[iVal] = m_prim[iVal][left_index];
                    prim_left_gradX[iVal] = GetPrimGrad(left_index, iVal, 0);
                    prim_left_gradY[iVal] = GetPrimGrad(left_index, iVal, 1);
                    prim_left_gradZ[iVal] = GetPrimGrad(left_index, iVal, 2);
                    limiterCoef_left[iVal] = GetLimiter(left_index, iVal);
                }
            }
            r[0] = coord_right.x() - coord_left.x();
            r[1] = coord_right.y() - coord_left.y();
            r[2] = coord_right.z() - coord_left.z();
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                grad(0) = prim_left_gradX[iVal];
                grad(1) = prim_left_gradY[iVal];
                grad(2) = prim_left_gradZ[iVal];
                riemann_para.prim_left(iVal) = prim_left[iVal] + 0.5 * limiterCoef_left[iVal] * grad.dot(r);
                grad(0) = prim_right_gradX[iVal];
                grad(1) = prim_right_gradY[iVal];
                grad(2) = prim_right_gradZ[iVal];
                riemann_para.prim_right(iVal) = prim_right[iVal] - 0.5 * limiterCoef_right[iVal] * grad.dot(r);
            }
            riemannSolver_->Solver(riemann_para);
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                GetResidual(iCell, iVal) -= riemann_para.flux(iVal) / jacobi;
            }
            // k direction
            riemann_para.norm(0) = GetMetricZeta(iCell)[0];
            riemann_para.norm(1) = GetMetricZeta(iCell)[1];
            riemann_para.norm(2) = GetMetricZeta(iCell)[2];
            riemann_para.nt = GetMetricZeta(iCell)[3];
            left_index = CellIndex(i, j, k);
            right_index = CellIndex(i, j, k + 1);
            coord_left = cell_center_coord[left_index];
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                prim_left[iVal] = m_prim[iVal][left_index];
                prim_left_gradX[iVal] = GetPrimGrad(left_index, iVal, 0);
                prim_left_gradY[iVal] = GetPrimGrad(left_index, iVal, 1);
                prim_left_gradZ[iVal] = GetPrimGrad(left_index, iVal, 2);
                limiterCoef_left[iVal] = GetLimiter(left_index, iVal);

            }
            if (cell_type[right_index] == CellType::Solid)
            {
                coord_right = bound_node_coord[iPatch];
                for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                {
                    prim_right[iVal] = m_prim_bound[iVal][iPatch];
                    prim_right_gradX[iVal] = GetPrimGradBound(iPatch, iVal, 0);
                    prim_right_gradY[iVal] = GetPrimGradBound(iPatch, iVal, 1);
                    prim_right_gradZ[iVal] = GetPrimGradBound(iPatch, iVal, 2);
                    limiterCoef_right[iVal] = GetLimiterBound(iPatch, iVal);
                }
            }
            else
            {
                coord_right = cell_center_coord[right_index];
                for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                {
                    prim_right[iVal] = m_prim[iVal][right_index];
                    prim_right_gradX[iVal] = GetPrimGrad(right_index, iVal, 0);
                    prim_right_gradY[iVal] = GetPrimGrad(right_index, iVal, 1);
                    prim_right_gradZ[iVal] = GetPrimGrad(right_index, iVal, 2);
                    limiterCoef_right[iVal] = GetLimiter(right_index, iVal);
                }
            }
            r[0] = coord_right.x() - coord_left.x();
            r[1] = coord_right.y() - coord_left.y();
            r[2] = coord_right.z() - coord_left.z();
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                grad(0) = prim_left_gradX[iVal];
                grad(1) = prim_left_gradY[iVal];
                grad(2) = prim_left_gradZ[iVal];
                riemann_para.prim_left(iVal) = prim_left[iVal] + 0.5 * limiterCoef_left[iVal] * grad.dot(r);
                grad(0) = prim_right_gradX[iVal];
                grad(1) = prim_right_gradY[iVal];
                grad(2) = prim_right_gradZ[iVal];
                riemann_para.prim_right(iVal) = prim_right[iVal] - 0.5 * limiterCoef_right[iVal] * grad.dot(r);
            }
            riemannSolver_->Solver(riemann_para);
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                GetResidual(iCell, iVal) += riemann_para.flux(iVal) / jacobi;
            }

            left_index = CellIndex(i, j, k - 1);
            right_index = CellIndex(i, j, k);
            coord_right = cell_center_coord[right_index];
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                prim_right[iVal] = m_prim[iVal][right_index];
                prim_right_gradX[iVal] = GetPrimGrad(right_index, iVal, 0);
                prim_right_gradY[iVal] = GetPrimGrad(right_index, iVal, 1);
                prim_right_gradZ[iVal] = GetPrimGrad(right_index, iVal, 2);
                limiterCoef_right[iVal] = GetLimiter(right_index, iVal);

            }
            if (cell_type[left_index] == CellType::Solid)
            {
                coord_left = bound_node_coord[iPatch];
                for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                {
                    prim_left[iVal] = m_prim_bound[iVal][iPatch];
                    prim_left_gradX[iVal] = GetPrimGradBound(iPatch, iVal, 0);
                    prim_left_gradY[iVal] = GetPrimGradBound(iPatch, iVal, 1);
                    prim_left_gradZ[iVal] = GetPrimGradBound(iPatch, iVal, 2);
                    limiterCoef_left[iVal] = GetLimiterBound(iPatch, iVal);
                }
            }
            else
            {
                coord_left = cell_center_coord[left_index];
                for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                {
                    prim_left[iVal] = m_prim[iVal][left_index];
                    prim_left_gradX[iVal] = GetPrimGrad(left_index, iVal, 0);
                    prim_left_gradY[iVal] = GetPrimGrad(left_index, iVal, 1);
                    prim_left_gradZ[iVal] = GetPrimGrad(left_index, iVal, 2);
                    limiterCoef_left[iVal] = GetLimiter(left_index, iVal);
                }
            }
            r[0] = coord_right.x() - coord_left.x();
            r[1] = coord_right.y() - coord_left.y();
            r[2] = coord_right.z() - coord_left.z();
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                grad(0) = prim_left_gradX[iVal];
                grad(1) = prim_left_gradY[iVal];
                grad(2) = prim_left_gradZ[iVal];
                riemann_para.prim_left(iVal) = prim_left[iVal] + 0.5 * limiterCoef_left[iVal] * grad.dot(r);
                grad(0) = prim_right_gradX[iVal];
                grad(1) = prim_right_gradY[iVal];
                grad(2) = prim_right_gradZ[iVal];
                riemann_para.prim_right(iVal) = prim_right[iVal] - 0.5 * limiterCoef_right[iVal] * grad.dot(r);
            }
            riemannSolver_->Solver(riemann_para);
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                GetResidual(iCell, iVal) -= riemann_para.flux(iVal) / jacobi;
            }
        }
    }

    void Solver_NS_3D_Zaran::ViscousFluxStruct()
    {
        //TODO
    }

    void Solver_NS_3D_Zaran::ViscousFluxMid()
    {
        //TODO
    }

    void Solver_NS_3D_Zaran::SourceFluxStruct()
    {
        //TODO
    }

    void Solver_NS_3D_Zaran::SourceFluxMid()
    {
        //TODO
    }

    void Solver_NS_3D_Zaran::BoundaryCondition()
    {
        auto& grid = GetGrid();
        BoundaryMapPtr& boundaryMapPtr = grid->GetBoundaryMap();
        auto& boundaryMap = boundaryMapPtr->GetBoundaryMap();
        auto& wallBound = boundaryMap["slipWall"];
        for (int iBound = 0; iBound < wallBound.size(); ++iBound)
            WallBC(wallBound[iBound]);
        auto& outletBound = boundaryMap["outlet"];
        for (int iBound = 0; iBound < outletBound.size(); ++iBound)
            OutletBC(outletBound[iBound]);
        auto& inletBound = boundaryMap["inlet"];
        for (int iBound = 0; iBound < inletBound.size(); ++iBound)
            InletBC(inletBound[iBound]);
        BoundPatchBC();
    }

    void Solver_NS_3D_Zaran::OutletBC(Boundary& bound)
    {
        int boundIndex = bound.GetIndex();
        int innerIndex = bound.GetInnerIndex();
        int ghost_index = bound.GetGhostIndex();
        m_prim[0][ghost_index] = m_prim[0][boundIndex];
        m_prim[1][ghost_index] = m_prim[1][boundIndex];
        m_prim[2][ghost_index] = m_prim[2][boundIndex];
        m_prim[3][ghost_index] = m_prim[3][boundIndex];
        m_prim[4][ghost_index] = m_prim[4][boundIndex];
        Prim2Cons(m_prim[0][ghost_index], m_prim[1][ghost_index], m_prim[2][ghost_index], m_prim[3][ghost_index], m_prim[4][ghost_index], GetCons(0, ghost_index), GetCons(1, ghost_index), GetCons(2, ghost_index), GetCons(3, ghost_index), GetCons(4, ghost_index));

    }

    void Solver_NS_3D_Zaran::InletBC(Boundary& bound)
    {
        FlowSolverPara* para = GetPara();
        int ghost_index = bound.GetGhostIndex();
        m_prim[0][ghost_index] = para->GetInflowDensity();
        m_prim[1][ghost_index] = para->GetInflowVelocityX();
        m_prim[2][ghost_index] = para->GetInflowVelocityY();
        m_prim[3][ghost_index] = para->GetInflowVelocityZ();
        m_prim[4][ghost_index] = para->GetInflowPressure();
        Prim2Cons(m_prim[0][ghost_index], m_prim[1][ghost_index], m_prim[2][ghost_index], m_prim[3][ghost_index], m_prim[4][ghost_index], GetCons(0, ghost_index), GetCons(1, ghost_index), GetCons(2, ghost_index), GetCons(3, ghost_index), GetCons(4, ghost_index));
    }

    void Solver_NS_3D_Zaran::WallBC(Boundary& bound)
    {
        int& inner_index = bound.GetInnerIndex();
        int bound_index = bound.GetIndex();
        auto& bound_norm = bound.GetNorm();
        m_prim[0][bound_index] = m_prim[0][inner_index];
        m_prim[1][bound_index] = m_prim[1][inner_index];
        m_prim[2][bound_index] = m_prim[2][inner_index];
        m_prim[3][bound_index] = m_prim[3][inner_index];
        m_prim[4][bound_index] = m_prim[4][inner_index];
        DVector3D inner_vel(m_prim[1][inner_index], m_prim[2][inner_index], m_prim[3][inner_index]);
        DVector3D bound_vel = inner_vel - (inner_vel.dot(bound_norm)) * bound_norm / (bound_norm.norm() * bound_norm.norm());
        m_prim[1][bound_index] = bound_vel(0);
        m_prim[2][bound_index] = bound_vel(1);
        m_prim[3][bound_index] = bound_vel(2);
        Prim2Cons(m_prim[0][bound_index], m_prim[1][bound_index], m_prim[2][bound_index], m_prim[3][bound_index], m_prim[4][bound_index], GetCons(0, bound_index), GetCons(1, bound_index), GetCons(2, bound_index), GetCons(3, bound_index), GetCons(4, bound_index));
    }

    void Solver_NS_3D_Zaran::BoundPatchBC()
    {
        auto& grid = GetGrid();
        auto& cell_topo = grid->GetCellTopo();
        auto& cell_type = cell_topo->GetType();
        auto& cell_center_coord = cell_topo->GetCenterCoord();
        auto& bound_patch = grid->GetBoundPatch();
        auto& bound_node_index = bound_patch.GetIndex();
        auto& bound_node_coord = bound_patch.GetCoordinate();
        auto& bound_norm = bound_patch.GetNormal();

        int i, j, k, iCell;
        auto CellIndex = [&](int i, int j, int k) {return grid->GetCellIndex(i, j, k); };
        DVector3D inner_vel, bound_vel;
        // #pragma omp parallel for private(iCell, i, j, k, inner_vel, bound_vel)
        for (int iPatch = 0;iPatch < bound_patch.GetPatchNum();iPatch++)
        {

            i = bound_node_index[iPatch][0];
            j = bound_node_index[iPatch][1];
            k = bound_node_index[iPatch][2];
            iCell = CellIndex(i, j, k);
            m_prim_bound[0][iPatch] = m_prim[0][iCell];
            m_prim_bound[1][iPatch] = m_prim[1][iCell];
            m_prim_bound[2][iPatch] = m_prim[2][iCell];
            m_prim_bound[3][iPatch] = m_prim[3][iCell];
            m_prim_bound[4][iPatch] = m_prim[4][iCell];
            inner_vel = { m_prim[1][iCell], m_prim[2][iCell], m_prim[3][iCell] };
            bound_vel = inner_vel - (inner_vel.dot(bound_norm[iPatch])) * bound_norm[iPatch] / (bound_norm[iPatch].norm() * bound_norm[iPatch].norm());
            m_prim_bound[1][iPatch] = bound_vel(0);
            m_prim_bound[2][iPatch] = bound_vel(1);
            m_prim_bound[3][iPatch] = bound_vel(2);
            m_prim_bound[1][iPatch] = bound_vel(0);
            m_prim_bound[2][iPatch] = bound_vel(1);
            m_prim_bound[3][iPatch] = bound_vel(2);

        }
    }


} // namespace zaran



