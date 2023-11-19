#include "Solver_NS_3D_Zaran.h"
#include "Vanleer.h"
namespace zaran
{

    void Solver_NS_3D_Zaran::CreateFieldData()
    {
        auto grid = GetGrid();
        int ni, nj, nk;
        grid->GetNodeNum(ni, nj, nk);
        int cell_num = (ni - 1) * (nj - 1) * (nk - 1);
        DArray empty_data(cell_num);
        auto& dataPtr = GetFieldData();
        auto& data = *dataPtr;
        //添加单元数据
        data.AddData("rho", empty_data);
        data.AddData("u", empty_data);
        data.AddData("v", empty_data);
        data.AddData("w", empty_data);
        data.AddData("p", empty_data);
        data.AddData("cons0", empty_data);
        data.AddData("cons1", empty_data);
        data.AddData("cons2", empty_data);
        data.AddData("cons3", empty_data);
        data.AddData("cons4", empty_data);
        data.AddData("res0", empty_data);
        data.AddData("res1", empty_data);
        data.AddData("res2", empty_data);
        data.AddData("res3", empty_data);
        data.AddData("res4", empty_data);
        data.AddData("dt", empty_data);
        data.AddData("limiterCoef0", empty_data);
        data.AddData("limiterCoef1", empty_data);
        data.AddData("limiterCoef2", empty_data);
        data.AddData("limiterCoef3", empty_data);
        data.AddData("limiterCoef4", empty_data);
        data.AddData("rhoGradX", empty_data);
        data.AddData("rhoGradY", empty_data);
        data.AddData("rhoGradZ", empty_data);
        data.AddData("uGradX", empty_data);
        data.AddData("uGradY", empty_data);
        data.AddData("uGradZ", empty_data);
        data.AddData("vGradX", empty_data);
        data.AddData("vGradY", empty_data);
        data.AddData("vGradZ", empty_data);
        data.AddData("wGradX", empty_data);
        data.AddData("wGradY", empty_data);
        data.AddData("wGradZ", empty_data);
        data.AddData("pGradX", empty_data);
        data.AddData("pGradY", empty_data);
        data.AddData("pGradZ", empty_data);
        data.AddData("consRK0", empty_data);
        data.AddData("consRK1", empty_data);
        data.AddData("consRK2", empty_data);
        data.AddData("consRK3", empty_data);
        data.AddData("consRK4", empty_data);
        data.AddData("coordTransXXi", empty_data);
        data.AddData("coordTransXEta", empty_data);
        data.AddData("coordTransXZeta", empty_data);
        data.AddData("coordTransXTau", empty_data);
        data.AddData("coordTransYXi", empty_data);
        data.AddData("coordTransYEta", empty_data);
        data.AddData("coordTransYZeta", empty_data);
        data.AddData("coordTransYTau", empty_data);
        data.AddData("coordTransZXi", empty_data);
        data.AddData("coordTransZEta", empty_data);
        data.AddData("coordTransZZeta", empty_data);
        data.AddData("coordTransZTau", empty_data);
        data.AddData("coordTransTXi", empty_data);
        data.AddData("coordTransTEta", empty_data);
        data.AddData("coordTransTZeta", empty_data);
        data.AddData("coordTransTTau", empty_data);
        data.AddData("coordTransXiX", empty_data);
        data.AddData("coordTransXiY", empty_data);
        data.AddData("coordTransXiZ", empty_data);
        data.AddData("coordTransXiT", empty_data);
        data.AddData("coordTransEtaX", empty_data);
        data.AddData("coordTransEtaY", empty_data);
        data.AddData("coordTransEtaZ", empty_data);
        data.AddData("coordTransEtaT", empty_data);
        data.AddData("coordTransZetaX", empty_data);
        data.AddData("coordTransZetaY", empty_data);
        data.AddData("coordTransZetaZ", empty_data);
        data.AddData("coordTransZetaT", empty_data);
        data.AddData("coordTransTauX", empty_data);
        data.AddData("coordTransTauY", empty_data);
        data.AddData("coordTransTauZ", empty_data);
        data.AddData("coordTransTauT", empty_data);
        data.AddData("coordTransJ", empty_data);
        //添加边界节点数据
        auto& bound_patch = grid->GetBoundPatch();
        int bound_node_num = bound_patch.GetPatchNum();
        empty_data.resize(bound_node_num);
        data.AddData("rhoB", empty_data);
        data.AddData("uB", empty_data);
        data.AddData("vB", empty_data);
        data.AddData("wB", empty_data);
        data.AddData("pB", empty_data);
        data.AddData("rhoGradXB", empty_data);
        data.AddData("rhoGradYB", empty_data);
        data.AddData("rhoGradZB", empty_data);
        data.AddData("uGradXB", empty_data);
        data.AddData("uGradYB", empty_data);
        data.AddData("uGradZB", empty_data);
        data.AddData("vGradXB", empty_data);
        data.AddData("vGradYB", empty_data);
        data.AddData("vGradZB", empty_data);
        data.AddData("wGradXB", empty_data);
        data.AddData("wGradYB", empty_data);
        data.AddData("wGradZB", empty_data);
        data.AddData("pGradXB", empty_data);
        data.AddData("pGradYB", empty_data);
        data.AddData("pGradZB", empty_data);
        data.AddData("limiterCoef0B", empty_data);
        data.AddData("limiterCoef1B", empty_data);
        data.AddData("limiterCoef2B", empty_data);
        data.AddData("limiterCoef3B", empty_data);
        data.AddData("limiterCoef4B", empty_data);

    }

    void Solver_NS_3D_Zaran::RegisterFieldData()
    {
        NSSolver::RegisterFieldData();
        auto data = GetFieldData();
        m_prim_bound.reserve(5);
        m_prim_bound_gradX.reserve(5);
        m_prim_bound_gradY.reserve(5);
        m_prim_bound_gradZ.reserve(5);
        m_limiter_bound.reserve(5);
        auto addDataToVector = [&data](Array<DArray*>& vec, const std::string& name) {
            vec.push_back(&data->GetData(name));
            };

        addDataToVector(m_prim_bound, "rhoB");
        addDataToVector(m_prim_bound, "uB");
        addDataToVector(m_prim_bound, "vB");
        addDataToVector(m_prim_bound, "wB");
        addDataToVector(m_prim_bound, "pB");
        addDataToVector(m_prim_bound_gradX, "rhoGradXB");
        addDataToVector(m_prim_bound_gradX, "uGradXB");
        addDataToVector(m_prim_bound_gradX, "vGradXB");
        addDataToVector(m_prim_bound_gradX, "wGradXB");
        addDataToVector(m_prim_bound_gradX, "pGradXB");
        addDataToVector(m_prim_bound_gradY, "rhoGradYB");
        addDataToVector(m_prim_bound_gradY, "uGradYB");
        addDataToVector(m_prim_bound_gradY, "vGradYB");
        addDataToVector(m_prim_bound_gradY, "wGradYB");
        addDataToVector(m_prim_bound_gradY, "pGradYB");
        addDataToVector(m_prim_bound_gradZ, "rhoGradZB");
        addDataToVector(m_prim_bound_gradZ, "uGradZB");
        addDataToVector(m_prim_bound_gradZ, "vGradZB");
        addDataToVector(m_prim_bound_gradZ, "wGradZB");
        addDataToVector(m_prim_bound_gradZ, "pGradZB");
        addDataToVector(m_limiter_bound, "limiterCoef0B");
        addDataToVector(m_limiter_bound, "limiterCoef1B");
        addDataToVector(m_limiter_bound, "limiterCoef2B");
        addDataToVector(m_limiter_bound, "limiterCoef3B");
        addDataToVector(m_limiter_bound, "limiterCoef4B");
    }
    void Solver_NS_3D_Zaran::InitField()
    {
        NSSolver::InitField();
        auto para = GetPara();
        DVector primInit = para->GetPrimitiveInflow();
        auto grid = GetGrid();
        int n_patch = grid->GetBoundPatch().GetPatchNum();
        for (int iPatch = 0;iPatch < n_patch;++iPatch)
        {
            (*m_prim_bound[0])[iPatch] = primInit[0];
            (*m_prim_bound[1])[iPatch] = primInit[1];
            (*m_prim_bound[2])[iPatch] = primInit[2];
            (*m_prim_bound[3])[iPatch] = primInit[3];
            (*m_prim_bound[4])[iPatch] = primInit[4];
        }
    }
    void Solver_NS_3D_Zaran::ComputeCoordTrans()
    {
        ComputeCoordTransStruct();
        ComputeCoordTransMid();
    }
    void Solver_NS_3D_Zaran::ComputeGradientWLS()
    {
        ComputeGradientWLSStruct();
        ComputeGradientWLSMid();
    }

    void Solver_NS_3D_Zaran::ComputeTimeStepLocal()
    {
        ComputeTimeStepLocalStruct();
        ComputeTimeStepLocalMid();
    }

    void Solver_NS_3D_Zaran::InviscidFlux()
    {
        InviscidFluxStruct();
        InviscidFluxMid();
    }

    void Solver_NS_3D_Zaran::ViscousFlux()
    {
        ViscousFluxStruct();
        ViscousFluxMid();
    }

    void Solver_NS_3D_Zaran::SourceFlux()
    {
        //TODO
    }

    void Solver_NS_3D_Zaran::ComputeLimiterCoefVK()
    {
        //TODO
    }

    void Solver_NS_3D_Zaran::ComputeLimiterCoefBJ()
    {
        //TODO
    }


    void Solver_NS_3D_Zaran::ComputeLimiterCoefNoLimiter()
    {
        //TODO
    }

    void Solver_NS_3D_Zaran::ComputeLimiterCoefOneOrder()
    {
        //TODO
    }

    void Solver_NS_3D_Zaran::ComputeBoundaryLimiterCoef()
    {
        //TODO
    }
    void Solver_NS_3D_Zaran::RungeKutta()
    {
        auto& grid = GetGrid();
        FlowSolverParaPtr para = GetPara();
        const DArray& rkCoef = para->GetRKCoef();
        int rkStage = rkCoef.size();
        auto& cons = m_Conservative;
        auto& cons_RK = m_ConservativeRK;
        auto& dt = *m_TimeStep;
        auto& res = m_Residual;
        auto& coordTrans = m_CoordTrans;
        int n_data = cons[0]->size();
        for (int iStage = 0; iStage < rkStage; ++iStage)
        {
            ComputeResidual();
#pragma omp parallel for
            for (int iNode = 0; iNode < n_data; ++iNode)
            {
                for (int iVal = 0; iVal < 5; ++iVal)
                {
                    auto& currentCons = (*cons[iVal])[iNode];
                    auto& currentCons_RK = (*cons_RK[iVal])[iNode];
                    auto& currentDt = dt[iNode];
                    auto& currentRes = (*res[iVal])[iNode];
                    currentCons = currentCons - rkCoef[iStage] * currentDt * currentRes * (*coordTrans[32])[iNode];
                }
            }
        }
    }
    void Solver_NS_3D_Zaran::ComputeCoordTransStruct()
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
        auto& coord_trans_coef = m_CoordTrans;
        CoordTrans coordTrans;
        for (int k = ks; k < ke - 1; k++)
        {
            for (int j = js; j < je - 1; j++)
            {
                for (int i = is; i < ie - 1; i++)
                {
                    iCell = CellIndex(i, j, k);
                    coordTrans.CalcCoordTrans(3, cell_center_coord[CellIndex(i + 1, j, k)], cell_center_coord[CellIndex(i - 1, j, k)], cell_center_coord[CellIndex(i, j + 1, k)], cell_center_coord[CellIndex(i, j - 1, k)], cell_center_coord[CellIndex(i, j, k + 1)], cell_center_coord[CellIndex(i, j, k - 1)]);
                    (*coord_trans_coef[0])[iCell] = coordTrans.GetX()[0];
                    (*coord_trans_coef[1])[iCell] = coordTrans.GetX()[1];
                    (*coord_trans_coef[2])[iCell] = coordTrans.GetX()[2];
                    (*coord_trans_coef[3])[iCell] = coordTrans.GetX()[3];
                    (*coord_trans_coef[4])[iCell] = coordTrans.GetY()[0];
                    (*coord_trans_coef[5])[iCell] = coordTrans.GetY()[1];
                    (*coord_trans_coef[6])[iCell] = coordTrans.GetY()[2];
                    (*coord_trans_coef[7])[iCell] = coordTrans.GetY()[3];
                    (*coord_trans_coef[8])[iCell] = coordTrans.GetZ()[0];
                    (*coord_trans_coef[9])[iCell] = coordTrans.GetZ()[1];
                    (*coord_trans_coef[10])[iCell] = coordTrans.GetZ()[2];
                    (*coord_trans_coef[11])[iCell] = coordTrans.GetZ()[3];
                    (*coord_trans_coef[12])[iCell] = coordTrans.GetT()[0];
                    (*coord_trans_coef[13])[iCell] = coordTrans.GetT()[1];
                    (*coord_trans_coef[14])[iCell] = coordTrans.GetT()[2];
                    (*coord_trans_coef[15])[iCell] = coordTrans.GetT()[3];
                    (*coord_trans_coef[16])[iCell] = coordTrans.GetXi()[0];
                    (*coord_trans_coef[17])[iCell] = coordTrans.GetXi()[1];
                    (*coord_trans_coef[18])[iCell] = coordTrans.GetXi()[2];
                    (*coord_trans_coef[19])[iCell] = coordTrans.GetXi()[3];
                    (*coord_trans_coef[20])[iCell] = coordTrans.GetEta()[0];
                    (*coord_trans_coef[21])[iCell] = coordTrans.GetEta()[1];
                    (*coord_trans_coef[22])[iCell] = coordTrans.GetEta()[2];
                    (*coord_trans_coef[23])[iCell] = coordTrans.GetEta()[3];
                    (*coord_trans_coef[24])[iCell] = coordTrans.GetZeta()[0];
                    (*coord_trans_coef[25])[iCell] = coordTrans.GetZeta()[1];
                    (*coord_trans_coef[26])[iCell] = coordTrans.GetZeta()[2];
                    (*coord_trans_coef[27])[iCell] = coordTrans.GetZeta()[3];
                    (*coord_trans_coef[28])[iCell] = coordTrans.GetTau()[0];
                    (*coord_trans_coef[29])[iCell] = coordTrans.GetTau()[1];
                    (*coord_trans_coef[30])[iCell] = coordTrans.GetTau()[2];
                    (*coord_trans_coef[31])[iCell] = coordTrans.GetTau()[3];
                    (*coord_trans_coef[32])[iCell] = coordTrans.J();
                }
            }
        }
    }

    void Solver_NS_3D_Zaran::ComputeCoordTransMid()
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
        auto& coord_trans_coef = m_CoordTrans;
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
            (*coord_trans_coef[0])[iCell] = coordTrans.GetX()[0];
            (*coord_trans_coef[1])[iCell] = coordTrans.GetX()[1];
            (*coord_trans_coef[2])[iCell] = coordTrans.GetX()[2];
            (*coord_trans_coef[3])[iCell] = coordTrans.GetX()[3];
            (*coord_trans_coef[4])[iCell] = coordTrans.GetY()[0];
            (*coord_trans_coef[5])[iCell] = coordTrans.GetY()[1];
            (*coord_trans_coef[6])[iCell] = coordTrans.GetY()[2];
            (*coord_trans_coef[7])[iCell] = coordTrans.GetY()[3];
            (*coord_trans_coef[8])[iCell] = coordTrans.GetZ()[0];
            (*coord_trans_coef[9])[iCell] = coordTrans.GetZ()[1];
            (*coord_trans_coef[10])[iCell] = coordTrans.GetZ()[2];
            (*coord_trans_coef[11])[iCell] = coordTrans.GetZ()[3];
            (*coord_trans_coef[12])[iCell] = coordTrans.GetT()[0];
            (*coord_trans_coef[13])[iCell] = coordTrans.GetT()[1];
            (*coord_trans_coef[14])[iCell] = coordTrans.GetT()[2];
            (*coord_trans_coef[15])[iCell] = coordTrans.GetT()[3];
            (*coord_trans_coef[16])[iCell] = coordTrans.GetXi()[0];
            (*coord_trans_coef[17])[iCell] = coordTrans.GetXi()[1];
            (*coord_trans_coef[18])[iCell] = coordTrans.GetXi()[2];
            (*coord_trans_coef[19])[iCell] = coordTrans.GetXi()[3];
            (*coord_trans_coef[20])[iCell] = coordTrans.GetEta()[0];
            (*coord_trans_coef[21])[iCell] = coordTrans.GetEta()[1];
            (*coord_trans_coef[22])[iCell] = coordTrans.GetEta()[2];
            (*coord_trans_coef[23])[iCell] = coordTrans.GetEta()[3];
            (*coord_trans_coef[24])[iCell] = coordTrans.GetZeta()[0];
            (*coord_trans_coef[25])[iCell] = coordTrans.GetZeta()[1];
            (*coord_trans_coef[26])[iCell] = coordTrans.GetZeta()[2];
            (*coord_trans_coef[27])[iCell] = coordTrans.GetZeta()[3];
            (*coord_trans_coef[28])[iCell] = coordTrans.GetTau()[0];
            (*coord_trans_coef[29])[iCell] = coordTrans.GetTau()[1];
            (*coord_trans_coef[30])[iCell] = coordTrans.GetTau()[2];
            (*coord_trans_coef[31])[iCell] = coordTrans.GetTau()[3];
            (*coord_trans_coef[32])[iCell] = coordTrans.J();
            if (isinf(coordTrans.J()))
                ZaranLog::error("i:{},j:{},k:{},J:{}", i, j, k, coordTrans.J());
        }
    }

    void Solver_NS_3D_Zaran::ComputeGradientWLSStruct()
    {
        //TODO
    }

    void Solver_NS_3D_Zaran::ComputeGradientWLSMid()
    {
        //TODO
    }

    void Solver_NS_3D_Zaran::ComputeTimeStepLocalStruct()
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
        auto& coord_trans_coef = m_CoordTrans;
        auto& rho = *m_Primitive[0];
        auto& u = *m_Primitive[1];
        auto& v = *m_Primitive[2];
        auto& w = *m_Primitive[3];
        auto& p = *m_Primitive[4];
        auto& coordTrans = m_CoordTrans;
        auto& dt = *m_TimeStep;
        auto para = GetPara();
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
                    double c = sqrt(gamma * p[iCell] / rho[iCell]);
                    double normXi = sqrt((*coordTrans[16])[iCell] * (*coordTrans[16])[iCell] + (*coordTrans[17])[iCell] * (*coordTrans[17])[iCell] + (*coordTrans[18])[iCell] * (*coordTrans[19])[iCell]);
                    double normEta = sqrt((*coordTrans[20])[iCell] * (*coordTrans[20])[iCell] + (*coordTrans[21])[iCell] * (*coordTrans[21])[iCell] + (*coordTrans[22])[iCell] * (*coordTrans[22])[iCell]);
                    double normZeta = sqrt((*coordTrans[24])[iCell] * (*coordTrans[24])[iCell] + (*coordTrans[25])[iCell] * (*coordTrans[25])[iCell] + (*coordTrans[26])[iCell] * (*coordTrans[26])[iCell]);
                    double uXi = u[iCell] * (*coordTrans[16])[iCell] + v[iCell] * (*coordTrans[17])[iCell] + w[iCell] * (*coordTrans[18])[iCell] + (*coordTrans[19])[iCell];
                    double uEta = u[iCell] * (*coordTrans[20])[iCell] + v[iCell] * (*coordTrans[21])[iCell] + w[iCell] * (*coordTrans[22])[iCell] + (*coordTrans[23])[iCell];
                    double uZeta = u[iCell] * (*coordTrans[24])[iCell] + v[iCell] * (*coordTrans[25])[iCell] + w[iCell] * (*coordTrans[26])[iCell] + (*coordTrans[27])[iCell];
                    double lambda = abs(uXi) + abs(uEta) + abs(uZeta) + c * (normXi + normEta + normZeta);
                    dt[iCell] = cfl / lambda;
                    min_dt = Min(min_dt, dt[iCell]);
                }
            }
        }
        GlobalData::Update("dt", min_dt);
    }

    void Solver_NS_3D_Zaran::ComputeTimeStepLocalMid()
    {
        //TODO
    }
    void Solver_NS_3D_Zaran::InviscidFluxStruct()
    {
        auto& grid = GetGrid();
        auto& cell_topo = grid->GetCellTopo();
        auto& cell_type = cell_topo->GetType();
        auto& cell_center_coord = cell_topo->GetCenterCoord();
        auto& prim = m_Primitive;
        auto& cons = m_Conservative;
        auto& primGradX = m_PrimGradX;
        auto& primGradY = m_PrimGradY;
        auto& primGradZ = m_PrimGradZ;
        auto& limiterCoef = m_LimiterCoef;
        auto& res = m_Residual;
        auto& coordTrans = m_CoordTrans;
        // 起始点和终止点的编号,s: start, e: end
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        int iCell;
        DVector3D r, grad;
        int left_index, right_index;
        RiemannSolverPara riemann_para;
        riemann_para.gammaL = riemann_para.gammaR = 1.4;
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
                    auto& jacobi = (*coordTrans[32])[iCell];
                    // i direction
                    riemann_para.norm(0) = (*coordTrans[16])[iCell];
                    riemann_para.norm(1) = (*coordTrans[17])[iCell];
                    riemann_para.norm(2) = (*coordTrans[18])[iCell];
                    riemann_para.nt = (*coordTrans[19])[iCell];
                    left_index = grid->GetCellIndex(i, j, k);
                    right_index = grid->GetCellIndex(i + 1, j, k);
                    r[0] = cell_center_coord[right_index].x() - cell_center_coord[left_index].x();
                    r[1] = cell_center_coord[right_index].y() - cell_center_coord[left_index].y();
                    r[2] = cell_center_coord[right_index].z() - cell_center_coord[left_index].z();
                    for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                    {
                        grad(0) = (*primGradX[iVal])[left_index];
                        grad(1) = (*primGradY[iVal])[left_index];
                        grad(2) = (*primGradZ[iVal])[left_index];
                        riemann_para.primL(iVal) = (*prim[iVal])[left_index] + 0.5 * (*limiterCoef[iVal])[left_index] * grad.dot(r);
                        grad(0) = (*primGradX[iVal])[right_index];
                        grad(1) = (*primGradY[iVal])[right_index];
                        grad(2) = (*primGradZ[iVal])[right_index];
                        riemann_para.primR(iVal) = (*prim[iVal])[right_index] - 0.5 * (*limiterCoef[iVal])[right_index] * grad.dot(r);
                    }
                    riemannSolver_->Solver(riemann_para);

                    for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                    {
                        (*res[iVal])[iCell] += riemann_para.flux(iVal) / jacobi;
                    }

                    left_index = grid->GetCellIndex(i - 1, j, k);
                    right_index = grid->GetCellIndex(i, j, k);

                    r[0] = cell_center_coord[right_index].x() - cell_center_coord[left_index].x();
                    r[1] = cell_center_coord[right_index].y() - cell_center_coord[left_index].y();
                    r[2] = cell_center_coord[right_index].z() - cell_center_coord[left_index].z();
                    for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                    {
                        grad(0) = (*primGradX[iVal])[left_index];
                        grad(1) = (*primGradY[iVal])[left_index];
                        grad(2) = (*primGradZ[iVal])[left_index];
                        riemann_para.primL(iVal) = (*prim[iVal])[left_index] + 0.5 * (*limiterCoef[iVal])[left_index] * grad.dot(r);
                        grad(0) = (*primGradX[iVal])[right_index];
                        grad(1) = (*primGradY[iVal])[right_index];
                        grad(2) = (*primGradZ[iVal])[right_index];
                        riemann_para.primR(iVal) = (*prim[iVal])[right_index] - 0.5 * (*limiterCoef[iVal])[right_index] * grad.dot(r);
                    }
                    riemannSolver_->Solver(riemann_para);
                    for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                    {
                        (*res[iVal])[iCell] -= riemann_para.flux(iVal) / jacobi;
                    }
                    // j direction
                    riemann_para.norm(0) = (*coordTrans[20])[iCell];
                    riemann_para.norm(1) = (*coordTrans[21])[iCell];
                    riemann_para.norm(2) = (*coordTrans[22])[iCell];
                    riemann_para.nt = (*coordTrans[23])[iCell];
                    left_index = grid->GetCellIndex(i, j, k);
                    right_index = grid->GetCellIndex(i, j + 1, k);
                    r[0] = cell_center_coord[right_index].x() - cell_center_coord[left_index].x();
                    r[1] = cell_center_coord[right_index].y() - cell_center_coord[left_index].y();
                    r[2] = cell_center_coord[right_index].z() - cell_center_coord[left_index].z();
                    for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                    {
                        grad(0) = (*primGradX[iVal])[left_index];
                        grad(1) = (*primGradY[iVal])[left_index];
                        grad(2) = (*primGradZ[iVal])[left_index];
                        riemann_para.primL(iVal) = (*prim[iVal])[left_index] + 0.5 * (*limiterCoef[iVal])[left_index] * grad.dot(r);
                        grad(0) = (*primGradX[iVal])[right_index];
                        grad(1) = (*primGradY[iVal])[right_index];
                        grad(2) = (*primGradZ[iVal])[right_index];
                        riemann_para.primR(iVal) = (*prim[iVal])[right_index] - 0.5 * (*limiterCoef[iVal])[right_index] * grad.dot(r);
                    }
                    riemannSolver_->Solver(riemann_para);
                    for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                    {
                        (*res[iVal])[iCell] += riemann_para.flux(iVal) / jacobi;
                    }
                    left_index = grid->GetCellIndex(i, j - 1, k);
                    right_index = grid->GetCellIndex(i, j, k);

                    r[0] = cell_center_coord[right_index].x() - cell_center_coord[left_index].x();
                    r[1] = cell_center_coord[right_index].y() - cell_center_coord[left_index].y();
                    r[2] = cell_center_coord[right_index].z() - cell_center_coord[left_index].z();
                    for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                    {
                        grad(0) = (*primGradX[iVal])[left_index];
                        grad(1) = (*primGradY[iVal])[left_index];
                        grad(2) = (*primGradZ[iVal])[left_index];
                        riemann_para.primL(iVal) = (*prim[iVal])[left_index] + 0.5 * (*limiterCoef[iVal])[left_index] * grad.dot(r);
                        grad(0) = (*primGradX[iVal])[right_index];
                        grad(1) = (*primGradY[iVal])[right_index];
                        grad(2) = (*primGradZ[iVal])[right_index];
                        riemann_para.primR(iVal) = (*prim[iVal])[right_index] - 0.5 * (*limiterCoef[iVal])[right_index] * grad.dot(r);
                    }
                    riemannSolver_->Solver(riemann_para);
                    for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                    {
                        (*res[iVal])[iCell] -= riemann_para.flux(iVal) / jacobi;
                    }
                    // k direction
                    riemann_para.norm(0) = (*coordTrans[24])[iCell];
                    riemann_para.norm(1) = (*coordTrans[25])[iCell];
                    riemann_para.norm(2) = (*coordTrans[26])[iCell];
                    riemann_para.nt = (*coordTrans[27])[iCell];
                    left_index = grid->GetCellIndex(i, j, k);
                    right_index = grid->GetCellIndex(i, j, k + 1);
                    r[0] = cell_center_coord[right_index].x() - cell_center_coord[left_index].x();
                    r[1] = cell_center_coord[right_index].y() - cell_center_coord[left_index].y();
                    r[2] = cell_center_coord[right_index].z() - cell_center_coord[left_index].z();
                    for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                    {
                        grad(0) = (*primGradX[iVal])[left_index];
                        grad(1) = (*primGradY[iVal])[left_index];
                        grad(2) = (*primGradZ[iVal])[left_index];
                        riemann_para.primL(iVal) = (*prim[iVal])[left_index] + 0.5 * (*limiterCoef[iVal])[left_index] * grad.dot(r);
                        grad(0) = (*primGradX[iVal])[right_index];
                        grad(1) = (*primGradY[iVal])[right_index];
                        grad(2) = (*primGradZ[iVal])[right_index];
                        riemann_para.primR(iVal) = (*prim[iVal])[right_index] - 0.5 * (*limiterCoef[iVal])[right_index] * grad.dot(r);
                    }
                    riemannSolver_->Solver(riemann_para);
                    for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                    {
                        (*res[iVal])[iCell] += riemann_para.flux(iVal) / jacobi;
                    }
                    left_index = grid->GetCellIndex(i, j, k - 1);
                    right_index = grid->GetCellIndex(i, j, k);

                    r[0] = cell_center_coord[right_index].x() - cell_center_coord[left_index].x();
                    r[1] = cell_center_coord[right_index].y() - cell_center_coord[left_index].y();
                    r[2] = cell_center_coord[right_index].z() - cell_center_coord[left_index].z();
                    for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                    {
                        grad(0) = (*primGradX[iVal])[left_index];
                        grad(1) = (*primGradY[iVal])[left_index];
                        grad(2) = (*primGradZ[iVal])[left_index];
                        riemann_para.primL(iVal) = (*prim[iVal])[left_index] + 0.5 * (*limiterCoef[iVal])[left_index] * grad.dot(r);
                        grad(0) = (*primGradX[iVal])[right_index];
                        grad(1) = (*primGradY[iVal])[right_index];
                        grad(2) = (*primGradZ[iVal])[right_index];
                        riemann_para.primR(iVal) = (*prim[iVal])[right_index] - 0.5 * (*limiterCoef[iVal])[right_index] * grad.dot(r);
                    }
                    riemannSolver_->Solver(riemann_para);
                    for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                    {
                        (*res[iVal])[iCell] -= riemann_para.flux(iVal) / jacobi;
                    }

                    if (isnan((*res[0])[iCell]) || isnan((*res[1])[iCell]) || isnan((*res[2])[iCell]) || isnan((*res[3])[iCell]) || isnan((*res[4])[iCell]))
                    {
#ifdef USE_OMP

                        ZaranLog::error("i:{},j:{},k:{},cpu_index:{}", i, j, k, omp_get_thread_num());
#endif // DEBUG
                        ZaranLog::error("primL:{},{},{},{},{}", riemann_para.primL(0), riemann_para.primL(1), riemann_para.primL(2), riemann_para.primL(3), riemann_para.primL(4));
                        ZaranLog::error("primR:{},{},{},{},{}", riemann_para.primR(0), riemann_para.primR(1), riemann_para.primR(2), riemann_para.primR(3), riemann_para.primR(4));
                        ZaranLog::error("flux:{},{},{},{},{}", riemann_para.flux(0), riemann_para.flux(1), riemann_para.flux(2), riemann_para.flux(3), riemann_para.flux(4));
                        ZaranLog::error("norm:{},{},{}", riemann_para.norm(0), riemann_para.norm(1), riemann_para.norm(2));
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

        auto& prim = m_Primitive;
        auto& cons = m_Conservative;
        auto& primGradX = m_PrimGradX;
        auto& primGradY = m_PrimGradY;
        auto& primGradZ = m_PrimGradZ;
        auto& limiterCoef = m_LimiterCoef;
        auto& res = m_Residual;
        auto& coordTrans = m_CoordTrans;

        auto& prim_bound = m_prim_bound;
        auto& prim_bound_gradX = m_prim_bound_gradX;
        auto& prim_bound_gradY = m_prim_bound_gradY;
        auto& prim_bound_gradZ = m_prim_bound_gradZ;
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
        riemann_para.gammaL = riemann_para.gammaR = 1.4;
        for (int iPatch = 0;iPatch < bound_patch.GetPatchNum();iPatch++)
        {
            i = bound_node_index[iPatch][0];
            j = bound_node_index[iPatch][1];
            k = bound_node_index[iPatch][2];
            iCell = CellIndex(i, j, k);
            auto& jacobi = (*coordTrans[32])[iCell];
            // i direction
            riemann_para.norm(0) = (*coordTrans[16])[iCell];
            riemann_para.norm(1) = (*coordTrans[17])[iCell];
            riemann_para.norm(2) = (*coordTrans[18])[iCell];
            riemann_para.nt = (*coordTrans[19])[iCell];
            left_index = CellIndex(i, j, k);
            right_index = CellIndex(i + 1, j, k);
            coord_left = cell_center_coord[left_index];
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                prim_left[iVal] = (*prim[iVal])[left_index];
                prim_left_gradX[iVal] = (*primGradX[iVal])[left_index];
                prim_left_gradY[iVal] = (*primGradY[iVal])[left_index];
                prim_left_gradZ[iVal] = (*primGradZ[iVal])[left_index];
                limiterCoef_left[iVal] = (*limiterCoef[iVal])[left_index];
            }
            if (cell_type[right_index] == CellType::Solid)
            {
                coord_right = bound_node_coord[iPatch];
                for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                {
                    prim_right[iVal] = (*prim_bound[iVal])[iPatch];
                    prim_right_gradX[iVal] = (*prim_bound_gradX[iVal])[iPatch];
                    prim_right_gradY[iVal] = (*prim_bound_gradY[iVal])[iPatch];
                    prim_right_gradZ[iVal] = (*prim_bound_gradZ[iVal])[iPatch];
                    limiterCoef_right[iVal] = (*m_limiter_bound[iVal])[iPatch];
                }
            }
            else
            {
                coord_right = cell_center_coord[right_index];
                for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                {
                    prim_right[iVal] = (*prim[iVal])[right_index];
                    prim_right_gradX[iVal] = (*primGradX[iVal])[right_index];
                    prim_right_gradY[iVal] = (*primGradY[iVal])[right_index];
                    prim_right_gradZ[iVal] = (*primGradZ[iVal])[right_index];
                    limiterCoef_right[iVal] = (*limiterCoef[iVal])[right_index];
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
                riemann_para.primL(iVal) = prim_left[iVal] + 0.5 * limiterCoef_left[iVal] * grad.dot(r);
                grad(0) = prim_right_gradX[iVal];
                grad(1) = prim_right_gradY[iVal];
                grad(2) = prim_right_gradZ[iVal];
                riemann_para.primR(iVal) = prim_right[iVal] - 0.5 * limiterCoef_right[iVal] * grad.dot(r);
            }
            riemannSolver_->Solver(riemann_para);
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                (*res[iVal])[iCell] += riemann_para.flux(iVal) / jacobi;
            }

            left_index = CellIndex(i - 1, j, k);
            right_index = CellIndex(i, j, k);
            coord_right = cell_center_coord[right_index];
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                prim_right[iVal] = (*prim[iVal])[right_index];
                prim_right_gradX[iVal] = (*primGradX[iVal])[right_index];
                prim_right_gradY[iVal] = (*primGradY[iVal])[right_index];
                prim_right_gradZ[iVal] = (*primGradZ[iVal])[right_index];
                limiterCoef_right[iVal] = (*limiterCoef[iVal])[right_index];

            }
            if (cell_type[left_index] == CellType::Solid)
            {
                coord_left = bound_node_coord[iPatch];
                for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                {
                    prim_left[iVal] = (*prim_bound[iVal])[iPatch];
                    prim_left_gradX[iVal] = (*prim_bound_gradX[iVal])[iPatch];
                    prim_left_gradY[iVal] = (*prim_bound_gradY[iVal])[iPatch];
                    prim_left_gradZ[iVal] = (*prim_bound_gradZ[iVal])[iPatch];
                    limiterCoef_left[iVal] = (*m_limiter_bound[iVal])[iPatch];
                }
            }
            else
            {
                coord_left = cell_center_coord[left_index];
                for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                {
                    prim_left[iVal] = (*prim[iVal])[left_index];
                    prim_left_gradX[iVal] = (*primGradX[iVal])[left_index];
                    prim_left_gradY[iVal] = (*primGradY[iVal])[left_index];
                    prim_left_gradZ[iVal] = (*primGradZ[iVal])[left_index];
                    limiterCoef_left[iVal] = (*limiterCoef[iVal])[left_index];
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
                riemann_para.primL(iVal) = prim_left[iVal] + 0.5 * limiterCoef_left[iVal] * grad.dot(r);
                grad(0) = prim_right_gradX[iVal];
                grad(1) = prim_right_gradY[iVal];
                grad(2) = prim_right_gradZ[iVal];
                riemann_para.primR(iVal) = prim_right[iVal] - 0.5 * limiterCoef_right[iVal] * grad.dot(r);
            }
            riemannSolver_->Solver(riemann_para);
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                (*res[iVal])[iCell] -= riemann_para.flux(iVal) / jacobi;
            }
            // j direction
            riemann_para.norm(0) = (*coordTrans[20])[iCell];
            riemann_para.norm(1) = (*coordTrans[21])[iCell];
            riemann_para.norm(2) = (*coordTrans[22])[iCell];
            riemann_para.nt = (*coordTrans[23])[iCell];
            left_index = CellIndex(i, j, k);
            right_index = CellIndex(i, j + 1, k);
            coord_left = cell_center_coord[left_index];
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                prim_left[iVal] = (*prim[iVal])[left_index];
                prim_left_gradX[iVal] = (*primGradX[iVal])[left_index];
                prim_left_gradY[iVal] = (*primGradY[iVal])[left_index];
                prim_left_gradZ[iVal] = (*primGradZ[iVal])[left_index];
                limiterCoef_left[iVal] = (*limiterCoef[iVal])[left_index];

            }
            if (cell_type[right_index] == CellType::Solid)
            {
                coord_right = bound_node_coord[iPatch];
                for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                {
                    prim_right[iVal] = (*prim_bound[iVal])[iPatch];
                    prim_right_gradX[iVal] = (*prim_bound_gradX[iVal])[iPatch];
                    prim_right_gradY[iVal] = (*prim_bound_gradY[iVal])[iPatch];
                    prim_right_gradZ[iVal] = (*prim_bound_gradZ[iVal])[iPatch];
                    limiterCoef_right[iVal] = (*m_limiter_bound[iVal])[iPatch];
                }
            }
            else
            {
                coord_right = cell_center_coord[right_index];
                for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                {
                    prim_right[iVal] = (*prim[iVal])[right_index];
                    prim_right_gradX[iVal] = (*primGradX[iVal])[right_index];
                    prim_right_gradY[iVal] = (*primGradY[iVal])[right_index];
                    prim_right_gradZ[iVal] = (*primGradZ[iVal])[right_index];
                    limiterCoef_right[iVal] = (*limiterCoef[iVal])[right_index];
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
                riemann_para.primL(iVal) = prim_left[iVal] + 0.5 * limiterCoef_left[iVal] * grad.dot(r);
                grad(0) = prim_right_gradX[iVal];
                grad(1) = prim_right_gradY[iVal];
                grad(2) = prim_right_gradZ[iVal];
                riemann_para.primR(iVal) = prim_right[iVal] - 0.5 * limiterCoef_right[iVal] * grad.dot(r);
            }
            riemannSolver_->Solver(riemann_para);
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                (*res[iVal])[iCell] += riemann_para.flux(iVal) / jacobi;
            }

            left_index = CellIndex(i, j - 1, k);
            right_index = CellIndex(i, j, k);
            coord_right = cell_center_coord[right_index];
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                prim_right[iVal] = (*prim[iVal])[right_index];
                prim_right_gradX[iVal] = (*primGradX[iVal])[right_index];
                prim_right_gradY[iVal] = (*primGradY[iVal])[right_index];
                prim_right_gradZ[iVal] = (*primGradZ[iVal])[right_index];
                limiterCoef_right[iVal] = (*limiterCoef[iVal])[right_index];

            }
            if (cell_type[left_index] == CellType::Solid)
            {
                coord_left = bound_node_coord[iPatch];
                for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                {
                    prim_left[iVal] = (*prim_bound[iVal])[iPatch];
                    prim_left_gradX[iVal] = (*prim_bound_gradX[iVal])[iPatch];
                    prim_left_gradY[iVal] = (*prim_bound_gradY[iVal])[iPatch];
                    prim_left_gradZ[iVal] = (*prim_bound_gradZ[iVal])[iPatch];
                    limiterCoef_left[iVal] = (*m_limiter_bound[iVal])[iPatch];
                }
            }
            else
            {
                coord_left = cell_center_coord[left_index];
                for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                {
                    prim_left[iVal] = (*prim[iVal])[left_index];
                    prim_left_gradX[iVal] = (*primGradX[iVal])[left_index];
                    prim_left_gradY[iVal] = (*primGradY[iVal])[left_index];
                    prim_left_gradZ[iVal] = (*primGradZ[iVal])[left_index];
                    limiterCoef_left[iVal] = (*limiterCoef[iVal])[left_index];
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
                riemann_para.primL(iVal) = prim_left[iVal] + 0.5 * limiterCoef_left[iVal] * grad.dot(r);
                grad(0) = prim_right_gradX[iVal];
                grad(1) = prim_right_gradY[iVal];
                grad(2) = prim_right_gradZ[iVal];
                riemann_para.primR(iVal) = prim_right[iVal] - 0.5 * limiterCoef_right[iVal] * grad.dot(r);
            }
            riemannSolver_->Solver(riemann_para);
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                (*res[iVal])[iCell] -= riemann_para.flux(iVal) / jacobi;
            }
            // k direction
            riemann_para.norm(0) = (*coordTrans[24])[iCell];
            riemann_para.norm(1) = (*coordTrans[25])[iCell];
            riemann_para.norm(2) = (*coordTrans[26])[iCell];
            riemann_para.nt = (*coordTrans[27])[iCell];
            left_index = CellIndex(i, j, k);
            right_index = CellIndex(i, j, k + 1);
            coord_left = cell_center_coord[left_index];
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                prim_left[iVal] = (*prim[iVal])[left_index];
                prim_left_gradX[iVal] = (*primGradX[iVal])[left_index];
                prim_left_gradY[iVal] = (*primGradY[iVal])[left_index];
                prim_left_gradZ[iVal] = (*primGradZ[iVal])[left_index];
                limiterCoef_left[iVal] = (*limiterCoef[iVal])[left_index];

            }
            if (cell_type[right_index] == CellType::Solid)
            {
                coord_right = bound_node_coord[iPatch];
                for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                {
                    prim_right[iVal] = (*prim_bound[iVal])[iPatch];
                    prim_right_gradX[iVal] = (*prim_bound_gradX[iVal])[iPatch];
                    prim_right_gradY[iVal] = (*prim_bound_gradY[iVal])[iPatch];
                    prim_right_gradZ[iVal] = (*prim_bound_gradZ[iVal])[iPatch];
                    limiterCoef_right[iVal] = (*m_limiter_bound[iVal])[iPatch];
                }
            }
            else
            {
                coord_right = cell_center_coord[right_index];
                for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                {
                    prim_right[iVal] = (*prim[iVal])[right_index];
                    prim_right_gradX[iVal] = (*primGradX[iVal])[right_index];
                    prim_right_gradY[iVal] = (*primGradY[iVal])[right_index];
                    prim_right_gradZ[iVal] = (*primGradZ[iVal])[right_index];
                    limiterCoef_right[iVal] = (*limiterCoef[iVal])[right_index];
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
                riemann_para.primL(iVal) = prim_left[iVal] + 0.5 * limiterCoef_left[iVal] * grad.dot(r);
                grad(0) = prim_right_gradX[iVal];
                grad(1) = prim_right_gradY[iVal];
                grad(2) = prim_right_gradZ[iVal];
                riemann_para.primR(iVal) = prim_right[iVal] - 0.5 * limiterCoef_right[iVal] * grad.dot(r);
            }
            riemannSolver_->Solver(riemann_para);
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                (*res[iVal])[iCell] += riemann_para.flux(iVal) / jacobi;
            }

            left_index = CellIndex(i, j, k - 1);
            right_index = CellIndex(i, j, k);
            coord_right = cell_center_coord[right_index];
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                prim_right[iVal] = (*prim[iVal])[right_index];
                prim_right_gradX[iVal] = (*primGradX[iVal])[right_index];
                prim_right_gradY[iVal] = (*primGradY[iVal])[right_index];
                prim_right_gradZ[iVal] = (*primGradZ[iVal])[right_index];
                limiterCoef_right[iVal] = (*limiterCoef[iVal])[right_index];

            }
            if (cell_type[left_index] == CellType::Solid)
            {
                coord_left = bound_node_coord[iPatch];
                for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                {
                    prim_left[iVal] = (*prim_bound[iVal])[iPatch];
                    prim_left_gradX[iVal] = (*prim_bound_gradX[iVal])[iPatch];
                    prim_left_gradY[iVal] = (*prim_bound_gradY[iVal])[iPatch];
                    prim_left_gradZ[iVal] = (*prim_bound_gradZ[iVal])[iPatch];
                    limiterCoef_left[iVal] = (*m_limiter_bound[iVal])[iPatch];
                }
            }
            else
            {
                coord_left = cell_center_coord[left_index];
                for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
                {
                    prim_left[iVal] = (*prim[iVal])[left_index];
                    prim_left_gradX[iVal] = (*primGradX[iVal])[left_index];
                    prim_left_gradY[iVal] = (*primGradY[iVal])[left_index];
                    prim_left_gradZ[iVal] = (*primGradZ[iVal])[left_index];
                    limiterCoef_left[iVal] = (*limiterCoef[iVal])[left_index];
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
                riemann_para.primL(iVal) = prim_left[iVal] + 0.5 * limiterCoef_left[iVal] * grad.dot(r);
                grad(0) = prim_right_gradX[iVal];
                grad(1) = prim_right_gradY[iVal];
                grad(2) = prim_right_gradZ[iVal];
                riemann_para.primR(iVal) = prim_right[iVal] - 0.5 * limiterCoef_right[iVal] * grad.dot(r);
            }
            riemannSolver_->Solver(riemann_para);
            for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
            {
                (*res[iVal])[iCell] -= riemann_para.flux(iVal) / jacobi;
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
            ComputeWallBC(wallBound[iBound]);
        auto& outletBound = boundaryMap["outlet"];
        for (int iBound = 0; iBound < outletBound.size(); ++iBound)
            ComputeOutletBC(outletBound[iBound]);
        auto& inletBound = boundaryMap["inlet"];
        for (int iBound = 0; iBound < inletBound.size(); ++iBound)
            ComputeInletBC(inletBound[iBound]);
        ComputeBoundPatchBC();
    }

    void Solver_NS_3D_Zaran::ComputeOutletBC(Boundary& bound)
    {
        auto& rho = *m_Primitive[0];
        auto& u = *m_Primitive[1];
        auto& v = *m_Primitive[2];
        auto& w = *m_Primitive[3];
        auto& p = *m_Primitive[4];
        auto& cons0 = *m_Conservative[0];
        auto& cons1 = *m_Conservative[1];
        auto& cons2 = *m_Conservative[2];
        auto& cons3 = *m_Conservative[3];
        auto& cons4 = *m_Conservative[4];
        int boundIndex = bound.GetIndex();
        int innerIndex = bound.GetInnerIndex();
        int ghost_index = bound.GetGhostIndex();
        rho[ghost_index] = rho[boundIndex];
        u[ghost_index] = u[boundIndex];
        v[ghost_index] = v[boundIndex];
        w[ghost_index] = w[boundIndex];
        p[ghost_index] = p[boundIndex];
        Primitive2Conservative(rho[ghost_index], u[ghost_index], v[ghost_index], w[ghost_index], p[ghost_index],
            cons0[ghost_index], cons1[ghost_index], cons2[ghost_index], cons3[ghost_index], cons4[ghost_index]);

    }

    void Solver_NS_3D_Zaran::ComputeInletBC(Boundary& bound)
    {
        FlowSolverParaPtr para = GetPara();
        auto& rho = *m_Primitive[0];
        auto& u = *m_Primitive[1];
        auto& v = *m_Primitive[2];
        auto& w = *m_Primitive[3];
        auto& p = *m_Primitive[4];
        auto& cons0 = *m_Conservative[0];
        auto& cons1 = *m_Conservative[1];
        auto& cons2 = *m_Conservative[2];
        auto& cons3 = *m_Conservative[3];
        auto& cons4 = *m_Conservative[4];
        int ghost_index = bound.GetGhostIndex();
        auto& inlet_para = para->GetPrimitiveInflow();
        rho[ghost_index] = inlet_para[0];
        u[ghost_index] = inlet_para[1];
        v[ghost_index] = inlet_para[2];
        w[ghost_index] = inlet_para[3];
        p[ghost_index] = inlet_para[4];
        Primitive2Conservative(rho[ghost_index], u[ghost_index], v[ghost_index], w[ghost_index], p[ghost_index],
            cons0[ghost_index], cons1[ghost_index], cons2[ghost_index], cons3[ghost_index], cons4[ghost_index]);
    }

    void Solver_NS_3D_Zaran::ComputeWallBC(Boundary& bound)
    {
        int& inner_index = bound.GetInnerIndex();
        int bound_index = bound.GetIndex();
        auto& rho = *m_Primitive[0];
        auto& u = *m_Primitive[1];
        auto& v = *m_Primitive[2];
        auto& w = *m_Primitive[3];
        auto& p = *m_Primitive[4];
        auto& cons0 = *m_Conservative[0];
        auto& cons1 = *m_Conservative[1];
        auto& cons2 = *m_Conservative[2];
        auto& cons3 = *m_Conservative[3];
        auto& cons4 = *m_Conservative[4];
        auto& bound_norm = bound.GetNorm();
        rho[bound_index] = rho[inner_index];
        u[bound_index] = u[inner_index];
        v[bound_index] = v[inner_index];
        w[bound_index] = w[inner_index];
        p[bound_index] = p[inner_index];
        DVector3D inner_vel(u[inner_index], v[inner_index], w[inner_index]);
        DVector3D bound_vel = inner_vel - (inner_vel.dot(bound_norm)) * bound_norm / (bound_norm.norm() * bound_norm.norm());
        u[bound_index] = bound_vel(0);
        v[bound_index] = bound_vel(1);
        w[bound_index] = bound_vel(2);
        Primitive2Conservative(rho[bound_index], u[bound_index], v[bound_index], w[bound_index], p[bound_index],
            cons0[bound_index], cons1[bound_index], cons2[bound_index], cons3[bound_index], cons4[bound_index]);
    }

    void Solver_NS_3D_Zaran::ComputeBoundPatchBC()
    {
        auto& grid = GetGrid();
        auto& cell_topo = grid->GetCellTopo();
        auto& cell_type = cell_topo->GetType();
        auto& cell_center_coord = cell_topo->GetCenterCoord();
        auto& bound_patch = grid->GetBoundPatch();
        auto& bound_node_index = bound_patch.GetIndex();
        auto& bound_node_coord = bound_patch.GetCoordinate();
        auto& bound_norm = bound_patch.GetNormal();

        auto& prim = m_Primitive;
        auto& cons = m_Conservative;
        auto& primGradX = m_PrimGradX;
        auto& primGradY = m_PrimGradY;
        auto& primGradZ = m_PrimGradZ;
        auto& limiterCoef = m_LimiterCoef;
        auto& res = m_Residual;
        auto& coordTrans = m_CoordTrans;

        auto& prim_bound = m_prim_bound;
        auto& prim_bound_gradX = m_prim_bound_gradX;
        auto& prim_bound_gradY = m_prim_bound_gradY;
        auto& prim_bound_gradZ = m_prim_bound_gradZ;
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
            // ZaranLog::info("iPatch:{},iCell:{},inner prim:{},{},{},{},{}", iPatch, iCell, (*prim[0])[iCell], (*prim[1])[iCell], (*prim[2])[iCell], (*prim[3])[iCell], (*prim[4])[iCell]);
            (*prim_bound[0])[iPatch] = (*prim[0])[iCell];
            (*prim_bound[1])[iPatch] = (*prim[1])[iCell];
            (*prim_bound[2])[iPatch] = (*prim[2])[iCell];
            (*prim_bound[3])[iPatch] = (*prim[3])[iCell];
            (*prim_bound[4])[iPatch] = (*prim[4])[iCell];
            inner_vel = { (*prim[1])[iCell],(*prim[2])[iCell],(*prim[3])[iCell] };
            bound_vel = inner_vel - (inner_vel.dot(bound_norm[iPatch])) * bound_norm[iPatch] / (bound_norm[iPatch].norm() * bound_norm[iPatch].norm());
            (*prim_bound[1])[iPatch] = bound_vel(0);
            (*prim_bound[2])[iPatch] = bound_vel(1);
            (*prim_bound[3])[iPatch] = bound_vel(2);
            (*prim_bound[1])[iPatch] = bound_vel(0);
            (*prim_bound[2])[iPatch] = bound_vel(1);
            (*prim_bound[3])[iPatch] = bound_vel(2);
        }
    }


} // namespace zaran



