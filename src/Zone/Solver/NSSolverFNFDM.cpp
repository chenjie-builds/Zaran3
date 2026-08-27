#include "NSSolverFNFDM.h"
#include "File.h"
#include "Log.h"
#include "ZaranError.h"
#include "MathBasic.h"
#include <fstream>
namespace zaran
{
NSSolverFNFDM::NSSolverFNFDM(index_type index, string name, shared_ptr<FlowSolverParam> para, shared_ptr<GridFN> grid,
                             shared_ptr<DataManagerNS> data_manager)
    : NSSolver(index, name, para, grid, data_manager)
{
    m_grid = grid;
    m_grad_wlsq = make_shared<GradWLSQ>(grid);
    m_node_metric = make_shared<Metric>(grid->GetTotalNodeNum());
}

NSSolverFNFDM::~NSSolverFNFDM()
{
}

void NSSolverFNFDM::InitFieldFarfield()
{
    auto grid = GetGrid();
    auto para = GetPara();
    auto data_manager = GetDataManager();
    double prim_far[5];
    prim_far[0] = para->GetInflowDensity();
    prim_far[1] = para->GetInflowVelocityX();
    prim_far[2] = para->GetInflowVelocityY();
    prim_far[3] = para->GetInflowVelocityZ();
    prim_far[4] = para->GetInflowPressure();
    int n_node = grid->GetTotalNodeNum();
    for (int iVal = 0; iVal < para->GetEqNum(); ++iVal)
    {
        for (int iNode = 0; iNode < n_node; ++iNode)
            data_manager->SetPrim(iVal, iNode, prim_far[iVal]);
    }
}
void NSSolverFNFDM::InitFieldFarFieldZeroVel()
{
    auto grid = GetGrid();
    auto para = GetPara();
    auto data_manager = GetDataManager();
    int node_num = grid->GetTotalNodeNum();
    for (int iNode = 0; iNode < node_num; ++iNode)
    {
        data_manager->SetPrim(0, iNode, para->GetInflowDensity());
        data_manager->SetPrim(1, iNode, 0.0);
        data_manager->SetPrim(2, iNode, 0.0);
        data_manager->SetPrim(3, iNode, 0.0);
        data_manager->SetPrim(4, iNode, para->GetInflowPressure());
    }
}
void NSSolverFNFDM::InitFieldBackup()
{
    auto grid = GetGrid();
    auto para = GetPara();
    auto data_manager = GetDataManager();
    std::string restart_file_name = para->GetBackupFieldFileName();
    std::ifstream fin(restart_file_name);
    if (!fin.is_open())
    {
        Log::warn("Backup file not found!");
        throw ZaranError("Backup field file not found: " + restart_file_name);
    }
    int node_num = grid->GetTotalNodeNum();
    double value;
    for (int iNode = 0; iNode < node_num; ++iNode)
    {
        for (int iVal = 0; iVal < para->GetEqNum(); ++iVal)
        {
            fin >> value;
            data_manager->SetPrim(iVal, iNode, value);
        }
    }
    fin.close();
}
void NSSolverFNFDM::InitFieldExplosion()
{
    auto grid = GetGrid();
    auto para = GetPara();
    auto data_manager = GetDataManager();
    double prim_far[5];
    prim_far[0] = para->GetInflowDensity();
    prim_far[1] = 0.0;
    prim_far[2] = 0.0;
    prim_far[3] = 0.0;
    prim_far[4] = para->GetInflowPressure();
    double explosion_pressure = GlobalData::GetDouble("init.explosion.pressure");
    double prim_explosion[5];
    prim_explosion[0] = prim_far[0];
    prim_explosion[1] = 0.0;
    prim_explosion[2] = 0.0;
    prim_explosion[3] = 0.0;
    prim_explosion[4] = explosion_pressure * prim_far[4];
    double explosion_center[3];
    explosion_center[0] = GlobalData::GetDouble("init.explosion.center_x");
    explosion_center[1] = GlobalData::GetDouble("init.explosion.center_y");
    explosion_center[2] = GlobalData::GetDouble("init.explosion.center_z");
    double explosion_radius = GlobalData::GetDouble("init.explosion.radius");
    double x, y, z;
    int n_node = grid->GetTotalNodeNum();
    for (int iNode = 0; iNode < n_node; ++iNode)
    {
        x = grid->GetNode().GetCoord(iNode)[0];
        y = grid->GetNode().GetCoord(iNode)[1];
        z = grid->GetNode().GetCoord(iNode)[2];
        double distance = sqrt((x - explosion_center[0]) * (x - explosion_center[0]) +
                               (y - explosion_center[1]) * (y - explosion_center[1]) +
                               (z - explosion_center[2]) * (z - explosion_center[2]));
        if (distance < explosion_radius)
        {
            for (int iVal = 0; iVal < para->GetEqNum(); ++iVal)
            {
                data_manager->SetPrim(iVal, iNode, prim_explosion[iVal]);
            }
        }
        else
        {
            for (int iVal = 0; iVal < para->GetEqNum(); ++iVal)
            {
                data_manager->SetPrim(iVal, iNode, prim_far[iVal]);
            }
        }
    }
}
void NSSolverFNFDM::CalcCoordTransCoef()
{
    Log::info("Compute NS 3D Coordination Transformation Coefficients");
    auto grid = GetGrid();
    auto &node = grid->GetNode();
    int node_num = grid->GetTotalNodeNum();
    double max_jacobian = -LARGE_NUMBER;
    double min_jacobian = LARGE_NUMBER;
    int max_jacobian_node = -1;
    int min_jacobian_node = -1;
    // #ifdef USE_OMP
    // #pragma omp parallel for reduction(max : max_jacobian) reduction(min : min_jacobian)
    // #endif // USE_OMP
    for (int iNode = 0; iNode < node_num; ++iNode)
    {
        if (node.GetType(iNode) != NodeType::inner)
        {
            continue;
        }
        auto neighbors = node.GetNeighborNode(iNode);
        auto xLeft = node.GetCoord(neighbors[0]);
        auto xRight = node.GetCoord(neighbors[1]);
        auto yLeft = node.GetCoord(neighbors[2]);
        auto yRight = node.GetCoord(neighbors[3]);
        if (grid->GetDim() == THREE_DIM)
        {
            auto zLeft = node.GetCoord(neighbors[4]);
            auto zRight = node.GetCoord(neighbors[5]);
            m_node_metric->CalcMetric(iNode, xRight.data(), xLeft.data(), yRight.data(), yLeft.data(), zRight.data(),
                                      zLeft.data());
        }
        else
        {
            m_node_metric->CalcMetric(iNode, xRight.data(), xLeft.data(), yRight.data(), yLeft.data());
        }
        double jacobian = m_node_metric->GetJacobian(iNode);
        if (jacobian > max_jacobian)
        {
            max_jacobian = jacobian;
            max_jacobian_node = iNode;
        }
        if (jacobian < min_jacobian)
        {
            min_jacobian = jacobian;
            min_jacobian_node = iNode;
        }
    }
    Log::info("NS 3D Coordination Transformation Coefficients are computed");
    Log::info("Max Jacobian = {} at node {}", max_jacobian, max_jacobian_node);
    Log::info("Min Jacobian = {} at node {}", min_jacobian, min_jacobian_node);
    auto neighbor = node.GetNeighborNode(min_jacobian_node);
    // Log::info("Neighbor node of min Jacobian node: ");
    // for (int i = 0; i < 6; ++i)
    //{
    //	Log::info("Node {} Coord: {} {} {}", neighbor[i], node.GetCoord(neighbor[i])[0], node.GetCoord(neighbor[i])[1],
    // node.GetCoord(neighbor[i])[2]);
    // }
    // neighbor = node.GetNeighborNode(max_jacobian_node);
    // Log::info("Neighbor node of max Jacobian node: {}, type: {}", max_jacobian_node,
    // (int)node.GetType(max_jacobian_node)); for (int i = 0; i < 6; ++i)
    //{
    //	Log::info("Node {} Coord: {} {} {}", neighbor[i], node.GetCoord(neighbor[i])[0], node.GetCoord(neighbor[i])[1],
    // node.GetCoord(neighbor[i])[2]);
    // }
}

void NSSolverFNFDM::BackupField(std::string &back_folder)
{
    auto grid = GetGrid();
    auto para = GetPara();
    auto data_manager = GetDataManager();
    string backup_file_name = para->GetBackupFieldFileName();
    backup_file_name = back_folder + "/" + backup_file_name;
    if (IsFileExist(backup_file_name))
    {
        RemoveFile(backup_file_name);
    }
    std::ofstream fout(backup_file_name);
    int n_node = grid->GetTotalNodeNum();
    for (int iNode = 0; iNode < n_node; ++iNode)
    {
        for (int iVal = 0; iVal < para->GetEqNum(); ++iVal)
        {
            fout << data_manager->GetPrim(iVal, iNode) << " ";
        }
        fout << std::endl;
    }
    fout.close();
}

void NSSolverFNFDM::RungeKutta()
{
    auto grid = GetGrid();
    auto data_manager = GetDataManager();
    auto para = GetPara();
    int rk_step = para->GetRkStep();
    auto &node = grid->GetNode();
    int n_node = grid->GetTotalNodeNum();
    for (int iVal = 0; iVal < 5; ++iVal)
    {
#ifdef USE_OMP
#pragma omp parallel for
#endif // USE_OMP
        for (int iNode = 0; iNode < n_node; ++iNode)
        {
            data_manager->SetConsOld(iVal, iNode, data_manager->GetCons(iVal, iNode));
        }
    }
    for (int iStep = 0; iStep < rk_step; ++iStep)
    {
        BoundaryCondition();
        CalcResidual();
        auto &rk_coef = para->GetRkCoef(iStep);
        for (int idx_eq = 0; idx_eq < 5; ++idx_eq)
        {
#ifdef USE_OMP
#pragma omp parallel for
#endif // USE_OMP
            for (int iNode = 0; iNode < n_node; ++iNode)
            {
                if (node.GetType(iNode) != NodeType::inner)
                    continue;
                double dt = data_manager->GetTimeStep(iNode);
                double jacobi = m_node_metric->GetJacobian(iNode);
                data_manager->SetCons(idx_eq, iNode,
                                      rk_coef[0] * data_manager->GetConsOld(idx_eq, iNode) +
                                          rk_coef[1] * data_manager->GetCons(idx_eq, iNode) +
                                          rk_coef[2] * dt * jacobi * data_manager->GetResidual(idx_eq, iNode));
            }
        }
    }
}

void NSSolverFNFDM::Prim2Cons()
{
    auto grid = GetGrid();
    auto gas = GetGas();
    auto data_manager = GetDataManager();
    int node_num = grid->GetTotalNodeNum();
    double prim[5], cons[5];
#ifdef USE_OMP
#pragma omp parallel for private(prim, cons)
#endif // USE_OMP
    for (int iNode = 0; iNode < node_num; ++iNode)
    {
        for (int iVal = 0; iVal < 5; ++iVal)
        {
            prim[iVal] = data_manager->GetPrim(iVal, iNode);
        }
        gas->Prim2Cons(prim, cons);
        data_manager->SetCons(iNode, cons);
    }
}

void NSSolverFNFDM::Cons2Prim()
{
    auto grid = GetGrid();
    auto data_manager = GetDataManager();
    int node_num = grid->GetTotalNodeNum();
    double prim[5], cons[5];
#ifdef USE_OMP
#pragma omp parallel for private(prim, cons)
#endif // USE_OMP
    for (int iNode = 0; iNode < node_num; ++iNode)
    {
        for (int iVal = 0; iVal < 5; ++iVal)
        {
            cons[iVal] = data_manager->GetCons(iVal, iNode);
        }
        GetGas()->Cons2Prim(cons, prim);
        data_manager->SetPrim(iNode, prim);
    }
}

void NSSolverFNFDM::MidPointReconstruct2ndOrder(int index_left, int index_right, double *value_rec_left,
                                                double *value_rec_right)
{
    int equ_num = GetPara()->GetEqNum();
    auto grid = GetGrid();
    NodeFN &node = grid->GetNode();
    auto data_manager = GetDataManager();
    auto coord_left = node.GetCoord(index_left);
    auto coord_right = node.GetCoord(index_right);
    for (int iVal = 0; iVal < equ_num; ++iVal)
    {
        value_rec_left[iVal] = data_manager->GetPrim(iVal, index_left) +
                               0.5 * data_manager->GetLimiter(iVal, index_left) *
                                   ((coord_right[0] - coord_left[0]) * data_manager->GetPrimGrad(iVal, 0, index_left) +
                                    (coord_right[1] - coord_left[1]) * data_manager->GetPrimGrad(iVal, 1, index_left) +
                                    (coord_right[2] - coord_left[2]) * data_manager->GetPrimGrad(iVal, 2, index_left));
        value_rec_right[iVal] =
            data_manager->GetPrim(iVal, index_right) -
            0.5 * data_manager->GetLimiter(iVal, index_right) *
                ((coord_right[0] - coord_left[0]) * data_manager->GetPrimGrad(iVal, 0, index_right) +
                 (coord_right[1] - coord_left[1]) * data_manager->GetPrimGrad(iVal, 1, index_right) +
                 (coord_right[2] - coord_left[2]) * data_manager->GetPrimGrad(iVal, 2, index_right));
    }
    // if (value_rec_left[0] < 0 || value_rec_left[4] < 0)
    //{
    //     for (int iVal = 0; iVal < equ_num; ++iVal)
    //     {
    //         value_rec_left[iVal] = data_manager->GetPrim(iVal, index_left);
    //     }
    // }
    // if (value_rec_right[0] < 0 || value_rec_right[4] < 0)
    //{
    //     for (int iVal = 0; iVal < equ_num; ++iVal)
    //     {
    //         value_rec_right[iVal] = data_manager->GetPrim(iVal, index_right);
    //     }
    // }
}

void NSSolverFNFDM::MidPointReconstruct1stOrder(int index_left, int index_right, double *value_rec_left,
                                                double *value_rec_right)
{
    int equ_num = GetPara()->GetEqNum();
    auto grid = GetGrid();
    NodeFN &node = grid->GetNode();
    auto data_manager = GetDataManager();
    for (int iVal = 0; iVal < equ_num; ++iVal)
    {
        value_rec_left[iVal] = data_manager->GetPrim(iVal, index_left);
        value_rec_right[iVal] = data_manager->GetPrim(iVal, index_right);
    }
}

void NSSolverFNFDM::BoundaryCondition()
{
    auto grid = GetGrid();
    BoundManagerFN &bound_map = grid->GetBoundaryMap();
    for (auto &boundary : bound_map.GetBoundMap())
    {
        auto &bound_name = boundary.first;
        auto &bound = boundary.second;
        if (bound_name == "hole")
            continue;
        if (bound_name == "riemann")
        {
            BCFarfield(bound);
        }
        else if (bound_name == "inlet")
        {
            BCInflow(bound);
        }
        else if (bound_name == "outlet")
        {
            BCOutflow(bound);
        }
        else if (bound_name == "wall")
        {
            BCWall(bound);
        }
    }
}
void NSSolverFNFDM::BCInlow(BoundFN &bound)
{
    auto para = GetPara();
    auto data_manager = GetDataManager();
    int bound_index = bound.GetIdxBound();
    double prim_far[5];
    prim_far[0] = para->GetInflowDensity();
    prim_far[1] = para->GetInflowVelocityX();
    prim_far[2] = para->GetInflowVelocityY();
    prim_far[3] = para->GetInflowVelocityZ();
    prim_far[4] = para->GetInflowPressure();
    double cons_far[5];
    GetGas()->Prim2Cons(prim_far, cons_far);
    for (int iVal = 0; iVal < 5; ++iVal)
    {
        data_manager->SetPrim(iVal, bound_index, prim_far[iVal]);
        data_manager->SetCons(iVal, bound_index, cons_far[iVal]);
    }
}
void NSSolverFNFDM::BCInflow(dynamic_array<BoundFN> &bound_list)
{
    auto para = GetPara();
    auto data_manager = GetDataManager();
    double prim_far[5];
    prim_far[0] = para->GetInflowDensity();
    prim_far[1] = para->GetInflowVelocityX();
    prim_far[2] = para->GetInflowVelocityY();
    prim_far[3] = para->GetInflowVelocityZ();
    prim_far[4] = para->GetInflowPressure();
    double cons_far[5];
    GetGas()->Prim2Cons(prim_far, cons_far);
    for (auto &bound : bound_list)
    {
        int bound_index = bound.GetIdxBound();
        for (int iVal = 0; iVal < 5; ++iVal)
        {
            data_manager->SetPrim(iVal, bound_index, prim_far[iVal]);
            data_manager->SetCons(iVal, bound_index, cons_far[iVal]);
        }
    }
}
void NSSolverFNFDM::BCOutflow(BoundFN &bound)
{
    auto data_manager = GetDataManager();
    int idx_bnd = bound.GetIdxBound();
    int inner_index = bound.GetIdxRef();
    double prim[5], cons[5];
    for (int iVal = 0; iVal < 5; ++iVal)
    {
        prim[iVal] = data_manager->GetPrim(iVal, inner_index);
    }
    GetGas()->Prim2Cons(prim, cons);
    for (int iVal = 0; iVal < 5; ++iVal)
    {
        data_manager->SetPrim(iVal, idx_bnd, prim[iVal]);
        data_manager->SetCons(iVal, idx_bnd, cons[iVal]);
    }
}

void NSSolverFNFDM::BCOutflow(dynamic_array<BoundFN> &bound_list)
{
    auto data_manager = GetDataManager();
    for (auto &bound : bound_list)
    {
        int idx_bnd = bound.GetIdxBound();
        int inner_index = bound.GetIdxRef();
        double prim[5], cons[5];
        for (int iVal = 0; iVal < 5; ++iVal)
        {
            prim[iVal] = data_manager->GetPrim(iVal, inner_index);
        }
        GetGas()->Prim2Cons(prim, cons);
        for (int iVal = 0; iVal < 5; ++iVal)
        {
            data_manager->SetPrim(iVal, idx_bnd, prim[iVal]);
            data_manager->SetCons(iVal, idx_bnd, cons[iVal]);
        }
    }
}

void NSSolverFNFDM::BCWall(BoundFN &bound)
{
    auto grid = GetGrid();
    auto node = grid->GetNode();
    auto data_manager = GetDataManager();
    int idx_ref = bound.GetIdxRef();
    int idx_bnd = bound.GetIdxBound();
    auto norm_bnd = bound.GetNormBound();
    auto coord_bnd = node.GetCoord(idx_bnd);
    auto coord_ref = node.GetCoord(idx_ref);
    double prim_bnd[5];
    for (int iVal = 0; iVal < 5; ++iVal)
    {
        prim_bnd[iVal] = data_manager->GetPrim(iVal, idx_ref);
        // for (int idx_dim = 0; idx_dim < 3; idx_dim++)
        //{
        //     prim_bnd[iVal] += (coord_bnd[idx_dim] - coord_ref[idx_dim]) *
        //                       data_manager->GetPrimGrad(iVal, idx_dim, idx_ref) *
        //                       data_manager->GetLimiter(iVal, idx_ref);
        // }
    }
    // if (prim_bnd[0] < 0)
    //     prim_bnd[0] = data_manager->GetPrim(0, idx_ref);
    // if (prim_bnd[4] < 0)
    //     prim_bnd[4] = data_manager->GetPrim(4, idx_ref);

    double vel_bnd[3] = {prim_bnd[1], prim_bnd[2], prim_bnd[3]};
    double vn = vel_bnd[0] * norm_bnd[0] + vel_bnd[1] * norm_bnd[1] + vel_bnd[2] * norm_bnd[2];
    prim_bnd[1] = prim_bnd[1] - 1.0 * vn * norm_bnd[0];
    prim_bnd[2] = prim_bnd[2] - 1.0 * vn * norm_bnd[1];
    prim_bnd[3] = prim_bnd[3] - 1.0 * vn * norm_bnd[2];
    double cons_bnd[5];
    GetGas()->Prim2Cons(prim_bnd, cons_bnd);
    for (int iVal = 0; iVal < 5; ++iVal)
    {
        data_manager->SetPrim(iVal, idx_bnd, prim_bnd[iVal]);
        data_manager->SetCons(iVal, idx_bnd, cons_bnd[iVal]);
    }
}
void NSSolverFNFDM::BCWall(dynamic_array<BoundFN> &bound_list)
{
    auto grid = GetGrid();
    auto data_manager = GetDataManager();
    for (size_t iBound = 0; iBound < bound_list.size(); ++iBound)
    {
        auto &bound = bound_list[iBound];
        int idx_ref = bound.GetIdxRef();
        int idx_bnd = bound.GetIdxBound();
        auto norm_bnd = bound.GetNormBound();
        double prim_bnd[5];
        for (int iVal = 0; iVal < 5; ++iVal)
        {
            prim_bnd[iVal] = data_manager->GetPrim(iVal, idx_ref);
        }
        double vel_bnd[3] = {prim_bnd[1], prim_bnd[2], prim_bnd[3]};
        double vn = vel_bnd[0] * norm_bnd[0] + vel_bnd[1] * norm_bnd[1] + vel_bnd[2] * norm_bnd[2];
        prim_bnd[1] = prim_bnd[1] - 1.0 * vn * norm_bnd[0];
        prim_bnd[2] = prim_bnd[2] - 1.0 * vn * norm_bnd[1];
        prim_bnd[3] = prim_bnd[3] - 1.0 * vn * norm_bnd[2];
        double cons_bnd[5];
        GetGas()->Prim2Cons(prim_bnd, cons_bnd);
        data_manager->SetPrim(idx_bnd, prim_bnd);
        data_manager->SetCons(idx_bnd, cons_bnd);
    }
}

void NSSolverFNFDM::BCFarfield(BoundFN &bound)
{
    auto data_manager = GetDataManager();
    int idx_ref = bound.GetIdxRef();
    int idx_bnd = bound.GetIdxBound();
    auto norm_bnd = bound.GetNormBound();
    double prim_in[5] = {data_manager->GetPrim(0, idx_ref), data_manager->GetPrim(1, idx_ref),
                         data_manager->GetPrim(2, idx_ref), data_manager->GetPrim(3, idx_ref),
                         data_manager->GetPrim(4, idx_ref)};
    double prim_bnd[5] = {data_manager->GetPrim(0, idx_bnd), data_manager->GetPrim(1, idx_bnd),
                          data_manager->GetPrim(2, idx_bnd), data_manager->GetPrim(3, idx_bnd),
                          data_manager->GetPrim(4, idx_bnd)};
    auto para = GetPara();
    auto gas = GetGas();
    double gamma = gas->GetGamma();
    double prim_far[5] = {para->GetInflowDensity(), para->GetInflowVelocityX(), para->GetInflowVelocityY(),
                          para->GetInflowVelocityZ(), para->GetInflowPressure()};
    double vel_in[3] = {prim_in[1], prim_in[2], prim_in[3]};
    // 速度在边界法向上的投影
    double vn_in = vel_in[0] * norm_bnd[0] + vel_in[1] * norm_bnd[1] + vel_in[2] * norm_bnd[2];
    double vel_far[3] = {para->GetInflowVelocityX(), para->GetInflowVelocityY(), para->GetInflowVelocityZ()};
    double vn_far = vel_far[0] * norm_bnd[0] + vel_far[1] * norm_bnd[1] + vel_far[2] * norm_bnd[2];
    double c_in = sqrt(gamma * prim_in[4] / prim_in[0]);
    double c_far = sqrt(gamma * para->GetInflowPressure() / para->GetInflowDensity());
    double mach = sqrt(vel_in[0] * vel_in[0] + vel_in[1] * vel_in[1] + vel_in[2] * vel_in[2]) / c_in;
    // 超声速出口
    if (mach >= 1.0)
    {
        if (vn_in >= 0) // 超声速出口
        {
            for (int iVal = 0; iVal < 5; ++iVal)
            {
                data_manager->SetPrim(iVal, idx_bnd, prim_in[iVal]);
            }
        }
        else // 超声速入口
        {
            data_manager->SetPrim(0, idx_bnd, para->GetInflowDensity());
            data_manager->SetPrim(1, idx_bnd, para->GetInflowVelocityX());
            data_manager->SetPrim(2, idx_bnd, para->GetInflowVelocityY());
            data_manager->SetPrim(3, idx_bnd, para->GetInflowVelocityZ());
            data_manager->SetPrim(4, idx_bnd, para->GetInflowPressure());
        }
    }
    else
    {
        double gamma1 = gamma - 1.0;
        double r_p = vn_in + 2.0 * c_in / gamma1;
        double r_m = vn_far - 2.0 * c_far / gamma1;
        double vn_bound = 0.5 * (r_p + r_m);
        double c_bound = 0.25 * gamma1 * (r_p - r_m);
        if (vn_bound <= 0) // 亚声速入口
        {
            double entropy = (prim_far[4] / pow(prim_far[0], gamma));
            data_manager->SetPrim(0, idx_bnd, pow(c_bound * c_bound / (entropy * gamma), 1.0 / gamma1));
            data_manager->SetPrim(1, idx_bnd, prim_far[1] + norm_bnd[0] * (vn_bound - vn_far));
            data_manager->SetPrim(2, idx_bnd, prim_far[2] + norm_bnd[1] * (vn_bound - vn_far));
            data_manager->SetPrim(3, idx_bnd, prim_far[3] + norm_bnd[2] * (vn_bound - vn_far));
            data_manager->SetPrim(4, idx_bnd, c_bound * c_bound * data_manager->GetPrim(0, idx_bnd) / gamma);
        }
        else // 亚声速出口
        {
            double entropy = prim_in[4] / pow(prim_in[0], gamma);
            data_manager->SetPrim(0, idx_bnd, pow(c_bound * c_bound / (entropy * gamma), 1.0 / gamma1));
            data_manager->SetPrim(1, idx_bnd, prim_in[1] + norm_bnd[0] * (vn_bound - vn_in));
            data_manager->SetPrim(2, idx_bnd, prim_in[2] + norm_bnd[1] * (vn_bound - vn_in));
            data_manager->SetPrim(3, idx_bnd, prim_in[3] + norm_bnd[2] * (vn_bound - vn_in));
            data_manager->SetPrim(4, idx_bnd, c_bound * c_bound * data_manager->GetPrim(0, idx_bnd) / gamma);
        }
    }
    for (int iVal = 0; iVal < 5; ++iVal)
        prim_bnd[iVal] = data_manager->GetPrim(iVal, idx_bnd);
    double cons_bound[5];
    GetGas()->Prim2Cons(prim_bnd, cons_bound);
    for (int iVal = 0; iVal < 5; ++iVal)
    {
        data_manager->SetCons(iVal, idx_bnd, cons_bound[iVal]);
    }
}
void NSSolverFNFDM::BCFarfield(dynamic_array<BoundFN> &bound_list)
{
    auto data_manager = GetDataManager();
    auto para = GetPara();
    auto gas = GetGas();
    double gamma = gas->GetGamma();
    double prim_far[5] = {para->GetInflowDensity(), para->GetInflowVelocityX(), para->GetInflowVelocityY(),
                          para->GetInflowVelocityZ(), para->GetInflowPressure()};
    for (auto &bound : bound_list)
    {
        int idx_ref = bound.GetIdxRef();
        int idx_bnd = bound.GetIdxBound();
        auto norm_bnd = bound.GetNormBound();
        double prim_in[5] = {data_manager->GetPrim(0, idx_ref), data_manager->GetPrim(1, idx_ref),
                             data_manager->GetPrim(2, idx_ref), data_manager->GetPrim(3, idx_ref),
                             data_manager->GetPrim(4, idx_ref)};
        double prim_bnd[5] = {data_manager->GetPrim(0, idx_bnd), data_manager->GetPrim(1, idx_bnd),
                              data_manager->GetPrim(2, idx_bnd), data_manager->GetPrim(3, idx_bnd),
                              data_manager->GetPrim(4, idx_bnd)};
        double vel_in[3] = {prim_in[1], prim_in[2], prim_in[3]};
        // 速度在边界法向上的投影
        double vn_in = vel_in[0] * norm_bnd[0] + vel_in[1] * norm_bnd[1] + vel_in[2] * norm_bnd[2];
        double vel_far[3] = {para->GetInflowVelocityX(), para->GetInflowVelocityY(), para->GetInflowVelocityZ()};
        double vn_far = vel_far[0] * norm_bnd[0] + vel_far[1] * norm_bnd[1] + vel_far[2] * norm_bnd[2];
        double c_in = sqrt(gamma * prim_in[4] / prim_in[0]);
        double c_far = sqrt(gamma * para->GetInflowPressure() / para->GetInflowDensity());
        double mach = sqrt(vel_in[0] * vel_in[0] + vel_in[1] * vel_in[1] + vel_in[2] * vel_in[2]) / c_in;
        // 超声速出口
        if (mach >= 1.0)
        {
            if (vn_in >= 0) // 超声速出口
            {
                for (int iVal = 0; iVal < 5; ++iVal)
                {
                    data_manager->SetPrim(iVal, idx_bnd, prim_in[iVal]);
                }
            }
            else // 超声速入口
            {
                data_manager->SetPrim(0, idx_bnd, para->GetInflowDensity());
                data_manager->SetPrim(1, idx_bnd, para->GetInflowVelocityX());
                data_manager->SetPrim(2, idx_bnd, para->GetInflowVelocityY());
                data_manager->SetPrim(3, idx_bnd, para->GetInflowVelocityZ());
                data_manager->SetPrim(4, idx_bnd, para->GetInflowPressure());
            }
        }
        else
        {
            double gamma1 = gamma - 1.0;
            double r_p = vn_in + 2.0 * c_in / gamma1;
            double r_m = vn_far - 2.0 * c_far / gamma1;
            double vn_bound = 0.5 * (r_p + r_m);
            double c_bound = 0.25 * gamma1 * (r_p - r_m);
            if (vn_bound <= 0) // 亚声速入口
            {
                double entropy = (prim_far[4] / pow(prim_far[0], gamma));
                data_manager->SetPrim(0, idx_bnd, pow(c_bound * c_bound / (entropy * gamma), 1.0 / gamma1));
                data_manager->SetPrim(1, idx_bnd, prim_far[1] + norm_bnd[0] * (vn_bound - vn_far));
                data_manager->SetPrim(2, idx_bnd, prim_far[2] + norm_bnd[1] * (vn_bound - vn_far));
                data_manager->SetPrim(3, idx_bnd, prim_far[3] + norm_bnd[2] * (vn_bound - vn_far));
                data_manager->SetPrim(4, idx_bnd, c_bound * c_bound * data_manager->GetPrim(0, idx_bnd) / gamma);
            }
            else // 亚声速出口
            {
                double entropy = prim_in[4] / pow(prim_in[0], gamma);
                data_manager->SetPrim(0, idx_bnd, pow(c_bound * c_bound / (entropy * gamma), 1.0 / gamma1));
                data_manager->SetPrim(1, idx_bnd, prim_in[1] + norm_bnd[0] * (vn_bound - vn_in));
                data_manager->SetPrim(2, idx_bnd, prim_in[2] + norm_bnd[1] * (vn_bound - vn_in));
                data_manager->SetPrim(3, idx_bnd, prim_in[3] + norm_bnd[2] * (vn_bound - vn_in));
                data_manager->SetPrim(4, idx_bnd, c_bound * c_bound * data_manager->GetPrim(0, idx_bnd) / gamma);
            }
        }
        for (int iVal = 0; iVal < 5; ++iVal)
            prim_bnd[iVal] = data_manager->GetPrim(iVal, idx_bnd);
        double cons_bound[5];
        GetGas()->Prim2Cons(prim_bnd, cons_bound);
        for (int iVal = 0; iVal < 5; ++iVal)
        {
            data_manager->SetCons(iVal, idx_bnd, cons_bound[iVal]);
        }
    }
}

void NSSolverFNFDM::CalcPrimGrad()
{
    auto para = GetPara();
    if (para->GetGradScheme() == GradScheme::WLS)
    {
        CalcGradWLS();
    }
    else if (para->GetGradScheme() == GradScheme::UFDM)
    {
        CalcGradFNFDM();
    }
    else if (para->GetGradScheme() == GradScheme::NONE)
    {
        NoGradient();
    }
    else
    {
        Log::warn("Unsupported Gradiend Scheme!");
    }
}
void NSSolverFNFDM::CalcPrimGradBound()
{
    int equ_num = GetPara()->GetEqNum();
    auto grid = GetGrid();
    auto &node = grid->GetNode();
    auto data_manager = GetDataManager();
    BoundManagerFN &bound_map = grid->GetBoundaryMap();
    auto &boundaryMap = bound_map.GetBoundMap();
    for (auto &boundary : boundaryMap)
    {
        auto &boundName = boundary.first;
        auto &bound = boundary.second;
        if (boundName == "hole")
            continue;
#ifdef USE_OMP
#pragma omp parallel for
#endif // USE_OMP
        for (int iVal = 0; iVal < equ_num; ++iVal)
        {
            for (auto &iter_bound : bound)
            {
                auto &idx_bound = iter_bound.GetIdxBound();
                auto &index_ref = iter_bound.GetIdxRef();
                for (int iter_dim = 0; iter_dim < 3; ++iter_dim)
                {
                    data_manager->SetPrimitiveGrad(iVal, iter_dim, idx_bound, 0.0);
                }
            }
        }
    }
}

void NSSolverFNFDM::CalcGradWLS()
{
    auto grid = GetGrid();
    auto node = grid->GetNode();
    auto data_manager = GetDataManager();
    int equ_num = GetPara()->GetEqNum();
    for (int iEqu = 0; iEqu < equ_num; ++iEqu)
    {
        m_grad_wlsq->CalcGradient(grid, data_manager->GetPrim(iEqu), data_manager->GetPrimitiveGrad(iEqu, 0),
                                  data_manager->GetPrimitiveGrad(iEqu, 1), data_manager->GetPrimitiveGrad(iEqu, 2));
    }
    CalcPrimGradBound();
}

void NSSolverFNFDM::CalcGradFNFDM()
{
    // TODO
}

void NSSolverFNFDM::NoGradient()
{
    // do nothing
}

void NSSolverFNFDM::CalcLimiter()
{
    LimiterType limiter_type = GetPara()->GetLimiterType();
    if (limiter_type != LimiterType::first_order)
    {
        int firstOrderSteps = GlobalData::GetInt("space.first_order_steps");
        int currentIter = GlobalData::GetInt("iteration.current_iter");
        if (currentIter < firstOrderSteps)
        {
            Log::info("First {}/{} iteration steps use 1st-order scheme", currentIter, firstOrderSteps);
            CalcLimiterFirstOrder();
        }
    }
    if (limiter_type == LimiterType::vk)
        CalcLimiterVK();
    else if (limiter_type == LimiterType::barth)
        CalcLimiterBJ();
    else if (limiter_type == LimiterType::none)
        CalcLimiterNone();
    else if (limiter_type == LimiterType::first_order)
        CalcLimiterFirstOrder();
    else
        Log::warn("Unsupported Limiter Type");
    // ComputeBoundaryLimiterCoef();
}

void NSSolverFNFDM::CalcLimiterVK()
{
    //CalcLimiterFirstOrder();
    //return;
    auto grid = GetGrid();
    NodeFN &node = grid->GetNode();
    auto para = GetPara();
    auto data_manager = GetDataManager();
    int node_num = grid->GetTotalNodeNum();
    int equ_num = GetPara()->GetEqNum();
    double max_limit = para->GetMaxLimit();
    double eps = 1e-6;
    double vk_coef = 10;
    double min_dis;
    for (int iVal = 0; iVal < equ_num; ++iVal)
    {
#ifdef USE_OMP
#pragma omp parallel for private(eps)
#endif // USE_OMP
        for (int iNode = 0; iNode < node_num; ++iNode)
        {
            if (node.GetType(iNode) != NodeType::inner && node.GetType(iNode) != NodeType::hole)
                continue;
            auto currentNodeCoord = node.GetCoord(iNode);
            auto neighborNode = node.GetNeighborNode(iNode);
            int nNeighbor = node.GetNeighborNodeNum(iNode);
            double current_val = data_manager->GetPrim(iVal, iNode);
            double maxVal = current_val;
            double minVal = current_val;
            min_dis = LARGE_NUMBER;
            for (int iNeighbor = 0; iNeighbor < nNeighbor; ++iNeighbor)
            {
                maxVal = Max(maxVal, data_manager->GetPrim(iVal, neighborNode[iNeighbor]));
                minVal = Min(minVal, data_manager->GetPrim(iVal, neighborNode[iNeighbor]));

                min_dis = Min(min_dis, DistanceOfTwoPoints(node.GetCoord(iNode).data(),
                                                           node.GetCoord(neighborNode[iNeighbor]).data()));
            }
            // if (node->GetType(iNode) == NodeType::inner)
            // {
            // 	auto jacobi = m_node_metric->GetJacobian(iNode);
            // 	eps = pow(venkatCoeff, 3) / jacobi;
            // }
            // else
            // {
            // eps = venkatCoeff * (maxVal - minVal) + SMALL_NUMBER;
            // eps = eps * eps;
            eps = pow(vk_coef * min_dis, 3);

            // }
            double delta_max = maxVal - current_val;
            double delta_min = minVal - current_val;
            double limite_coef = LARGE_NUMBER;
            for (int iNeighbor = 0; iNeighbor < nNeighbor; ++iNeighbor)
            {
                auto neighborNodeCoord = node.GetCoord(neighborNode[iNeighbor]);
                double delta2 = 0.0;
                for (int iDim = 0; iDim < 3; ++iDim)
                    delta2 += (neighborNodeCoord[iDim] - currentNodeCoord[iDim]) *
                              data_manager->GetPrimGrad(iVal, iDim, iNode);
                delta2 *= 0.5;
                if (delta2 > 0)
                {
                    limite_coef = LimiterVK(delta_max, delta2, eps);
                }
                else if (delta2 < 0)
                {
                    limite_coef = LimiterVK(delta_min, delta2, eps);
                }
                else
                {
                    limite_coef = 1.0;
                }
                limite_coef = Min(limite_coef, max_limit);
            }
            data_manager->SetLimiter(iVal, iNode, limite_coef);
        }
    }
}

void NSSolverFNFDM::CalcLimiterBJ()
{
    auto grid = GetGrid();
    NodeFN &node = grid->GetNode();
    auto para = GetPara();
    auto data_manager = GetDataManager();
    double max_limit = para->GetMaxLimit();
    int node_num = grid->GetTotalNodeNum();
    int equ_num = GetPara()->GetEqNum();
    double maxVal, minVal;
    for (int iVal = 0; iVal < equ_num; ++iVal)
    {
#ifdef USE_OMP
#pragma omp parallel for private(maxVal, minVal)
#endif // USE_OMP
        for (int iNode = 0; iNode < node_num; ++iNode)
        {
            if (node.GetType(iNode) != NodeType::inner && node.GetType(iNode) != NodeType::hole)
                continue;
            auto currentNodeCoord = node.GetCoord(iNode);
            auto neighborNode = node.GetNeighborNode(iNode);
            auto nNeighbor = node.GetNeighborNodeNum(iNode);
            maxVal = data_manager->GetPrim(iVal, iNode);
            minVal = data_manager->GetPrim(iVal, iNode);
            for (int iNeighbor = 0; iNeighbor < nNeighbor; ++iNeighbor)
            {
                maxVal = Max(maxVal, data_manager->GetPrim(iVal, neighborNode[iNeighbor]));
                minVal = Min(minVal, data_manager->GetPrim(iVal, neighborNode[iNeighbor]));
            }
            double deltaMax = maxVal - data_manager->GetPrim(iVal, iNode);
            double deltaMin = minVal - data_manager->GetPrim(iVal, iNode);
            double temp_limit = LARGE_NUMBER;
            data_manager->SetLimiter(iVal, iNode, LARGE_NUMBER);
            for (int iNeighbor = 0; iNeighbor < nNeighbor; ++iNeighbor)
            {
                auto neighborNodeCoord = node.GetCoord(neighborNode[iNeighbor]);
                double delta2 = 0.0;
                for (int iDim = 0; iDim < 3; ++iDim)
                    delta2 += (neighborNodeCoord[iDim] - currentNodeCoord[iDim]) *
                              data_manager->GetPrimGrad(iVal, iDim, iNode);
                delta2 *= 0.5;
                if (delta2 > 0)
                {
                    temp_limit = LimiterBarth(maxVal - data_manager->GetPrim(iVal, iNode), delta2);
                }
                else if (delta2 < 0)
                {
                    temp_limit = LimiterBarth(minVal - data_manager->GetPrim(iVal, iNode), delta2);
                }
                else
                {
                    temp_limit = 1.0;
                }
                temp_limit = Min(temp_limit, max_limit);
                data_manager->SetLimiter(iVal, iNode, Min(data_manager->GetLimiter(iVal, iNode), temp_limit));
            }
        }
    }
}

void NSSolverFNFDM::CalcLimiterNone()
{
    auto para = GetPara();
    auto grid = GetGrid();
    NodeFN &node = grid->GetNode();
    auto data_manager = GetDataManager();
    int node_num = grid->GetTotalNodeNum();
    for (int iVal = 0; iVal < para->GetEqNum(); ++iVal)
    {
#ifdef USE_OMP
#pragma omp parallel for
#endif // USE_OMP
        for (int iNode = 0; iNode < node_num; ++iNode)
        {
            //if (node.GetType(iNode) != NodeType::inner && node.GetType(iNode) != NodeType::hole)
            //    continue;
            data_manager->SetLimiter(iVal, iNode, 1.0);
        }
    }
}

void NSSolverFNFDM::CalcLimiterFirstOrder()
{
    auto para = GetPara();
    auto grid = GetGrid();
    NodeFN &node = grid->GetNode();
    auto data_manager = GetDataManager();
    int node_num = grid->GetTotalNodeNum();
    for (int iVal = 0; iVal < para->GetEqNum(); ++iVal)
    {
#ifdef USE_OMP
#pragma omp parallel for
#endif // USE_OMP
        for (int iNode = 0; iNode < node_num; ++iNode)
        {
            data_manager->SetLimiter(iVal, iNode, 0.0);
        }
    }
}

void NSSolverFNFDM::CalcLimiterBound()
{
    int equ_num = GetPara()->GetEqNum();
    auto grid = GetGrid();
    auto data_manager = GetDataManager();
    BoundManagerFN &boundaryMapPtr = grid->GetBoundaryMap();
    auto &boundaryMap = boundaryMapPtr.GetBoundMap();
    for (auto &boundary : boundaryMap)
    {
        auto &boundName = boundary.first;
        auto &bound = boundary.second;
        if (boundName == "hole")
            continue;
        for (int iVal = 0; iVal < equ_num; ++iVal)
        {
#ifdef USE_OMP
#pragma omp parallel for
#endif // USE_OMP
            for (int iBound = 0; iBound < bound.size(); ++iBound)
            {
                auto &bound_index = bound[iBound].GetIdxBound();
                auto &inner_index = bound[iBound].GetIdxRef();
                data_manager->SetLimiter(iVal, bound_index, 0.0);
            }
        }
    }
}

void NSSolverFNFDM::CheckPrimtive()
{
    return;
    auto grid = GetGrid();
    NodeFN &node = grid->GetNode();
    auto data_manager = GetDataManager();
    int node_num = grid->GetTotalNodeNum();
    int equ_num = GetPara()->GetEqNum();
    int nonphysical_node_num = 0;
    auto density = data_manager->GetPrim(ID_DENSITY);
    auto pressure = data_manager->GetPrim(ID_PRESSURE);
#ifdef USE_OMP
#pragma omp parallel for reduction(+ : nonphysical_node_num)
#endif // USE_OMP
    for (int iNode = 0; iNode < node_num; ++iNode)
    {
        bool exist_nonphysical = false;
        data_manager->SetNonPhysical(iNode, -1);
        if (node.GetType(iNode) != NodeType::inner && node.GetType(iNode) != NodeType::hole)
            continue;
        if (density[iNode] < 0 || pressure[iNode] < 0)
        {
            exist_nonphysical = true;
        }
        if (!exist_nonphysical)
        {
            for (int iVal = 0; iVal < equ_num; ++iVal)
            {
                if (std::isnan(data_manager->GetPrim(iVal, iNode)) || std::isinf(data_manager->GetPrim(iVal, iNode)))
                {
                    exist_nonphysical = true;
                    break;
                }
            }
        }
        if (exist_nonphysical)
        {
            data_manager->SetNonPhysical(iNode, 1);
            nonphysical_node_num++;
            Log::info("Step: {}, Non-physical Node: {}", GlobalData::GetInt("iteration.current_iter"), iNode);
            Log::warn("Non-physical Node: {}, neighbor num: {}, prim: {:6E},{:6E},{:6E},{:6E},{:6E}", iNode,
                      node.GetNeighborNodeNum(iNode), data_manager->GetPrim(0, iNode), data_manager->GetPrim(1, iNode),
                      data_manager->GetPrim(2, iNode), data_manager->GetPrim(3, iNode),
                      data_manager->GetPrim(4, iNode));
            Log::warn("Non-physical Node: {}, coord: {:6E}, {:6E}, {:6E}", iNode, node.GetCoord(iNode)[0],
                      node.GetCoord(iNode)[1], node.GetCoord(iNode)[2]);
            int neighbor_num = node.GetNeighborNodeNum(iNode);
            for (int iNeighbor = 0; iNeighbor < neighbor_num; ++iNeighbor)
            {
                Log::warn("Neighbor Node: {}, index:{}, neighbor_num: {}, Jacobian: {:6E},\n prim: "
                          "{:6E},{:6E},{:6E},{:6E},{:6E}\ndensity_grad: {:6E},{:6E},{:6E},\npressure_grad: "
                          "{:6E},{:6E},{:6E},\n limiter: {:6E},{:6E},{:6E},{:6E},{:6E}\n",
                          iNeighbor, node.GetNeighborNode(iNode)[iNeighbor],
                          node.GetNeighborNodeNum(node.GetNeighborNode(iNode)[iNeighbor]),
                          m_node_metric->GetJacobian(node.GetNeighborNode(iNode)[iNeighbor]),
                          data_manager->GetPrim(0, node.GetNeighborNode(iNode)[iNeighbor]),
                          data_manager->GetPrim(1, node.GetNeighborNode(iNode)[iNeighbor]),
                          data_manager->GetPrim(2, node.GetNeighborNode(iNode)[iNeighbor]),
                          data_manager->GetPrim(3, node.GetNeighborNode(iNode)[iNeighbor]),
                          data_manager->GetPrim(4, node.GetNeighborNode(iNode)[iNeighbor]),
                          data_manager->GetPrimGrad(0, 0, node.GetNeighborNode(iNode)[iNeighbor]),
                          data_manager->GetPrimGrad(0, 1, node.GetNeighborNode(iNode)[iNeighbor]),
                          data_manager->GetPrimGrad(0, 2, node.GetNeighborNode(iNode)[iNeighbor]),
                          data_manager->GetPrimGrad(4, 0, node.GetNeighborNode(iNode)[iNeighbor]),
                          data_manager->GetPrimGrad(4, 1, node.GetNeighborNode(iNode)[iNeighbor]),
                          data_manager->GetPrimGrad(4, 2, node.GetNeighborNode(iNode)[iNeighbor]),
                          data_manager->GetLimiter(0, node.GetNeighborNode(iNode)[iNeighbor]),
                          data_manager->GetLimiter(1, node.GetNeighborNode(iNode)[iNeighbor]),
                          data_manager->GetLimiter(2, node.GetNeighborNode(iNode)[iNeighbor]),
                          data_manager->GetLimiter(3, node.GetNeighborNode(iNode)[iNeighbor]),
                          data_manager->GetLimiter(4, node.GetNeighborNode(iNode)[iNeighbor]));
                // for (int jNeighbor = 0; jNeighbor <
                // node->GetNeighborNodeNum(node->GetNeighborNode(iNode)[iNeighbor]); jNeighbor++)
                // {
                // 	Log::warn("Neighbor Node: {}, Neighbor: {}, index: {}, Jacobian: {:6E}",
                // 			  iNeighbor, jNeighbor,
                // node->GetNeighborNode(node->GetNeighborNode(iNode)[iNeighbor])[jNeighbor]),
                // 		m_node_metric->GetJacobian(node->GetNeighborNode(node->GetNeighborNode(iNode)[iNeighbor])[jNeighbor]);
                // }
            }
        }
    }
    if (nonphysical_node_num > 0)
    {

        Log::warn("Non-physical Node Num: {}", nonphysical_node_num);
        auto para = GetPara();
        double cfl = para->GetCflNumber();
        para->ReduceCflNumber();
        Log::warn("CFL Number is reduced to {}", cfl);
    }
}

void NSSolverFNFDM::CheckResidual()
{
}

void NSSolverFNFDM::FixPrimtive()
{
    return;
    auto grid = GetGrid();
    NodeFN &node = grid->GetNode();
    auto data_manager = GetDataManager();
    int total_node_num = grid->GetTotalNodeNum();
    int equ_num = GetPara()->GetEqNum();
    dynamic_array<double> weight, distance;
    dynamic_array<int> physical_neighbor;
    double sum = 0;
#ifdef USE_OMP
#pragma omp parallel for private(physical_neighbor, weight, distance, sum)
#endif // USE_OMP
    for (int iNode = 0; iNode < total_node_num; ++iNode)
    {

        if (node.GetType(iNode) != NodeType::inner && node.GetType(iNode) != NodeType::hole)
            continue;
        if (data_manager->GetNonPhysical(iNode) == -1)
            continue;
        physical_neighbor.clear();
        auto current_node_coord = node.GetCoord(iNode);
        auto neighborNode = node.GetNeighborNode(iNode);
        int nNeighbor = node.GetNeighborNodeNum(iNode);
        for (int iNeighbor = 0; iNeighbor < nNeighbor; ++iNeighbor)
        {
            if (data_manager->GetNonPhysical(neighborNode[iNeighbor]) == 1)
                continue;
            physical_neighbor.push_back(neighborNode[iNeighbor]);
        }
        weight.resize(physical_neighbor.size());
        distance.resize(physical_neighbor.size());
        sum = 0;
        for (int iNeighbor = 0; iNeighbor < physical_neighbor.size(); ++iNeighbor)
        {
            if (data_manager->GetNonPhysical(physical_neighbor[iNeighbor]) == 1)
                continue;
            auto neighbor_node_coord = node.GetCoord(physical_neighbor[iNeighbor]);
            distance[iNeighbor] = DistanceOfTwoPoints(current_node_coord.data(), neighbor_node_coord.data());
            sum += 1.0 / distance[iNeighbor];
        }
        for (int iNeighbor = 0; iNeighbor < physical_neighbor.size(); ++iNeighbor)
        {
            if (data_manager->GetNonPhysical(physical_neighbor[iNeighbor]) == 1)
                continue;
            weight[iNeighbor] = 1.0 / (distance[iNeighbor] * sum);
        }
        for (int iVal = 0; iVal < equ_num; ++iVal)
        {
            if (std::isnan(data_manager->GetPrim(iVal, iNode)) || std::isinf(data_manager->GetPrim(iVal, iNode)))
            {
                data_manager->SetPrim(iVal, iNode, 0);
                for (int iNeighbor = 0; iNeighbor < physical_neighbor.size(); ++iNeighbor)
                {
                    if (data_manager->GetNonPhysical(physical_neighbor[iNeighbor]) == 1)
                        continue;
                    data_manager->SetPrim(
                        iVal, iNode, data_manager->GetPrim(iVal, physical_neighbor[iNeighbor]) * weight[iNeighbor]);
                }
            }
        }
        if (data_manager->GetPrim(0, iNode) < 0)
            data_manager->SetPrim(0, iNode, 0);
        for (int iNeighbor = 0; iNeighbor < physical_neighbor.size(); ++iNeighbor)
        {
            if (data_manager->GetNonPhysical(physical_neighbor[iNeighbor]) == 1)
                continue;
            data_manager->SetPrim(0, iNode,
                                  data_manager->GetPrim(0, iNode) +
                                      data_manager->GetPrim(0, physical_neighbor[iNeighbor]) * weight[iNeighbor]);
        }
        if (data_manager->GetPrim(0, iNode) < 0)
            data_manager->SetPrim(0, iNode, 0);
        for (int iNeighbor = 0; iNeighbor < physical_neighbor.size(); ++iNeighbor)
        {
            if (data_manager->GetNonPhysical(physical_neighbor[iNeighbor]) == 1)
                continue;
            data_manager->SetPrim(4, iNode,
                                  data_manager->GetPrim(1, iNode) +
                                      data_manager->GetPrim(4, physical_neighbor[iNeighbor]) * weight[iNeighbor]);
        }
        double prim[5] = {data_manager->GetPrim(0, iNode), data_manager->GetPrim(1, iNode),
                          data_manager->GetPrim(2, iNode), data_manager->GetPrim(3, iNode),
                          data_manager->GetPrim(4, iNode)};
        double cons[5];
        GetGas()->Prim2Cons(prim, cons);
        for (int iVal = 0; iVal < 5; ++iVal)
        {
            data_manager->SetCons(iVal, iNode, cons[iVal]);
            data_manager->SetResidual(iVal, iNode, 0.0);
        }
    }
}

void NSSolverFNFDM::Preprocess()
{
    NSSolver::Preprocess();
    CalcPrimGrad();
    CalcLimiter();
}

void NSSolverFNFDM::CalcTimeStepLocal()
{
    auto grid = GetGrid();
    auto node = grid->GetNode();
    auto para = GetPara();
    auto data_manager = GetDataManager();
    double gamma = para->GetGas()->GetGamma();
    double cfl = para->GetCflNumber();
    int inner_node_num = grid->GetInnerNodeNum();
    auto &inner_node = grid->GetInnerNode();
    auto density = data_manager->GetPrim(ID_DENSITY);
    auto pressure = data_manager->GetPrim(ID_PRESSURE);
    auto velocity_x = data_manager->GetPrim(ID_VELOCITY_X);
    auto velocity_y = data_manager->GetPrim(ID_VELOCITY_Y);
    auto velocity_z = data_manager->GetPrim(ID_VELOCITY_Z);
#ifdef USE_OMP
#pragma omp parallel for
#endif // USE_OMP
    for (int iNode = 0; iNode < inner_node_num; ++iNode)
    {
        int idx = inner_node[iNode];
        auto xi = m_node_metric->GetXi(idx);
        auto eta = m_node_metric->GetEta(idx);
        auto zeta = m_node_metric->GetZeta(idx);
        auto jacobi = m_node_metric->GetJacobian(idx);
        double c = sqrt(gamma * pressure[idx] / density[idx]);
        double norm_xi = sqrt(xi[0] * xi[0] + xi[1] * xi[1] + xi[2] * xi[2]);
        double norm_eta = sqrt(eta[0] * eta[0] + eta[1] * eta[1] + eta[2] * eta[2]);
        double norm_zeta = sqrt(zeta[0] * zeta[0] + zeta[1] * zeta[1] + zeta[2] * zeta[2]);
        double u_xi = velocity_x[idx] * xi[0] + velocity_y[idx] * xi[1] + velocity_z[idx] * xi[2];
        double u_eta = velocity_x[idx] * eta[0] + velocity_y[idx] * eta[1] + velocity_z[idx] * eta[2];
        double u_zeta = velocity_x[idx] * zeta[0] + velocity_y[idx] * zeta[1] + velocity_z[idx] * zeta[2];
        double lamda = abs(u_xi) + abs(u_eta) + abs(u_zeta) + c * (norm_xi + norm_eta + norm_zeta);
        lamda = lamda * jacobi;
        data_manager->SetTimeStep(idx, cfl / lamda);
    }
}
void NSSolverFNFDM::CalcMinTimeStep(double &dt)
{
    auto grid = GetGrid();
    auto node = grid->GetNode();
    auto para = GetPara();
    auto data_manager = GetDataManager();
    double gamma = para->GetGas()->GetGamma();
    double cfl = para->GetCflNumber();
    int inner_node_num = grid->GetInnerNodeNum();
    auto &inner_node = grid->GetInnerNode();
    double min_dt = LARGE_NUMBER;
#ifdef USE_OMP
#pragma omp parallel for reduction(min : min_dt)
#endif // USE_OMP
    for (int iNode = 0; iNode < inner_node_num; ++iNode)
    {
        int idx = inner_node[iNode];
        if (min_dt > data_manager->GetTimeStep(idx))
        {
            min_dt = data_manager->GetTimeStep(idx);
        }
    }
    dt = min_dt;
}

void NSSolverFNFDM::ReduceTimeStep(double &dt)
{
    auto grid = GetGrid();
    auto para = GetPara();
    auto data_manager = GetDataManager();
#ifdef USE_OMP
#pragma omp parallel for
#endif // USE_OMP
    for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
    {
        data_manager->SetTimeStep(iNode, dt);
    }
}

void NSSolverFNFDM::ZeroResidual()
{
    int equ_num = GetPara()->GetEqNum();
    auto grid = GetGrid();
    auto node = grid->GetNode();
    auto data_manager = GetDataManager();
    int node_num = grid->GetTotalNodeNum();
    double zero_residual[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
#ifdef USE_OMP
#pragma omp parallel for
#endif // USE_OMP
    for (int iNode = 0; iNode < node_num; ++iNode)
    {
        data_manager->SetResidual(iNode, zero_residual);
    }
}

void NSSolverFNFDM::CalcConvectionResidual()
{
    auto grid = GetGrid();
    NodeFN &node = grid->GetNode();
    auto data_manager = GetDataManager();
    int equ_num = GetPara()->GetEqNum();
    int inner_node_num = grid->GetInnerNodeNum();
    auto &inner_node = grid->GetInnerNode();
    RiemannSolverPara riemann_para[2];
    double res_tmp[5];
//#ifdef USE_OMP
//#pragma omp parallel for private(riemann_para, res_tmp)
//#endif // USE_OMP
    for (int iNode = 0; iNode < inner_node_num; ++iNode)
    {
        for (int i = 0; i < 2; ++i)
        {
            riemann_para[i].gamma_left = riemann_para[i].gamma_right = GetGas()->GetGamma();
        }

        bool exist_negative = false;
        int idx = inner_node[iNode];
        auto coord = node.GetCoord(idx);
        for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
        {
            res_tmp[idx_eq] = data_manager->GetResidual(idx_eq, idx);
        }
        double jacobi = m_node_metric->GetJacobian(idx);
        auto neighbor = node.GetNeighborNode(idx);
        for (dimension_type dim = 0; dim < grid->GetDim(); dim++)
        {
            riemann_para[0].norm(0) = m_node_metric->GetMetrics(dim, idx)[0];
            riemann_para[0].norm(1) = m_node_metric->GetMetrics(dim, idx)[1];
            riemann_para[0].norm(2) = m_node_metric->GetMetrics(dim, idx)[2];
            riemann_para[0].nt = m_node_metric->GetMetrics(dim, idx)[3];
            riemann_para[1].norm(0) = m_node_metric->GetMetrics(dim, idx)[0];
            riemann_para[1].norm(1) = m_node_metric->GetMetrics(dim, idx)[1];
            riemann_para[1].norm(2) = m_node_metric->GetMetrics(dim, idx)[2];
            riemann_para[1].nt = m_node_metric->GetMetrics(dim, idx)[3];
            MidPointReconstruct2ndOrder(idx, neighbor[2 * dim + 1], riemann_para[0].prim_left.data(),
                                        riemann_para[0].prim_right.data());
            MidPointReconstruct2ndOrder(neighbor[2 * dim], idx, riemann_para[1].prim_left.data(),
                                        riemann_para[1].prim_right.data());
            m_riemann_solver->Solver(riemann_para[0]);
            m_riemann_solver->Solver(riemann_para[1]);
            for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
            {
                res_tmp[idx_eq] -= riemann_para[0].flux[idx_eq] - riemann_para[1].flux[idx_eq];
            }
        }
        data_manager->SetResidual(idx, res_tmp);

        // check negative density and pressure
        // for (int i = 0; i < 6; ++i)
        //{
        //    if (riemann_para[i].prim_left(0) < 0 || riemann_para[i].prim_right(0) < 0 ||
        //        riemann_para[i].prim_left(4) < 0 || riemann_para[i].prim_right(4) < 0)
        //    {
        //        exist_negative = true;
        //        break;
        //    }
        //}
        // if (exist_negative || data_manager->GetNonPhysical(idx) > 0)
        //{
        //    MidPointReconstruct1stOrder(idx, neighbor[1], riemann_para[0].prim_left.data(),
        //                                riemann_para[0].prim_right.data());
        //    MidPointReconstruct1stOrder(neighbor[0], idx, riemann_para[1].prim_left.data(),
        //                                riemann_para[1].prim_right.data());
        //    MidPointReconstruct1stOrder(idx, neighbor[3], riemann_para[2].prim_left.data(),
        //                                riemann_para[2].prim_right.data());
        //    MidPointReconstruct1stOrder(neighbor[2], idx, riemann_para[3].prim_left.data(),
        //                                riemann_para[3].prim_right.data());
        //    MidPointReconstruct1stOrder(idx, neighbor[5], riemann_para[4].prim_left.data(),
        //                                riemann_para[4].prim_right.data());
        //    MidPointReconstruct1stOrder(neighbor[4], idx, riemann_para[5].prim_left.data(),
        //                                riemann_para[5].prim_right.data());
        //}
    }
}

void NSSolverFNFDM::CalcViscousResidual()
{
    return;
    CalcViscousFlux();
    CalcViscousFluxGrad();
    auto grid = GetGrid();
    NodeFN &node = grid->GetNode();
    auto data_manager = GetDataManager();
    int node_num = grid->GetTotalNodeNum();
    int inner_node_num = grid->GetInnerNodeNum();
    auto &inner_node = grid->GetInnerNode();
    auto para = GetPara();
    double res_vis;
    for (int iNode = 0; iNode < inner_node_num; ++iNode)
    {
        int idx = inner_node[iNode];
        for (int iVar = 0; iVar < para->GetEqNum(); ++iVar)
        {
            res_vis = data_manager->GetViscousFluxGrad(iVar, 0, 0, idx) +
                      data_manager->GetViscousFluxGrad(iVar, 1, 1, idx) +
                      data_manager->GetViscousFluxGrad(iVar, 2, 2, idx);
            res_vis /= para->GetDimensionless().GetRe();
            data_manager->SetResidual(iVar, idx, data_manager->GetResidual(iVar, idx) - res_vis);
        }
    }
}
void NSSolverFNFDM::CalcViscousFlux()
{
    auto para = GetPara();
    auto gas = GetGas();
    auto grid = GetGrid();
    auto node = grid->GetNode();
    auto data_manager = GetDataManager();
    int node_num = grid->GetTotalNodeNum();
    double *viscous_flux_x, *viscous_flux_y, *viscous_flux_z;
    double tau_xx, tau_yy, tau_zz, tau_xy, tau_xz, tau_yz;
    double grad_rho[3], grad_u[3], grad_v[3], grad_w[3], grad_p[3];
    double vel[3];
    double tempeture;
    double vis_coef;
    double therm_coef;
    double lamda;
    for (int iNode = 0; iNode < node_num; ++iNode)
    {
        vel[0] = data_manager->GetPrim(1, iNode);
        vel[1] = data_manager->GetPrim(2, iNode);
        vel[2] = data_manager->GetPrim(3, iNode);
        for (int iDim = 0; iDim < 3; ++iDim)
        {
            grad_rho[iDim] = data_manager->GetPrimGrad(0, iDim, iNode);
            grad_u[iDim] = data_manager->GetPrimGrad(1, iDim, iNode);
            grad_v[iDim] = data_manager->GetPrimGrad(2, iDim, iNode);
            grad_w[iDim] = data_manager->GetPrimGrad(3, iDim, iNode);
            grad_p[iDim] = data_manager->GetPrimGrad(4, iDim, iNode);
        }
        tempeture =
            gas->CalcTemperature(data_manager->GetPrim(ID_DENSITY, iNode), data_manager->GetPrim(ID_PRESSURE, iNode));
        vis_coef = gas->CalcMu(tempeture);
        lamda = 2.0 / 3.0 * vis_coef;
        therm_coef = gas->CalcK(tempeture);
        tau_xx = 2 * vis_coef * grad_u[0] - lamda * (grad_u[0] + grad_v[1] + grad_w[2]);
        tau_yy = 2 * vis_coef * grad_v[1] - lamda * (grad_u[0] + grad_v[1] + grad_w[2]);
        tau_zz = 2 * vis_coef * grad_w[2] - lamda * (grad_u[0] + grad_v[1] + grad_w[2]);
        tau_xy = vis_coef * (grad_u[1] + grad_v[0]);
        tau_xz = vis_coef * (grad_u[2] + grad_w[0]);
        tau_yz = vis_coef * (grad_v[2] + grad_w[1]);
        data_manager->SetViscousFlux(0, 0, iNode, 0);
        data_manager->SetViscousFlux(1, 0, iNode, tau_xx);
        data_manager->SetViscousFlux(2, 0, iNode, tau_xy);
        data_manager->SetViscousFlux(3, 0, iNode, tau_xz);
        data_manager->SetViscousFlux(4, 0, iNode,
                                     vel[0] * tau_xx + vel[1] * tau_xy + vel[2] * tau_xz + therm_coef * grad_rho[0]);
        data_manager->SetViscousFlux(0, 1, iNode, 0);
        data_manager->SetViscousFlux(1, 1, iNode, tau_xy);
        data_manager->SetViscousFlux(2, 1, iNode, tau_yy);
        data_manager->SetViscousFlux(3, 1, iNode, tau_yz);
        data_manager->SetViscousFlux(4, 1, iNode,
                                     vel[0] * tau_xy + vel[1] * tau_yy + vel[2] * tau_yz + therm_coef * grad_rho[1]);
        data_manager->SetViscousFlux(0, 2, iNode, 0);
        data_manager->SetViscousFlux(1, 2, iNode, tau_xz);
        data_manager->SetViscousFlux(2, 2, iNode, tau_yz);
        data_manager->SetViscousFlux(3, 2, iNode, tau_zz);
        data_manager->SetViscousFlux(4, 2, iNode,
                                     vel[0] * tau_xz + vel[1] * tau_yz + vel[2] * tau_zz + therm_coef * grad_rho[2]);
    }
}
void NSSolverFNFDM::CalcViscousFluxGrad()
{
    return;
    int equ_num = GetPara()->GetEqNum();
    auto grid = GetGrid();
    auto data_manager = GetDataManager();
    for (int iEqu = 0; iEqu < equ_num; ++iEqu)
    {
        m_grad_wlsq->CalcGradient(
            grid, data_manager->GetViscousFlux(iEqu, 0), data_manager->GetViscousFluxGrad(iEqu, 0, 0),
            data_manager->GetViscousFluxGrad(iEqu, 0, 1), data_manager->GetViscousFluxGrad(iEqu, 0, 2));
        m_grad_wlsq->CalcGradient(
            grid, data_manager->GetViscousFlux(iEqu, 1), data_manager->GetViscousFluxGrad(iEqu, 1, 0),
            data_manager->GetViscousFluxGrad(iEqu, 1, 1), data_manager->GetViscousFluxGrad(iEqu, 1, 2));
        m_grad_wlsq->CalcGradient(
            grid, data_manager->GetViscousFlux(iEqu, 2), data_manager->GetViscousFluxGrad(iEqu, 2, 0),
            data_manager->GetViscousFluxGrad(iEqu, 2, 1), data_manager->GetViscousFluxGrad(iEqu, 2, 2));
    }
}
void NSSolverFNFDM::CalcSourceResidual()
{
}

void NSSolverFNFDM::CalcForce()
{
    return;
    auto grid = GetGrid();
    FaceFN &face = grid->GetFace();
    auto para = GetPara();
    auto data_manager = GetDataManager();
    auto dimensionless = para->GetDimensionless();
    int face_num = face.GetFaceNum();
    double face_pressure;
    double force[3] = {0};
    for (int iFace = 0; iFace < face_num; ++iFace)
    {
        face_pressure = 0;
        auto face2node = face.GetFace2Node(iFace);
        for (int iNode = 0; iNode < face.GetFaceNodeNum(iFace); ++iNode)
        {
            face_pressure += data_manager->GetPrim(ID_PRESSURE, face2node[iNode]);
        }
        face_pressure /= face.GetFaceNodeNum(iFace);
        face_pressure -= dimensionless.GetPressureDL(0);
        force[0] -= face_pressure * face.GetArea(iFace) * face.GetNormal(iFace)[0];
        force[1] -= face_pressure * face.GetArea(iFace) * face.GetNormal(iFace)[1];
        force[2] -= face_pressure * face.GetArea(iFace) * face.GetNormal(iFace)[2];
    }
    force[0] = dimensionless.GetForceReal(force[0]);
    force[1] = dimensionless.GetForceReal(force[1]);
    force[2] = dimensionless.GetForceReal(force[2]);
    Log::info("Force: {}, {}, {}", force[0], force[1], force[2]);
}

shared_ptr<GridFN> NSSolverFNFDM::GetGrid()
{
    return m_grid;
}
} // namespace zaran