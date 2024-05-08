#include "NSSolver.h"
#include "File.h"
#include <corecrt_math.h>
#include <fstream>
namespace zaran
{
void NSSolver::Init()
{
    InitSolver();
    InitField();
    CalcMetric();
}
void NSSolver::InitField()
{
    Grid *grid = GetGrid();
    FlowSolverPara *para = GetPara();
    const InitFieldType &init_type = para->GetInitFieldType();
    if (init_type == InitFieldType::FarFlowNoVelocity)
    {
        InitFieldFarFieldNoVelocity();
    }
    else if (init_type == InitFieldType::FarFlow)
    {
        InitFieldFarFlow();
    }
    else if (init_type == InitFieldType::Backup)
    {
        InitFieldBackup();
    }
    else
    {
        Log::warn("Initialize Failed!");
        Log::warn("Wrong Flow field init parameter:{}", GlobalData::GetString("initFieldType"));
        exit(0);
    }
    Prim2Cons();
    Log::info("Flow Field Initialize Finished!");
}

void NSSolver::InitFieldFarFlow()
{
    Grid *grid = GetGrid();
    FlowSolverPara *para = GetPara();
    double inflow_prim[5];
    inflow_prim[0] = para->GetInflowDensity();
    inflow_prim[1] = para->GetInflowVelocityX();
    inflow_prim[2] = para->GetInflowVelocityY();
    inflow_prim[3] = para->GetInflowVelocityZ();
    inflow_prim[4] = para->GetInflowPressure();
    int n_node = grid->GetTotalNodeNum();
    for (int iVal = 0; iVal < GetEquNum(); ++iVal)
    {
        for (int iNode = 0; iNode < n_node; ++iNode)
            m_prim[iVal][iNode] = inflow_prim[iVal];
    }
}

void NSSolver::InitFieldFarFieldNoVelocity()
{
    Grid *grid = GetGrid();
    FlowSolverPara *para = GetPara();
    double inflow_prim[5];
    inflow_prim[0] = para->GetInflowDensity();
    inflow_prim[1] = 0.0;
    inflow_prim[2] = 0.0;
    inflow_prim[3] = 0.0;
    inflow_prim[4] = para->GetInflowPressure();
    int n_node = grid->GetTotalNodeNum();
    NodeTopo *node_topo = grid->GetNodeTopo();
    auto &node_coord = node_topo->GetCoordinate();
    for (int iNode = 0; iNode < n_node; ++iNode)
    {
        double x = node_coord[iNode].x();
        double y = node_coord[iNode].y();
        double z = node_coord[iNode].z();
        m_prim[0][iNode] = para->GetInflowDensity();
        m_prim[4][iNode] = para->GetInflowPressure();
        if (x > -5.0 && x < 500.0 && y > -90.0 && y < 90.0 && z > -90.0 && z < 90.0)
        {
            m_prim[1][iNode] = 0.0;
            m_prim[2][iNode] = 0.0;
            m_prim[3][iNode] = 0.0;
        }
        else
        {
            m_prim[1][iNode] = para->GetInflowVelocityX();
            m_prim[2][iNode] = para->GetInflowVelocityY();
            m_prim[3][iNode] = para->GetInflowVelocityZ();
        }
    }
}

void NSSolver::InitFieldBackup()
{
    Grid *grid = GetGrid();
    FlowSolverPara *para = GetPara();
    std::string restartFileName = para->GetBackupFieldFileName();
    std::ifstream fin(restartFileName);
    if (!fin.is_open())
    {
        Log::warn("Backup file not found!");
        exit(0);
    }
    int n_node = grid->GetTotalNodeNum();
    for (int iNode = 0; iNode < n_node; ++iNode)
    {
        for (int iVal = 0; iVal < GetEquNum(); ++iVal)
            fin >> m_prim[iVal][iNode];
    }
    fin.close();
}

void NSSolver::InitSolver()
{
    Log::info("Initialize NS Solver!");
    SetEquNum(5);
    FlowSolver::InitSolver();
    std::string riemannSolverType = GlobalData::GetString("riemannSolver");
    RiemannSolverFactory riemannSolverFactory;
    riemannSolverFactory.Create(riemannSolver_, riemannSolverType);
    Log::info("NS Solver Initialize Finished!");
}

void NSSolver::CreateFieldData()
{
    Grid *grid = GetGrid();
    int nTotalNodeNum = grid->GetTotalNodeNum();
    FieldDataType type = FieldDataType::real;
    FieldData *data = GetFieldData();
    data->AddData("density", type, nTotalNodeNum);
    data->AddData("velocity_x", type, nTotalNodeNum);
    data->AddData("velocity_y", type, nTotalNodeNum);
    data->AddData("velocity_w", type, nTotalNodeNum);
    data->AddData("pressure", type, nTotalNodeNum);
    data->AddData("temperture", type, nTotalNodeNum);
    data->AddData("conservative", type, nTotalNodeNum * GetEquNum());
    data->AddData("dt", type, nTotalNodeNum);
    data->AddData("metrics", type, nTotalNodeNum * 17);
    data->AddData("nonPhysical", FieldDataType::integer, nTotalNodeNum);
    data->AddData("residual", type, nTotalNodeNum * GetEquNum());
    data->AddData("limiter_coef", type, nTotalNodeNum * GetEquNum());
    data->AddData("prim_grad", type, nTotalNodeNum * GetEquNum() * 3);
    data->AddData("viscous_flux", type, nTotalNodeNum * GetEquNum());
    data->AddData("viscous_flux_grad", type, nTotalNodeNum * GetEquNum() * 3);
}
void NSSolver::RegisterFieldData()
{
    FlowSolverPara *para = GetPara();
    FieldData *data = GetFieldData();
    m_prim = new double *[5];
    data->GetData("density", m_prim[0]);
    data->GetData("velocity_x", m_prim[1]);
    data->GetData("velocity_y", m_prim[2]);
    data->GetData("velocity_w", m_prim[3]);
    data->GetData("pressure", m_prim[4]);
    data->GetData("conservative", m_cons);
    data->GetData("residual", m_residual);
    data->GetData("limiter_coef", m_limiter);
    data->GetData("prim_grad", m_prim_grad);
    data->GetData("metrics", m_metric);
    data->GetData("dt", m_dt);
    data->GetData("nonPhysical", m_non_physical);
}

void NSSolver::Preprocess()
{
    CalcTimeStep();
    CalcPrimGrad();
    CalcLimiter();
}

void NSSolver::Postprocess()
{
    UpdateField();
    CheckPrimtive();
    FixPrimtive();
    CalcForce();
}

void NSSolver::Solve()
{

    BoundaryCondition();
    TimeAdvance();
}
double NSSolver::ComputeMaxResidual()
{
    double maxRes = 0;
    int n_node = GetGrid()->GetTotalNodeNum();
#pragma omp parallel
    for (int iNode = 0; iNode < n_node; ++iNode)
    {
        double res = fabs(GetResidual(iNode, 0));
        if (res > maxRes)
            maxRes = res;
    }
    return maxRes;
}

void NSSolver::BackupField(std::string &back_folder)
{
    Grid *grid = GetGrid();
    FlowSolverPara *para = GetPara();
    std::string backupFileName = para->GetBackupFieldFileName();
    backupFileName = back_folder + "/" + backupFileName;
    if (IsFileExist(backupFileName))
    {
        DeleteFile(backupFileName);
    }
    std::ofstream fout(backupFileName);
    int n_node = grid->GetTotalNodeNum();
    for (int iNode = 0; iNode < n_node; ++iNode)
    {
        for (int iVal = 0; iVal < GetEquNum(); ++iVal)
        {
            fout << m_prim[iVal][iNode] << " ";
        }
        fout << std::endl;
    }
    fout.close();
}

void NSSolver::CalcTimeStep()
{
    GlobalData::Update("dt", LARGE_NUMBER);
    CalcTimeStepLocal();
    double dt = GlobalData::GetDouble("dt");
    double current_time = GlobalData::GetDouble("currentTime");
    double endTime = GlobalData::GetDouble("endTime");
    if (current_time + dt > endTime)
    {
        dt = endTime - current_time;
        current_time = endTime;
    }
    else
        current_time += dt;
    GlobalData::Update("currentTime", current_time);
    int isSteady = GlobalData::GetInt("isSteady");
    if (isSteady == 0)
    {
        SnycTimeStepWithGlobal(dt);
    }
}

void NSSolver::SnycTimeStepWithGlobal(double &dt)
{
    Grid *grid = GetGrid();
    FlowSolverPara *para = GetPara();
#pragma omp parallel for
    for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
    {
        m_dt[iNode] = dt;
    }
}

void NSSolver::TimeAdvance()
{
    RungeKutta();
}

void NSSolver::MidPointReconstruct(int index_left, int index_right, double *value_rec_left, double *value_rec_right)
{
    Grid *grid = GetGrid();
    NodeTopo *nodeTopo = grid->GetNodeTopo();
    auto &nodeCoord = nodeTopo->GetCoordinate();
    DVector3D r = nodeCoord[index_right] - nodeCoord[index_left];
    for (int iVal = 0; iVal < GetEquNum(); ++iVal)
    {
        value_rec_left[iVal] = m_prim[iVal][index_left] + 0.5 * GetLimiter(index_left, iVal) *
                                                              (r.x() * GetPrimGrad(index_left, iVal, 0) +
                                                               r.y() * GetPrimGrad(index_left, iVal, 1) +
                                                               r.z() * GetPrimGrad(index_left, iVal, 2));
        value_rec_right[iVal] = m_prim[iVal][index_right] - 0.5 * GetLimiter(index_right, iVal) *
                                                                (r.x() * GetPrimGrad(index_right, iVal, 0) +
                                                                 r.y() * GetPrimGrad(index_right, iVal, 1) +
                                                                 r.z() * GetPrimGrad(index_right, iVal, 2));
    }
#if 0
		auto OutputError = [&](int iNode)
			{
				DVector3D xRight, xLeft, yRight, yLeft, zRight, zLeft;
				auto& tempI = nodeTopo->GetTemplateI();
				auto& tempJ = nodeTopo->GetTemplateJ();
				auto& tempK = nodeTopo->GetTemplateK();
				xLeft = nodeCoord[tempI[iNode][0]];
				xRight = nodeCoord[tempI[iNode][2]];
				yLeft = nodeCoord[tempJ[iNode][0]];
				yRight = nodeCoord[tempJ[iNode][2]];
				zLeft = nodeCoord[tempK[iNode][0]];
				zRight = nodeCoord[tempK[iNode][2]];
				// ZaranLog::warn("Node {}: {},{},{}", iNode, nodeCoord[iNode].x(), nodeCoord[iNode].y(), nodeCoord[iNode].z());
				// ZaranLog::info("xLeft index={}: {},{},{}", tempI[iNode][0], xLeft.x(), xLeft.y(), xLeft.z());
				// ZaranLog::info("xRight index={}: {},{},{}", tempI[iNode][2], xRight.x(), xRight.y(), xRight.z());
				// ZaranLog::info("yLeft index={}: {},{},{}", tempJ[iNode][0], yLeft.x(), yLeft.y(), yLeft.z());
				// ZaranLog::info("yRight index={}: {},{},{}", tempJ[iNode][2], yRight.x(), yRight.y(), yRight.z());
				// ZaranLog::info("zLeft index={}: {},{},{}", tempK[iNode][0], zLeft.x(), zLeft.y(), zLeft.z());
				// ZaranLog::info("zRight index={}: {},{},{}", tempK[iNode][2], zRight.x(), zRight.y(), zRight.z());
			};
		//check negative pressure and density
		if (value_rec_left[0] < 0 || value_rec_left[4] < 0)
		{
			for (int iVal = 0;iVal < GetNumberOfEquations();++iVal)
			{
				value_rec_left[iVal] = m_prim[iVal][index_left];
				value_rec_right[iVal] = m_prim[iVal][index_right];
			}
			OutputError(index_left);
		}
		if (value_rec_right[0] < 0 || value_rec_right[4] < 0)
		{
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				value_rec_left[iVal] = m_prim[iVal][index_left];
				value_rec_right[iVal] = m_prim[iVal][index_right];
			}
			OutputError(index_right);
		}
#endif
}

void NSSolver::MidPointReconstructFirstOrder(int index_left, int index_right, double *value_rec_left,
                                             double *value_rec_right)
{
    Grid *grid = GetGrid();
    NodeTopo *nodeTopo = grid->GetNodeTopo();
    auto &nodeCoord = nodeTopo->GetCoordinate();
    DVector3D r = nodeCoord[index_right] - nodeCoord[index_left];
    for (int iVal = 0; iVal < GetEquNum(); ++iVal)
    {
        value_rec_left[iVal] = m_prim[iVal][index_left];
        value_rec_right[iVal] = m_prim[iVal][index_right];
    }
}

void NSSolver::BoundaryCondition()
{
    Grid *grid = GetGrid();
    BoundaryMap *boundaryMap = grid->GetBoundaryMap();
    auto &wallBound = boundaryMap->GetBoundary("slipWall");
#pragma omp parallel for
    for (int iBound = 0; iBound < wallBound.size(); ++iBound)
        WallBC(wallBound[iBound]);
    auto &outletBound = boundaryMap->GetBoundary("outlet");
#pragma omp parallel for
    for (int iBound = 0; iBound < outletBound.size(); ++iBound)
        RiemannBC(outletBound[iBound]);
    auto &inletBound = boundaryMap->GetBoundary("inlet");
#pragma omp parallel for
    for (int iBound = 0; iBound < inletBound.size(); ++iBound)
        RiemannBC(inletBound[iBound]);
}

void NSSolver::CalcPrimGrad()
{
    FlowSolverPara *para = GetPara();
    if (para->GetGradScheme() == GradScheme::wls)
    {
        CalcGradWLS();
    }
    else if (para->GetGradScheme() == GradScheme::ufdm)
    {
        CalcGradUFDM();
    }
    else if (para->GetGradScheme() == GradScheme::noGrad)
    {
        NoGradient();
    }
    else
    {
        Log::warn("Unsupported Gradiend Scheme!");
    }
    CalcPrimGradBound();
}
void NSSolver::CalcPrimGradBound()
{
    Grid *grid = GetGrid();
    NodeTopo *nodeTopo = grid->GetNodeTopo();
    BoundaryMap *boundaryMapPtr = grid->GetBoundaryMap();
    auto &boundaryMap = boundaryMapPtr->GetBoundaryMap();
    for (auto &boundary : boundaryMap)
    {
        auto &boundName = boundary.first;
        auto &bound = boundary.second;
        if (boundName == "hole")
            continue;
#pragma omp parallel for
        for (int iVal = 0; iVal < GetEquNum(); ++iVal)
        {
            for (int iBound = 0; iBound < bound.size(); ++iBound)
            {
                auto &boundIndex = bound[iBound].GetIndex();
                auto &innerIndex = bound[iBound].GetInnerIndex();
                int boundNeighborNum = nodeTopo->GetNeighborCloud()[boundIndex].size();
                GetPrimGrad(boundIndex, iVal, 0) = GetPrimGrad(boundIndex, iVal, 1) = GetPrimGrad(boundIndex, iVal, 2) =
                    0.0;
            }
        }
    }
}

void NSSolver::RungeKutta()
{
    Grid *grid = GetGrid();
    FlowSolverPara *para = GetPara();
    const DArray &rk_coef = para->GetRKCoef();
    int rkStage = rk_coef.size();
    NodeTopo *nodeTopo = grid->GetNodeTopo();
    int nInnerNode = grid->GetInnerNodeNum();
    int nBoundNode = grid->GetBoundNodeNum();
    for (int iStage = 0; iStage < rkStage; ++iStage)
    {
        CalcResidual();
#pragma omp parallel for
        for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
        {
            for (int iVal = 0; iVal < 5; ++iVal)
            {
                GetCons(iNode, iVal) -=
                    rk_coef[iStage] * m_dt[iNode] * GetResidual(iNode, iVal) * GetMetricJacob(iNode);
            }
        }
    }
}

void NSSolver::Prim2Cons()
{
    double gamma = 1.4;
    Grid *grid = GetGrid();
    int n_node = grid->GetTotalNodeNum();
#pragma omp parallel for
    for (int iNode = 0; iNode < n_node; ++iNode)
        Prim2Cons(m_prim[0][iNode], m_prim[1][iNode], m_prim[2][iNode], m_prim[3][iNode], m_prim[4][iNode],
                  GetCons(iNode, 0), GetCons(iNode, 1), GetCons(iNode, 2), GetCons(iNode, 3), GetCons(iNode, 4));
}

void NSSolver::Prim2Cons(double &rho, double &u, double &v, double &w, double &p, double &cons0, double &cons1,
                         double &cons2, double &cons3, double &cons4)
{
    double gamma = 1.4;
    double v2 = u * u + v * v + w * w;
    cons0 = rho;
    cons1 = rho * u;
    cons2 = rho * v;
    cons3 = rho * w;
    cons4 = p / (gamma - 1) + 0.5 * rho * v2;
}

void NSSolver::Cons2Prim()
{

    Grid *grid = GetGrid();
    int n_node = grid->GetTotalNodeNum();
#pragma omp parallel for
    for (int iNode = 0; iNode < n_node; ++iNode)
    {
        Cons2Prim(GetCons(iNode, 0), GetCons(iNode, 1), GetCons(iNode, 2), GetCons(iNode, 3), GetCons(iNode, 4),
                  m_prim[0][iNode], m_prim[1][iNode], m_prim[2][iNode], m_prim[3][iNode], m_prim[4][iNode]);
    }
}

void NSSolver::Cons2Prim(double &cons0, double &cons1, double &cons2, double &cons3, double &cons4, double &rho,
                         double &u, double &v, double &w, double &p)
{
    double gamma = 1.4;
    rho = cons0;
    u = cons1 / rho;
    v = cons2 / rho;
    w = cons3 / rho;
    double v2 = u * u + v * v + w * w;
    p = (gamma - 1) * (cons4 - 0.5 * rho * v2);
}

void NSSolver::CalcResidual()
{
    ZeroResidual();
    ConvectiveResidual();
    ViscousResidual();
    SourceTermResidual();
}

void NSSolver::CalcGradUFDM()
{
    Log::warn("TO DO Gradient Function UFDM!");
}

void NSSolver::UpdateField()
{
    Cons2Prim();
}

void NSSolver::InletBC(Boundary &bound)
{
    FlowSolverPara *para = GetPara();
    int boundIndex = bound.GetIndex();
    m_prim[0][boundIndex] = para->GetInflowDensity();
    m_prim[1][boundIndex] = para->GetInflowVelocityX();
    m_prim[2][boundIndex] = para->GetInflowVelocityY();
    m_prim[3][boundIndex] = para->GetInflowVelocityZ();
    m_prim[4][boundIndex] = para->GetInflowPressure();
    Prim2Cons(m_prim[0][boundIndex], m_prim[1][boundIndex], m_prim[2][boundIndex], m_prim[3][boundIndex],
              m_prim[4][boundIndex], GetCons(boundIndex, 0), GetCons(boundIndex, 1), GetCons(boundIndex, 2),
              GetCons(boundIndex, 3), GetCons(boundIndex, 4));
}

void NSSolver::OutletBC(Boundary &bound)
{
    int boundIndex = bound.GetIndex();
    int innerIndex = bound.GetInnerIndex();
    m_prim[0][boundIndex] = m_prim[0][innerIndex];
    m_prim[1][boundIndex] = m_prim[1][innerIndex];
    m_prim[2][boundIndex] = m_prim[2][innerIndex];
    m_prim[3][boundIndex] = m_prim[3][innerIndex];
    m_prim[4][boundIndex] = m_prim[4][innerIndex];
    Prim2Cons(m_prim[0][boundIndex], m_prim[1][boundIndex], m_prim[2][boundIndex], m_prim[3][boundIndex],
              m_prim[4][boundIndex], GetCons(boundIndex, 0), GetCons(boundIndex, 1), GetCons(boundIndex, 2),
              GetCons(boundIndex, 3), GetCons(boundIndex, 4));
}

void NSSolver::WallBC(Boundary &bound)
{
    int &innerIndex = bound.GetInnerIndex();
    int boundIndex = bound.GetIndex();
    auto &boundNorm = bound.GetNorm();
    m_prim[0][boundIndex] = m_prim[0][innerIndex];
    m_prim[1][boundIndex] = m_prim[1][innerIndex];
    m_prim[2][boundIndex] = m_prim[2][innerIndex];
    m_prim[3][boundIndex] = m_prim[3][innerIndex];
    m_prim[4][boundIndex] = m_prim[4][innerIndex];
    DVector3D innerVel(m_prim[1][innerIndex], m_prim[2][innerIndex], m_prim[3][innerIndex]);
    DVector3D boundVel = innerVel - (innerVel.dot(boundNorm)) * boundNorm / (boundNorm.norm() * boundNorm.norm());
    m_prim[1][boundIndex] = boundVel(0);
    m_prim[2][boundIndex] = boundVel(1);
    m_prim[3][boundIndex] = boundVel(2);
    Prim2Cons(m_prim[0][boundIndex], m_prim[1][boundIndex], m_prim[2][boundIndex], m_prim[3][boundIndex],
              m_prim[4][boundIndex], GetCons(boundIndex, 0), GetCons(boundIndex, 1), GetCons(boundIndex, 2),
              GetCons(boundIndex, 3), GetCons(boundIndex, 4));
}

void NSSolver::RiemannBC(Boundary &bound)
{
    int id_in = bound.GetInnerIndex();
    int id_bound = bound.GetIndex();
    auto &norm_bound = bound.GetNorm();
    double prim_in[5] = {m_prim[0][id_in], m_prim[1][id_in], m_prim[2][id_in], m_prim[3][id_in], m_prim[4][id_in]};
    double prim_bound[5] = {m_prim[0][id_bound], m_prim[1][id_bound], m_prim[2][id_bound], m_prim[3][id_bound],
                            m_prim[4][id_bound]};
    FlowSolverPara *para = GetPara();
    double prim_far[5] = {para->GetInflowDensity(), para->GetInflowVelocityX(), para->GetInflowVelocityY(),
                          para->GetInflowVelocityZ(), para->GetInflowPressure()};
    double vel_in[3] = {m_prim[1][id_in], m_prim[2][id_in], m_prim[3][id_in]};
    // 速度在边界法向上的投影
    double vn_in = vel_in[0] * norm_bound.x() + vel_in[1] * norm_bound.y() + vel_in[2] * norm_bound.z();
    double vel_far[3] = {para->GetInflowVelocityX(), para->GetInflowVelocityY(), para->GetInflowVelocityZ()};
    double vn_far = vel_far[0] * norm_bound.x() + vel_far[1] * norm_bound.y() + vel_far[2] * norm_bound.z();
    double c_in = sqrt(1.4 * m_prim[4][id_in] / m_prim[0][id_in]);
    double c_far = sqrt(1.4 * para->GetInflowPressure() / para->GetInflowDensity());

    double mach = sqrt(vel_in[0] * vel_in[0] + vel_in[1] * vel_in[1] + vel_in[2] * vel_in[2]) / c_in;
    // 超声速出口
    if (mach >= 1.0)
    {
        if (vn_in >= 0) // 超声速出口
        {
            m_prim[0][id_bound] = m_prim[0][id_in];
            m_prim[1][id_bound] = m_prim[1][id_in];
            m_prim[2][id_bound] = m_prim[2][id_in];
            m_prim[3][id_bound] = m_prim[3][id_in];
            m_prim[4][id_bound] = m_prim[4][id_in];
        }
        else // 超声速入口
        {
            m_prim[0][id_bound] = para->GetInflowDensity();
            m_prim[1][id_bound] = para->GetInflowVelocityX();
            m_prim[2][id_bound] = para->GetInflowVelocityY();
            m_prim[3][id_bound] = para->GetInflowVelocityZ();
            m_prim[4][id_bound] = para->GetInflowPressure();
        }
    }
    else
    {
        double r_p = vn_in + 2.0 * c_in / (1.4 - 1.0);
        double r_m = vn_far - 2.0 * c_far / (1.4 - 1.0);
        double vn_bound = 0.5 * (r_p + r_m);
        double c_bound = 0.25 * (1.4 - 1.0) * (r_p - r_m);
        if (vn_bound <= 0) // 亚声速入口
        {
            double entropy = (prim_far[4] / pow(prim_far[0], 1.4));
            m_prim[0][id_bound] = prim_far[0] * pow(c_bound * c_bound / (entropy * 1.4), 1.0 / 1.4);
            m_prim[1][id_bound] = prim_far[1] + norm_bound.x() * (vn_bound - vn_far);
            m_prim[2][id_bound] = prim_far[2] + norm_bound.y() * (vn_bound - vn_far);
            m_prim[3][id_bound] = prim_far[3] + norm_bound.z() * (vn_bound - vn_far);
            m_prim[4][id_bound] = c_bound * c_bound * m_prim[0][id_bound] / 1.4;
        }
        else // 亚声速出口
        {
            double entropy = prim_in[4] / pow(prim_in[0], 1.4);
            m_prim[0][id_bound] = prim_in[0] * pow(c_bound * c_bound / (entropy * 1.4), 1.0 / 1.4);
            m_prim[1][id_bound] = prim_in[1] + norm_bound.x() * (vn_bound - vn_in);
            m_prim[2][id_bound] = prim_in[2] + norm_bound.y() * (vn_bound - vn_in);
            m_prim[3][id_bound] = prim_in[3] + norm_bound.z() * (vn_bound - vn_in);
            m_prim[4][id_bound] = c_bound * c_bound * m_prim[0][id_bound] / 1.4;
        }
    }
}

void NSSolver::CalcLimiter()
{
    string limiterType = GlobalData::GetString("limiterType");
    if (limiterType != "1st-order")
    {
        int firstOrderSteps = GlobalData::GetInt("firstOrderSteps");
        int currentIter = GlobalData::GetInt("currentIter");
        if (currentIter < firstOrderSteps)
        {
            limiterType = "1st-order";
            Log::info("First {}/{} iteration steps use 1st-order scheme", currentIter, firstOrderSteps);
        }
    }
    if (limiterType == "vk")
        CalcLimiterVK();
    else if (limiterType == "barth")
        CalcLimiterBJ();
    else if (limiterType == "noLimiter")
        CalcLimiterNone();
    else if (limiterType == "1st-order")
        CalcLimiterFirstOrder();
    else
        Log::warn("Unsupported Limiter Type: {}", limiterType);
    // ComputeBoundaryLimiterCoef();
}
void NSSolver::CalcLimiterVK()
{
    Grid *grid = GetGrid();
    NodeTopo *nodeTopo = grid->GetNodeTopo();
    auto &nodeType = nodeTopo->GetType();
    auto &nodeCoord = nodeTopo->GetCoordinate();
    auto &nodeNeighbor = nodeTopo->GetNeighborCloud();
    int nTotalNodeNum = grid->GetTotalNodeNum();
    double maxVal, minVal;
    double eps = 1e-6;
    double venkatCoeff = 1.0e-5;
#pragma omp parallel for private(maxVal, minVal, eps)
    for (int iVal = 0; iVal < GetEquNum(); ++iVal)
    {
        for (int iNode = 0; iNode < nTotalNodeNum; ++iNode)
        {
            if (nodeType[iNode] != NodeType::inner && nodeType[iNode] != NodeType::hole)
                continue;
            auto &neighborNode = nodeNeighbor[iNode];
            maxVal = m_prim[iVal][iNode];
            minVal = m_prim[iVal][iNode];
            for (int iNeighbor = 0; iNeighbor < neighborNode.size(); ++iNeighbor)
            {
                maxVal = Max(maxVal, m_prim[iVal][neighborNode[iNeighbor]]);
                minVal = Min(minVal, m_prim[iVal][neighborNode[iNeighbor]]);
            }
            eps = venkatCoeff * (maxVal - minVal);
            eps = eps * eps + SMALL_NUMBER;
            // eps = venkatCoeff * (maxVal - minVal) + SMALL_NUMBER;
            double deltaMax = maxVal - m_prim[iVal][iNode];
            double deltaMin = minVal - m_prim[iVal][iNode];
            double tempCoef = LARGE_NUMBER;
            GetLimiter(iNode, iVal) = LARGE_NUMBER;
            for (int iNeighbor = 0; iNeighbor < neighborNode.size(); ++iNeighbor)
            {
                auto current2Neighbor = nodeCoord[neighborNode[iNeighbor]] - nodeCoord[iNode];
                double delta2 = current2Neighbor(0) * GetPrimGrad(iNode, iVal, 0) +
                                current2Neighbor(1) * GetPrimGrad(iNode, iVal, 1) +
                                current2Neighbor(2) * GetPrimGrad(iNode, iVal, 2);
                delta2 *= 0.5;
                if (delta2 > 0)
                {
                    tempCoef = LimiterVK(maxVal - m_prim[iVal][iNode], delta2, eps);
                }
                else if (delta2 < 0)
                {
                    tempCoef = LimiterVK(minVal - m_prim[iVal][iNode], delta2, eps);
                }
                else
                {
                    tempCoef = 1.0;
                }
                GetLimiter(iNode, iVal) = Min(GetLimiter(iNode, iVal), tempCoef);
            }
        }
    }
}
void NSSolver::CalcLimiterBJ()
{
    Grid *grid = GetGrid();
    NodeTopo *nodeTopo = grid->GetNodeTopo();
    auto &nodeType = nodeTopo->GetType();
    auto &nodeCoord = nodeTopo->GetCoordinate();
    auto &nodeNeighbor = nodeTopo->GetNeighborCloud();
    int nTotalNodeNum = grid->GetTotalNodeNum();
    double maxVal, minVal;
#pragma omp parallel for private(maxVal, minVal)
    for (int iVal = 0; iVal < GetEquNum(); ++iVal)
    {
        for (int iNode = 0; iNode < nTotalNodeNum; ++iNode)
        {
            if (nodeType[iNode] != NodeType::inner && nodeType[iNode] != NodeType::hole)
                continue;
            auto &neighborNode = nodeNeighbor[iNode];
            maxVal = m_prim[iVal][iNode];
            minVal = m_prim[iVal][iNode];
            for (int iNeighbor = 0; iNeighbor < neighborNode.size(); ++iNeighbor)
            {
                maxVal = Max(maxVal, m_prim[iVal][neighborNode[iNeighbor]]);
                minVal = Min(minVal, m_prim[iVal][neighborNode[iNeighbor]]);
            }
            double deltaMax = maxVal - m_prim[iVal][iNode];
            double deltaMin = minVal - m_prim[iVal][iNode];
            double tempCoef = LARGE_NUMBER;
            GetLimiter(iNode, iVal) = LARGE_NUMBER;
            for (int iNeighbor = 0; iNeighbor < neighborNode.size(); ++iNeighbor)
            {
                auto current2Neighbor = nodeCoord[neighborNode[iNeighbor]] - nodeCoord[iNode];
                double delta2 = current2Neighbor(0) * GetPrimGrad(iNode, iVal, 0) +
                                current2Neighbor(1) * GetPrimGrad(iNode, iVal, 1) +
                                current2Neighbor(2) * GetPrimGrad(iNode, iVal, 2);
                delta2 *= 0.5;
                if (delta2 > 0)
                {
                    tempCoef = LimiterBarth(maxVal - m_prim[iVal][iNode], delta2);
                }
                else if (delta2 < 0)
                {
                    tempCoef = LimiterBarth(minVal - m_prim[iVal][iNode], delta2);
                }
                else
                {
                    tempCoef = 1.0;
                }
                GetLimiter(iNode, iVal) = Min(GetLimiter(iNode, iVal), tempCoef);
            }
        }
    }
}
void NSSolver::CalcLimiterNone()
{
    Grid *grid = GetGrid();
    NodeTopo *nodeTopo = grid->GetNodeTopo();
    auto &nodeType = nodeTopo->GetType();
    auto &limiterCoef = m_limiter;
    int nTotalNodeNum = grid->GetTotalNodeNum();
    for (int iVal = 0; iVal < GetEquNum(); ++iVal)
    {
        for (int iNode = 0; iNode < nTotalNodeNum; ++iNode)
        {
            if (nodeType[iNode] != NodeType::inner && nodeType[iNode] != NodeType::hole)
                continue;
            GetLimiter(iNode, iVal) = 1.0;
        }
    }
}
void NSSolver::CalcLimiterFirstOrder()
{
    Grid *grid = GetGrid();
    NodeTopo *nodeTopo = grid->GetNodeTopo();

    auto &nodeType = nodeTopo->GetType();
    auto &nodeCoord = nodeTopo->GetCoordinate();
    auto &nodeNeighbor = nodeTopo->GetNeighborCloud();
    int nTotalNodeNum = grid->GetTotalNodeNum();
    double maxVal, minVal;
    for (int iVal = 0; iVal < GetEquNum(); ++iVal)
    {
        for (int iNode = 0; iNode < nTotalNodeNum; ++iNode)
        {
            if (nodeType[iNode] != NodeType::inner && nodeType[iNode] != NodeType::hole)
                continue;
            GetLimiter(iNode, iVal) = 0.0;
        }
    }
}
void NSSolver::CalcLimiterBound()
{
    Grid *grid = GetGrid();
    BoundaryMap *boundaryMapPtr = grid->GetBoundaryMap();
    auto &boundaryMap = boundaryMapPtr->GetBoundaryMap();
    for (auto &boundary : boundaryMap)
    {
        auto &boundName = boundary.first;
        auto &bound = boundary.second;
        if (boundName == "hole")
            continue;
        for (int iVal = 0; iVal < GetEquNum(); ++iVal)
        {
            for (int iBound = 0; iBound < bound.size(); ++iBound)
            {
                auto &boundIndex = bound[iBound].GetIndex();
                auto &innerIndex = bound[iBound].GetInnerIndex();
                GetLimiter(boundIndex, iVal) = 0.0;
            }
        }
    }
}
void NSSolver::CheckPrimtive()
{
    Grid *grid = GetGrid();
    NodeTopo *nodeTopo = grid->GetNodeTopo();

    auto &nodeType = nodeTopo->GetType();
    auto &nodeCoord = nodeTopo->GetCoordinate();
    auto &nodeNeighbor = nodeTopo->GetNeighborCloud();
    int nTotalNodeNum = grid->GetTotalNodeNum();
    int equation_num = GetEquNum();
    DArray ave_prim(equation_num, 0.0);
    int nonphysical_node_num = 0;
#pragma omp parallel for private(ave_prim) reduction(+ : nonphysical_node_num)
    for (int iNode = 0; iNode < nTotalNodeNum; ++iNode)
    {
        bool exist_nonphysical = false;
        m_non_physical[iNode] = -1;
        if (nodeType[iNode] != NodeType::inner && nodeType[iNode] != NodeType::hole)
            continue;
        if (m_prim[0][iNode] < 0 || m_prim[4][iNode] < 0)
        {
            exist_nonphysical = true;
        }
        if (!exist_nonphysical)
        {
            for (int iVal = 0; iVal < equation_num; ++iVal)
            {
                if (isnan(m_prim[iVal][iNode]) || isinf(m_prim[iVal][iNode]))
                {
                    exist_nonphysical = true;
                    break;
                }
            }
        }
        if (exist_nonphysical)
        {
            m_non_physical[iNode] = 1;
            nonphysical_node_num++;
            Log::info("Non-physical Node: {}, neighbor num: {}, prim: {},{},{},{},{}", iNode,
                      nodeNeighbor[iNode].size(), m_prim[0][iNode], m_prim[1][iNode], m_prim[2][iNode],
                      m_prim[3][iNode], m_prim[4][iNode]);
            Log::info("Non-physical Node: {}, coord: {},{},{},", iNode, nodeCoord[iNode].x(), nodeCoord[iNode].y(),
                      nodeCoord[iNode].z());
        }
    }
    if (nonphysical_node_num > 0)
    {
        Log::warn("Non-physical Node Num: {}", nonphysical_node_num);
        FlowSolverPara *para = GetPara();
        double cfl = para->GetCflNumber();
        cfl = cfl / 5.0;
        para->SetCflNumber(cfl);
        Log::warn("CFL Number is reduced to {}", cfl);
    }
    else
    {
        FlowSolverPara *para = GetPara();
        double cfl = para->GetCflNumber();
        double cfl_max = GlobalData::GetDouble("cflNumber");
        if (cfl < cfl_max)
        {
            cfl = Min(cfl_max, cfl * 1.6);
            para->SetCflNumber(cfl);
            Log::info("CFL Number is increased to {}", cfl);
        }
        else
            para->SetCflNumber(cfl_max);
    }
}

void NSSolver::CheckResidual()
{
}

void NSSolver::FixPrimtive()
{
    Grid *grid = GetGrid();
    NodeTopo *node_topo = grid->GetNodeTopo();
    auto &node_type = node_topo->GetType();
    auto &node_coord = node_topo->GetCoordinate();
    auto &node_neighbor = node_topo->GetNeighborCloud();
    int total_node_num = grid->GetTotalNodeNum();
    int equation_num = GetEquNum();
    DArray weight, distance;
    IArray physical_neighbor;
    double sum = 0;
#pragma omp parallel for private(physical_neighbor, weight, distance, sum)
    for (int iNode = 0; iNode < total_node_num; ++iNode)
    {
        if (node_type[iNode] != NodeType::inner && node_type[iNode] != NodeType::hole)
            continue;
        if (m_non_physical[iNode] == -1)
            continue;
        physical_neighbor.clear();
        auto &neighborNode = node_neighbor[iNode];
        for (int iNeighbor = 0; iNeighbor < neighborNode.size(); ++iNeighbor)
        {
            if (m_non_physical[neighborNode[iNeighbor]] == 1)
                continue;
            physical_neighbor.push_back(neighborNode[iNeighbor]);
        }
        weight.resize(physical_neighbor.size());
        distance.resize(physical_neighbor.size());
        sum = 0;
        for (int iNeighbor = 0; iNeighbor < physical_neighbor.size(); ++iNeighbor)
        {
            if (m_non_physical[physical_neighbor[iNeighbor]] == 1)
                continue;
            distance[iNeighbor] = (node_coord[physical_neighbor[iNeighbor]] - node_coord[iNode]).norm();
            sum += 1.0 / distance[iNeighbor];
        }
        for (int iNeighbor = 0; iNeighbor < physical_neighbor.size(); ++iNeighbor)
        {
            if (m_non_physical[physical_neighbor[iNeighbor]] == 1)
                continue;
            weight[iNeighbor] = 1.0 / (distance[iNeighbor] * sum);
        }
        for (int iVal = 0; iVal < GetEquNum(); ++iVal)
        {
            if (isnan(m_prim[iVal][iNode]) || isinf(m_prim[iVal][iNode]))
            {
                m_prim[iVal][iNode] = 0;
                for (int iNeighbor = 0; iNeighbor < physical_neighbor.size(); ++iNeighbor)
                {
                    if (m_non_physical[physical_neighbor[iNeighbor]] == 1)
                        continue;
                    m_prim[iVal][iNode] += m_prim[iVal][physical_neighbor[iNeighbor]] * weight[iNeighbor];
                }
            }
        }
        if (m_prim[0][iNode] < 0)
            m_prim[0][iNode] = 0;
        for (int iNeighbor = 0; iNeighbor < physical_neighbor.size(); ++iNeighbor)
        {
            if (m_non_physical[physical_neighbor[iNeighbor]] == 1)
                continue;
            m_prim[0][iNode] += m_prim[0][physical_neighbor[iNeighbor]] * weight[iNeighbor];
        }
        if (m_prim[4][iNode] < 0)
            m_prim[4][iNode] = 0;
        for (int iNeighbor = 0; iNeighbor < physical_neighbor.size(); ++iNeighbor)
        {
            if (m_non_physical[physical_neighbor[iNeighbor]] == 1)
                continue;
            m_prim[4][iNode] += m_prim[4][physical_neighbor[iNeighbor]] * weight[iNeighbor];
        }
        Prim2Cons(m_prim[0][iNode], m_prim[1][iNode], m_prim[2][iNode], m_prim[3][iNode], m_prim[4][iNode],
                  GetCons(iNode, 0), GetCons(iNode, 1), GetCons(iNode, 2), GetCons(iNode, 3), GetCons(iNode, 4));
        GetResidual(iNode, 0) = 0;
    }
}
void NSSolver::NoGradient()
{
    // do nothing
}

} // namespace zaran