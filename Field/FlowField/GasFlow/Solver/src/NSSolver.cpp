#include "NSSolver.h"
#include<fstream>
namespace zaran {
	void NSSolver::Init()
	{
		SetNumberOfEquations(5);
		InitSolver();
		InitData();
		CalcMetric();
	}
	void NSSolver::InitData()
	{
		GridPtr grid = GetGrid();
		FlowSolverPara* para = GetPara();
		int initType = para->GetInitFieldType();
		if (initType == 0)
		{
			InitFieldFarFieldNoVelocity();
		}
		else if (initType == 1)
		{
			InitFieldFarField();
		}
		else if (initType == 2)
		{
			InitFieldRestart();
		}
		else
		{
			Log::warn("Initialize Failed!");
			Log::warn("Wrong Flow field init parameter:{}", initType);
			exit(0);
		}
		Prim2Cons();
		Log::info("Flow Field Initialize Finished!");
	}

	void NSSolver::InitFieldFarField()
	{
		GridPtr grid = GetGrid();
		FlowSolverPara* para = GetPara();
		double inflow_prim[5];
		inflow_prim[0] = para->GetInflowDensity();
		inflow_prim[1] = para->GetInflowVelocityX();
		inflow_prim[2] = para->GetInflowVelocityY();
		inflow_prim[3] = para->GetInflowVelocityZ();
		inflow_prim[4] = para->GetInflowPressure();
		int n_node = grid->GetTotalNodeNum();
		for (int iVal = 0;iVal < GetNumberOfEquations();++iVal)
		{
			for (int iNode = 0; iNode < n_node; ++iNode)
				m_prim[iVal][iNode] = inflow_prim[iVal];
		}
	}

	void NSSolver::InitFieldFarFieldNoVelocity()
	{
		GridPtr grid = GetGrid();
		FlowSolverPara* para = GetPara();
		double* inflow_prim = new double[5];
		inflow_prim[0] = para->GetInflowDensity();
		inflow_prim[1] = 0.0;
		inflow_prim[2] = 0.0;
		inflow_prim[3] = 0.0;
		inflow_prim[4] = para->GetInflowPressure();
		int n_node = grid->GetTotalNodeNum();
		for (int iVal = 0;iVal < GetNumberOfEquations();++iVal)
		{
			for (int iNode = 0; iNode < n_node; ++iNode)
			{
				m_prim[iVal][iNode] = inflow_prim[iVal];
			}
		}
	}

	void NSSolver::InitFieldRestart()
	{
		GridPtr grid = GetGrid();
		FlowSolverPara* para = GetPara();
		std::string restartFileName = "backup.dat";
		std::ifstream fin(restartFileName);
		if (!fin.is_open())
		{
			Log::warn("Restart file not found!");
			exit(0);
		}
		int n_node = grid->GetTotalNodeNum();
		for (int iNode = 0; iNode < n_node; ++iNode)
		{
			for (int iVal = 0;iVal < GetNumberOfEquations();++iVal)
				fin >> m_prim[iVal][iNode];
		}
		fin.close();
	}

	void NSSolver::InitSolver()
	{
		Log::info("Initialize NS Solver!");
		FlowSolver::InitSolver();
		SetNumberOfEquations(5);
		std::string riemannSolverType = GlobalData::GetString("riemannSolver");
		RiemannSolverFactory riemannSolverFactory;
		riemannSolverFactory.Create(riemannSolver_, riemannSolverType);
		Log::info("NS Solver Initialize Finished!");
	}

	void NSSolver::CreateData()
	{

		auto& grid = *GetGrid();
		int nTotalNodeNum = grid.GetTotalNodeNum();
		FieldDataType type = FieldDataType::real;
		auto& dataPtr = GetFieldData();
		// dataPtr = std::make_shared<FieldData>();
		auto& data = *dataPtr;
		data.AddData("rho", type, nTotalNodeNum);
		data.AddData("u", type, nTotalNodeNum);
		data.AddData("v", type, nTotalNodeNum);
		data.AddData("w", type, nTotalNodeNum);
		data.AddData("p", type, nTotalNodeNum);
		data.AddData("cons0", type, nTotalNodeNum);
		data.AddData("cons1", type, nTotalNodeNum);
		data.AddData("cons2", type, nTotalNodeNum);
		data.AddData("cons3", type, nTotalNodeNum);
		data.AddData("cons4", type, nTotalNodeNum);
		data.AddData("res0", type, nTotalNodeNum);
		data.AddData("res1", type, nTotalNodeNum);
		data.AddData("res2", type, nTotalNodeNum);
		data.AddData("res3", type, nTotalNodeNum);
		data.AddData("res4", type, nTotalNodeNum);
		data.AddData("dt", type, nTotalNodeNum);
		data.AddData("limiterCoef0", type, nTotalNodeNum);
		data.AddData("limiterCoef1", type, nTotalNodeNum);
		data.AddData("limiterCoef2", type, nTotalNodeNum);
		data.AddData("limiterCoef3", type, nTotalNodeNum);
		data.AddData("limiterCoef4", type, nTotalNodeNum);
		data.AddData("rhoGradX", type, nTotalNodeNum);
		data.AddData("rhoGradY", type, nTotalNodeNum);
		data.AddData("rhoGradZ", type, nTotalNodeNum);
		data.AddData("uGradX", type, nTotalNodeNum);
		data.AddData("uGradY", type, nTotalNodeNum);
		data.AddData("uGradZ", type, nTotalNodeNum);
		data.AddData("vGradX", type, nTotalNodeNum);
		data.AddData("vGradY", type, nTotalNodeNum);
		data.AddData("vGradZ", type, nTotalNodeNum);
		data.AddData("wGradX", type, nTotalNodeNum);
		data.AddData("wGradY", type, nTotalNodeNum);
		data.AddData("wGradZ", type, nTotalNodeNum);
		data.AddData("pGradX", type, nTotalNodeNum);
		data.AddData("pGradY", type, nTotalNodeNum);
		data.AddData("pGradZ", type, nTotalNodeNum);
		data.AddData("coordTransXXi", type, nTotalNodeNum);
		data.AddData("coordTransXEta", type, nTotalNodeNum);
		data.AddData("coordTransXZeta", type, nTotalNodeNum);
		data.AddData("coordTransXTau", type, nTotalNodeNum);
		data.AddData("coordTransYXi", type, nTotalNodeNum);
		data.AddData("coordTransYEta", type, nTotalNodeNum);
		data.AddData("coordTransYZeta", type, nTotalNodeNum);
		data.AddData("coordTransYTau", type, nTotalNodeNum);
		data.AddData("coordTransZXi", type, nTotalNodeNum);
		data.AddData("coordTransZEta", type, nTotalNodeNum);
		data.AddData("coordTransZZeta", type, nTotalNodeNum);
		data.AddData("coordTransZTau", type, nTotalNodeNum);
		data.AddData("coordTransTXi", type, nTotalNodeNum);
		data.AddData("coordTransTEta", type, nTotalNodeNum);
		data.AddData("coordTransTZeta", type, nTotalNodeNum);
		data.AddData("coordTransTTau", type, nTotalNodeNum);
		data.AddData("coordTransXiX", type, nTotalNodeNum);
		data.AddData("coordTransXiY", type, nTotalNodeNum);
		data.AddData("coordTransXiZ", type, nTotalNodeNum);
		data.AddData("coordTransXiT", type, nTotalNodeNum);
		data.AddData("coordTransEtaX", type, nTotalNodeNum);
		data.AddData("coordTransEtaY", type, nTotalNodeNum);
		data.AddData("coordTransEtaZ", type, nTotalNodeNum);
		data.AddData("coordTransEtaT", type, nTotalNodeNum);
		data.AddData("coordTransZetaX", type, nTotalNodeNum);
		data.AddData("coordTransZetaY", type, nTotalNodeNum);
		data.AddData("coordTransZetaZ", type, nTotalNodeNum);
		data.AddData("coordTransZetaT", type, nTotalNodeNum);
		data.AddData("coordTransTauX", type, nTotalNodeNum);
		data.AddData("coordTransTauY", type, nTotalNodeNum);
		data.AddData("coordTransTauZ", type, nTotalNodeNum);
		data.AddData("coordTransTauT", type, nTotalNodeNum);
		data.AddData("coordTransJ", type, nTotalNodeNum);
		data.AddData("nonPhysical", FieldDataType::integer, nTotalNodeNum);
	}
	void NSSolver::RegisterFieldData()
	{
		FlowSolverPara* para = GetPara();
		auto& data = *GetFieldData();
		m_prim = new double* [5];
		m_cons = new double* [5];
		m_residual = new double* [5];
		m_limiter = new double* [5];
		m_PrimGradX = new double* [5];
		m_PrimGradY = new double* [5];
		m_PrimGradZ = new double* [5];
		m_metric = new double* [33];
		data.GetData("rho", m_prim[0]);
		data.GetData("u", m_prim[1]);
		data.GetData("v", m_prim[2]);
		data.GetData("w", m_prim[3]);
		data.GetData("p", m_prim[4]);
		data.GetData("cons0", m_cons[0]);
		data.GetData("cons1", m_cons[1]);
		data.GetData("cons2", m_cons[2]);
		data.GetData("cons3", m_cons[3]);
		data.GetData("cons4", m_cons[4]);
		data.GetData("res0", m_residual[0]);
		data.GetData("res1", m_residual[1]);
		data.GetData("res2", m_residual[2]);
		data.GetData("res3", m_residual[3]);
		data.GetData("res4", m_residual[4]);
		data.GetData("limiterCoef0", m_limiter[0]);
		data.GetData("limiterCoef1", m_limiter[1]);
		data.GetData("limiterCoef2", m_limiter[2]);
		data.GetData("limiterCoef3", m_limiter[3]);
		data.GetData("limiterCoef4", m_limiter[4]);
		data.GetData("rhoGradX", m_PrimGradX[0]);
		data.GetData("rhoGradY", m_PrimGradY[0]);
		data.GetData("rhoGradZ", m_PrimGradZ[0]);
		data.GetData("uGradX", m_PrimGradX[1]);
		data.GetData("uGradY", m_PrimGradY[1]);
		data.GetData("uGradZ", m_PrimGradZ[1]);
		data.GetData("vGradX", m_PrimGradX[2]);
		data.GetData("vGradY", m_PrimGradY[2]);
		data.GetData("vGradZ", m_PrimGradZ[2]);
		data.GetData("wGradX", m_PrimGradX[3]);
		data.GetData("wGradY", m_PrimGradY[3]);
		data.GetData("wGradZ", m_PrimGradZ[3]);
		data.GetData("pGradX", m_PrimGradX[4]);
		data.GetData("pGradY", m_PrimGradY[4]);
		data.GetData("pGradZ", m_PrimGradZ[4]);
		data.GetData("coordTransXXi", m_metric[0]);
		data.GetData("coordTransXEta", m_metric[1]);
		data.GetData("coordTransXZeta", m_metric[2]);
		data.GetData("coordTransXTau", m_metric[3]);
		data.GetData("coordTransYXi", m_metric[4]);
		data.GetData("coordTransYEta", m_metric[5]);
		data.GetData("coordTransYZeta", m_metric[6]);
		data.GetData("coordTransYTau", m_metric[7]);
		data.GetData("coordTransZXi", m_metric[8]);
		data.GetData("coordTransZEta", m_metric[9]);
		data.GetData("coordTransZZeta", m_metric[10]);
		data.GetData("coordTransZTau", m_metric[11]);
		data.GetData("coordTransTXi", m_metric[12]);
		data.GetData("coordTransTEta", m_metric[13]);
		data.GetData("coordTransTZeta", m_metric[14]);
		data.GetData("coordTransTTau", m_metric[15]);
		data.GetData("coordTransXiX", m_metric[16]);
		data.GetData("coordTransXiY", m_metric[17]);
		data.GetData("coordTransXiZ", m_metric[18]);
		data.GetData("coordTransXiT", m_metric[19]);
		data.GetData("coordTransEtaX", m_metric[20]);
		data.GetData("coordTransEtaY", m_metric[21]);
		data.GetData("coordTransEtaZ", m_metric[22]);
		data.GetData("coordTransEtaT", m_metric[23]);
		data.GetData("coordTransZetaX", m_metric[24]);
		data.GetData("coordTransZetaY", m_metric[25]);
		data.GetData("coordTransZetaZ", m_metric[26]);
		data.GetData("coordTransZetaT", m_metric[27]);
		data.GetData("coordTransTauX", m_metric[28]);
		data.GetData("coordTransTauY", m_metric[29]);
		data.GetData("coordTransTauZ", m_metric[30]);
		data.GetData("coordTransTauT", m_metric[31]);
		data.GetData("coordTransJ", m_metric[32]);
		data.GetData("dt", m_dt);
		data.GetData("nonPhysical", m_non_physical);
	}

	void NSSolver::Solve()
	{
		int iNode = 9587;
		CalcTimeStep();
		CalcPrimGrad();
		CalcLimiter();
		BoundaryCondition();
		TimeAdvance();
		UpdateField();
		CheckPrimtive();
		FixPrimtive();
		CalcForce();
	}
	double NSSolver::ComputeMaxResidual()
	{
		double maxRes = 0;
		auto& rho_res = m_residual[0];
		int n_node = GetGrid()->GetTotalNodeNum();
#pragma omp parallel
		for (int i = 0; i < n_node; ++i)
		{
			double res = abs(rho_res[i]);
			if (res > maxRes)
				maxRes = res;
		}
		return maxRes;
	}

	void NSSolver::BackupField()
	{
		GridPtr grid = GetGrid();
		FlowSolverPara* para = GetPara();
		std::string backupFileName = "backup.dat";
		std::ofstream fout(backupFileName);
		int n_node = grid->GetTotalNodeNum();
		for (int iNode = 0; iNode < n_node; ++iNode)
		{
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
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
		int isSteady = GlobalData::GetInt("isSteady");
		if (isSteady == 0)
		{
			double delta_t = GlobalData::GetDouble("dt");
			double current_time = GlobalData::GetDouble("currentTime");
			double endTime = GlobalData::GetDouble("endTime");
			if (current_time + delta_t > endTime)
			{
				delta_t = endTime - current_time;
				current_time = endTime;
			}
			else
				current_time += delta_t;
			GlobalData::Update("currentTime", current_time);
			SnycTimeStepWithGlobal(delta_t);
		}
	}

	void NSSolver::SnycTimeStepWithGlobal(double& dt)
	{
		GridPtr grid = GetGrid();
		FlowSolverPara* para = GetPara();
#pragma omp parallel for
		for (int iNode = 0;iNode < grid->GetTotalNodeNum();++iNode)
		{
			m_dt[iNode] = dt;
		}
	}

	void NSSolver::TimeAdvance()
	{
		RungeKutta();
	}

	void NSSolver::MidPointReconstruct(int index_left, int index_right, double* value_rec_left, double* value_rec_right)
	{
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		DVector3D r = nodeCoord[index_right] - nodeCoord[index_left];
		for (int iVal = 0;iVal < GetNumberOfEquations();++iVal)
		{
			value_rec_left[iVal] = m_prim[iVal][index_left] +
				0.5 * m_limiter[iVal][index_left] *
				(r.x() * m_PrimGradX[iVal][index_left] +
					r.y() * m_PrimGradY[iVal][index_left] +
					r.z() * m_PrimGradZ[iVal][index_left]);
			value_rec_right[iVal] = m_prim[iVal][index_right] +
				0.5 * m_limiter[iVal][index_right] *
				(r.x() * m_PrimGradX[iVal][index_right] +
					r.y() * m_PrimGradY[iVal][index_right] +
					r.z() * m_PrimGradZ[iVal][index_right]);
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

	void NSSolver::MidPointReconstructFirstOrder(int index_left, int index_right, double* value_rec_left, double* value_rec_right)
	{
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		DVector3D r = nodeCoord[index_right] - nodeCoord[index_left];
		for (int iVal = 0;iVal < GetNumberOfEquations();++iVal)
		{
			value_rec_left[iVal] = m_prim[iVal][index_left];
			value_rec_right[iVal] = m_prim[iVal][index_right];
		}
	}

	void NSSolver::BoundaryCondition()
	{
		GridPtr& grid = GetGrid();
		BoundaryMapPtr& boundaryMapPtr = grid->GetBoundaryMap();
		auto& boundaryMap = boundaryMapPtr->GetBoundaryMap();
		auto& wallBound = boundaryMap["slipWall"];
#pragma omp parallel for
		for (int iBound = 0; iBound < wallBound.size(); ++iBound)
			WallBC(wallBound[iBound]);
		auto& outletBound = boundaryMap["outlet"];
#pragma omp parallel for
		for (int iBound = 0; iBound < outletBound.size(); ++iBound)
			OutletBC(outletBound[iBound]);
		auto& inletBound = boundaryMap["inlet"];
#pragma omp parallel for
		for (int iBound = 0; iBound < inletBound.size(); ++iBound)
			InletBC(inletBound[iBound]);

	}


	void NSSolver::CalcPrimGrad()
	{
		FlowSolverPara* para = GetPara();
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
		auto& grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		BoundaryMapPtr& boundaryMapPtr = grid->GetBoundaryMap();
		auto& boundaryMap = boundaryMapPtr->GetBoundaryMap();
		for (auto& boundary : boundaryMap)
		{
			auto& boundName = boundary.first;
			auto& bound = boundary.second;
			if (boundName == "hole")
				continue;
#pragma omp parallel for
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				for (int iBound = 0; iBound < bound.size(); ++iBound)
				{
					auto& boundIndex = bound[iBound].GetIndex();
					auto& innerIndex = bound[iBound].GetInnerIndex();
					int boundNeighborNum = nodeTopo->GetNeighborCloud()[boundIndex].size();
					// if (boundNeighborNum <= 6)
					// 	ZaranLog::info("Boundary {} has {} neighbors", boundIndex, boundNeighborNum);
					m_PrimGradX[iVal][boundIndex] = 0;
					m_PrimGradY[iVal][boundIndex] = 0;
					m_PrimGradZ[iVal][boundIndex] = 0;
					// m_PrimGradX[iVal][boundIndex] = m_PrimGradX[iVal][innerIndex];
					// m_PrimGradY[iVal][boundIndex] = m_PrimGradY[iVal][innerIndex];
					// m_PrimGradZ[iVal][boundIndex] = m_PrimGradZ[iVal][innerIndex];
				}
			}
		}
	}

	void NSSolver::RungeKutta()
	{
		GridPtr grid = GetGrid();
		FlowSolverPara* para = GetPara();
		const DArray& rk_coef = para->GetRKCoef();
		int rkStage = rk_coef.size();
		auto& nodeTopo = grid->GetNodeTopo();
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
					m_cons[iVal][iNode] = m_cons[iVal][iNode] - rk_coef[iStage] * m_dt[iNode] * m_residual[iVal][iNode] * m_metric[32][iNode];
				}
			}

		}
	}

	void NSSolver::Prim2Cons()
	{
		double gamma = 1.4;
		GridPtr grid = GetGrid();
		int n_node = grid->GetTotalNodeNum();
#pragma omp parallel for
		for (int iNode = 0; iNode < n_node; ++iNode)
			Prim2Cons(m_prim[0][iNode], m_prim[1][iNode], m_prim[2][iNode], m_prim[3][iNode], m_prim[4][iNode],
				m_cons[0][iNode], m_cons[1][iNode], m_cons[2][iNode], m_cons[3][iNode], m_cons[4][iNode]);
	}

	void NSSolver::Prim2Cons(double& rho, double& u, double& v, double& w, double& p, double& cons0, double& cons1, double& cons2, double& cons3, double& cons4)
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

		GridPtr grid = GetGrid();
		int n_node = grid->GetTotalNodeNum();
#pragma omp parallel for
		for (int iNode = 0; iNode < n_node; ++iNode)
		{
			Cons2Prim(m_cons[0][iNode], m_cons[1][iNode], m_cons[2][iNode], m_cons[3][iNode], m_cons[4][iNode],
				m_prim[0][iNode], m_prim[1][iNode], m_prim[2][iNode], m_prim[3][iNode], m_prim[4][iNode]);
		}
	}

	void NSSolver::Cons2Prim(double& cons0, double& cons1, double& cons2, double& cons3, double& cons4, double& rho, double& u, double& v, double& w, double& p)
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
		InviscidFlux();
		ViscousFlux();
		SourceFlux();
	}

	void NSSolver::CalcGradUFDM()
	{
		Log::warn("TO DO Gradient Function UFDM!");
	}

	void NSSolver::UpdateField()
	{
		Cons2Prim();
	}

	void NSSolver::InletBC(Boundary& bound)
	{
		FlowSolverPara* para = GetPara();
		int boundIndex = bound.GetIndex();
		m_prim[0][boundIndex] = para->GetInflowDensity();
		m_prim[1][boundIndex] = para->GetInflowVelocityX();
		m_prim[2][boundIndex] = para->GetInflowVelocityY();
		m_prim[3][boundIndex] = para->GetInflowVelocityZ();
		m_prim[4][boundIndex] = para->GetInflowPressure();
		Prim2Cons(m_prim[0][boundIndex], m_prim[1][boundIndex], m_prim[2][boundIndex], m_prim[3][boundIndex], m_prim[4][boundIndex],
			m_cons[0][boundIndex], m_cons[1][boundIndex], m_cons[2][boundIndex], m_cons[3][boundIndex], m_cons[4][boundIndex]);
	}

	void NSSolver::OutletBC(Boundary& bound)
	{
		int boundIndex = bound.GetIndex();
		int innerIndex = bound.GetInnerIndex();
		m_prim[0][boundIndex] = m_prim[0][innerIndex];
		m_prim[1][boundIndex] = m_prim[1][innerIndex];
		m_prim[2][boundIndex] = m_prim[2][innerIndex];
		m_prim[3][boundIndex] = m_prim[3][innerIndex];
		m_prim[4][boundIndex] = m_prim[4][innerIndex];
		Prim2Cons(m_prim[0][boundIndex], m_prim[1][boundIndex], m_prim[2][boundIndex], m_prim[3][boundIndex], m_prim[4][boundIndex],
			m_cons[0][boundIndex], m_cons[1][boundIndex], m_cons[2][boundIndex], m_cons[3][boundIndex], m_cons[4][boundIndex]);


	}

	void NSSolver::WallBC(Boundary& bound)
	{
		int& innerIndex = bound.GetInnerIndex();
		int boundIndex = bound.GetIndex();
		auto& boundNorm = bound.GetNorm();
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
		Prim2Cons(m_prim[0][boundIndex], m_prim[1][boundIndex], m_prim[2][boundIndex], m_prim[3][boundIndex], m_prim[4][boundIndex],
			m_cons[0][boundIndex], m_cons[1][boundIndex], m_cons[2][boundIndex], m_cons[3][boundIndex], m_cons[4][boundIndex]);
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
		//ComputeBoundaryLimiterCoef();
	}
	void NSSolver::CalcLimiterVK()
	{
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& nodeNeighbor = nodeTopo->GetNeighborCloud();
		int nTotalNodeNum = grid->GetTotalNodeNum();
		double maxVal, minVal;
		double eps = 1e-6;
		double venkatCoeff = 1.0e-5;
#pragma omp parallel for private(maxVal, minVal, eps)
		for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
		{
			for (int iNode = 0; iNode < nTotalNodeNum; ++iNode)
			{
				if (nodeType[iNode] != NodeType::inner && nodeType[iNode] != NodeType::hole)
					continue;
				auto& neighborNode = nodeNeighbor[iNode];
				maxVal = m_prim[iVal][iNode];
				minVal = m_prim[iVal][iNode];
				for (int iNeighbor = 0; iNeighbor < neighborNode.size(); ++iNeighbor)
				{
					maxVal = Max(maxVal, m_prim[iVal][neighborNode[iNeighbor]]);
					minVal = Min(minVal, m_prim[iVal][neighborNode[iNeighbor]]);
				}
				eps = venkatCoeff * (maxVal - minVal);
				eps = eps * eps + SMALL_NUMBER;
				//eps = venkatCoeff * (maxVal - minVal) + SMALL_NUMBER;
				double gradx = m_PrimGradX[iVal][iNode];
				double grady = m_PrimGradY[iVal][iNode];
				double gradz = m_PrimGradZ[iVal][iNode];
				double deltaMax = maxVal - m_prim[iVal][iNode];
				double deltaMin = minVal - m_prim[iVal][iNode];
				double tempCoef = LARGE_NUMBER;
				m_limiter[iVal][iNode] = LARGE_NUMBER;
				for (int iNeighbor = 0; iNeighbor < neighborNode.size(); ++iNeighbor)
				{
					auto current2Neighbor = nodeCoord[neighborNode[iNeighbor]] - nodeCoord[iNode];
					double delta2 = current2Neighbor(0) * gradx + current2Neighbor(1) * grady + current2Neighbor(2) * gradz;
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
					m_limiter[iVal][iNode] = Min(m_limiter[iVal][iNode], tempCoef);
				}
			}
		}
	}
	void NSSolver::CalcLimiterBJ()
	{
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& nodeNeighbor = nodeTopo->GetNeighborCloud();
		int nTotalNodeNum = grid->GetTotalNodeNum();
		double maxVal, minVal;
#pragma omp parallel for private(maxVal, minVal)
		for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
		{
			for (int iNode = 0; iNode < nTotalNodeNum; ++iNode)
			{
				if (nodeType[iNode] != NodeType::inner && nodeType[iNode] != NodeType::hole)
					continue;
				auto& neighborNode = nodeNeighbor[iNode];
				maxVal = m_prim[iVal][iNode];
				minVal = m_prim[iVal][iNode];
				for (int iNeighbor = 0; iNeighbor < neighborNode.size(); ++iNeighbor)
				{
					maxVal = Max(maxVal, m_prim[iVal][neighborNode[iNeighbor]]);
					minVal = Min(minVal, m_prim[iVal][neighborNode[iNeighbor]]);
				}
				double gradx = m_PrimGradX[iVal][iNode];
				double grady = m_PrimGradY[iVal][iNode];
				double gradz = m_PrimGradZ[iVal][iNode];
				double deltaMax = maxVal - m_prim[iVal][iNode];
				double deltaMin = minVal - m_prim[iVal][iNode];
				double tempCoef = LARGE_NUMBER;
				m_limiter[iVal][iNode] = LARGE_NUMBER;
				for (int iNeighbor = 0; iNeighbor < neighborNode.size(); ++iNeighbor)
				{
					auto current2Neighbor = nodeCoord[neighborNode[iNeighbor]] - nodeCoord[iNode];
					double delta2 = current2Neighbor(0) * gradx + current2Neighbor(1) * grady + current2Neighbor(2) * gradz;
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
					m_limiter[iVal][iNode] = Min(m_limiter[iVal][iNode], tempCoef);
				}
			}
		}
	}
	void NSSolver::CalcLimiterNone()
	{
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		auto& limiterCoef = m_limiter;
		int nTotalNodeNum = grid->GetTotalNodeNum();
		for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
		{
			for (int iNode = 0; iNode < nTotalNodeNum; ++iNode)
			{
				if (nodeType[iNode] != NodeType::inner && nodeType[iNode] != NodeType::hole)
					continue;
				m_limiter[iVal][iNode] = 1.0;
			}
		}
	}
	void NSSolver::CalcLimiterFirstOrder()
	{
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& nodeNeighbor = nodeTopo->GetNeighborCloud();
		int nTotalNodeNum = grid->GetTotalNodeNum();
		double maxVal, minVal;
		for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
		{
			for (int iNode = 0; iNode < nTotalNodeNum; ++iNode)
			{
				if (nodeType[iNode] != NodeType::inner && nodeType[iNode] != NodeType::hole)
					continue;
				m_limiter[iVal][iNode] = 0.0;
			}
		}
	}
	void NSSolver::CalcLimiterBound()
	{
		auto& grid = GetGrid();
		BoundaryMapPtr& boundaryMapPtr = grid->GetBoundaryMap();
		auto& boundaryMap = boundaryMapPtr->GetBoundaryMap();
		for (auto& boundary : boundaryMap)
		{
			auto& boundName = boundary.first;
			auto& bound = boundary.second;
			if (boundName == "hole")
				continue;
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				for (int iBound = 0; iBound < bound.size(); ++iBound)
				{
					auto& boundIndex = bound[iBound].GetIndex();
					auto& innerIndex = bound[iBound].GetInnerIndex();
					m_limiter[iVal][boundIndex] = 0.0;
				}
			}
		}
	}
	void NSSolver::CheckPrimtive()
	{
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& nodeNeighbor = nodeTopo->GetNeighborCloud();
		int nTotalNodeNum = grid->GetTotalNodeNum();
		int equation_num = GetNumberOfEquations();
		DArray ave_prim(equation_num, 0.0);
		int nonphysical_node_num = 0;
#pragma omp parallel for private(ave_prim) reduction(+:nonphysical_node_num)
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
				Log::info("Non-physical Node: {}, neighbor num: {}, prim: {},{},{},{},{}", iNode, nodeNeighbor[iNode].size(), m_prim[0][iNode], m_prim[1][iNode], m_prim[2][iNode], m_prim[3][iNode], m_prim[4][iNode]);
				Log::info("Non-physical Node: {}, coord: {},{},{},", iNode, nodeCoord[iNode].x(), nodeCoord[iNode].y(), nodeCoord[iNode].z());
				exit(0);
			}
		}
		if (nonphysical_node_num > 0)
		{
			Log::warn("Non-physical Node Num: {}", nonphysical_node_num);
			FlowSolverPara* para = GetPara();
			double cfl = para->GetCflNumber();
			cfl = cfl / 5.0;
			para->SetCflNumber(cfl);
			Log::warn("CFL Number is reduced to {}", cfl);
		}
		else
		{
			FlowSolverPara* para = GetPara();
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
		GridPtr grid = GetGrid();
		auto& node_topo = grid->GetNodeTopo();
		auto& node_type = node_topo->GetType();
		auto& node_coord = node_topo->GetCoordinate();
		auto& node_neighbor = node_topo->GetNeighborCloud();
		int total_node_num = grid->GetTotalNodeNum();
		int equation_num = GetNumberOfEquations();
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
			auto& neighborNode = node_neighbor[iNode];
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
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
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
				m_cons[0][iNode], m_cons[1][iNode], m_cons[2][iNode], m_cons[3][iNode], m_cons[4][iNode]);
			m_residual[0][iNode] = 0;
		}

	}
	void NSSolver::NoGradient()
	{
		// do nothing
	}

}