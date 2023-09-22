#include "NSSolver.h"
#include<fstream>
namespace zaran {
	void NSSolver::Init()
	{
		InitSolver();
		InitField();
		ComputeCoordTrans();
	}
	void NSSolver::InitField()
	{
		GridPtr grid = GetGrid();
		auto& rho = *m_Primtive[0];
		auto& u = *m_Primtive[1];
		auto& v = *m_Primtive[2];
		auto& w = *m_Primtive[3];
		auto& p = *m_Primtive[4];
		FlowSolverParaPtr para = GetPara();
		int initType = para->GetInitFieldType();
		DVector primInit = para->GetPrimitiveInflow();
		ZaranLog::info("primitive inlet:({},{},{},{},{})", primInit(0), primInit(1), primInit(2), primInit(3), primInit(4));
		if (initType == 0)
		{
			primInit(1) = primInit(2) = primInit(3) = 0;
		}
		else if (initType == 1)
		{
			//do nothing
		}
		else
		{
			ZaranLog::warn("流场初始化失败!");
			ZaranLog::warn("初始化流场参数设置错误:{}", initType);
			exit(0);
		}
		int nTotalNodeNum = grid->GetTotalNodeNum();
		double x, y, z;
		for (int iNode = 0; iNode < nTotalNodeNum; ++iNode)
		{
			rho[iNode] = primInit(0);
			u[iNode] = primInit(1);
			v[iNode] = primInit(2);
			w[iNode] = primInit(3);
			p[iNode] = primInit(4);
		}
		Primitive2Conservative();
	}

	void NSSolver::InitSolver()
	{
		FlowSolver::InitSolver();
		SetNumberOfEquations(5);
		std::string riemannSolverType = GlobalData::GetString("riemannSolver");
		RiemannSolverFactory riemannSolverFactory;
		riemannSolverFactory.Create(riemannSolver_, riemannSolverType);
		CreateLimiter();
	}

	void NSSolver::CreateFieldData()
	{

		auto& grid = *GetGrid();
		int nTotalNodeNum = grid.GetTotalNodeNum();
		DArray emptyData(nTotalNodeNum);
		auto& dataPtr = GetFieldData();
		dataPtr = std::make_shared<FieldData>();
		auto& data = *dataPtr;
		data.AddData("rho", emptyData);
		data.AddData("u", emptyData);
		data.AddData("v", emptyData);
		data.AddData("w", emptyData);
		data.AddData("p", emptyData);
		data.AddData("cons0", emptyData);
		data.AddData("cons1", emptyData);
		data.AddData("cons2", emptyData);
		data.AddData("cons3", emptyData);
		data.AddData("cons4", emptyData);
		data.AddData("res0", emptyData);
		data.AddData("res1", emptyData);
		data.AddData("res2", emptyData);
		data.AddData("res3", emptyData);
		data.AddData("res4", emptyData);
		data.AddData("dt", emptyData);
		data.AddData("limiterCoef0", emptyData);
		data.AddData("limiterCoef1", emptyData);
		data.AddData("limiterCoef2", emptyData);
		data.AddData("limiterCoef3", emptyData);
		data.AddData("limiterCoef4", emptyData);
		data.AddData("rhoGradX", emptyData);
		data.AddData("rhoGradY", emptyData);
		data.AddData("rhoGradZ", emptyData);
		data.AddData("uGradX", emptyData);
		data.AddData("uGradY", emptyData);
		data.AddData("uGradZ", emptyData);
		data.AddData("vGradX", emptyData);
		data.AddData("vGradY", emptyData);
		data.AddData("vGradZ", emptyData);
		data.AddData("wGradX", emptyData);
		data.AddData("wGradY", emptyData);
		data.AddData("wGradZ", emptyData);
		data.AddData("pGradX", emptyData);
		data.AddData("pGradY", emptyData);
		data.AddData("pGradZ", emptyData);
		data.AddData("consRK0", emptyData);
		data.AddData("consRK1", emptyData);
		data.AddData("consRK2", emptyData);
		data.AddData("consRK3", emptyData);
		data.AddData("consRK4", emptyData);
		data.AddData("coordTransXXi", emptyData);
		data.AddData("coordTransXEta", emptyData);
		data.AddData("coordTransXZeta", emptyData);
		data.AddData("coordTransXTau", emptyData);
		data.AddData("coordTransYXi", emptyData);
		data.AddData("coordTransYEta", emptyData);
		data.AddData("coordTransYZeta", emptyData);
		data.AddData("coordTransYTau", emptyData);
		data.AddData("coordTransZXi", emptyData);
		data.AddData("coordTransZEta", emptyData);
		data.AddData("coordTransZZeta", emptyData);
		data.AddData("coordTransZTau", emptyData);
		data.AddData("coordTransTXi", emptyData);
		data.AddData("coordTransTEta", emptyData);
		data.AddData("coordTransTZeta", emptyData);
		data.AddData("coordTransTTau", emptyData);
		data.AddData("coordTransXiX", emptyData);
		data.AddData("coordTransXiY", emptyData);
		data.AddData("coordTransXiZ", emptyData);
		data.AddData("coordTransXiT", emptyData);
		data.AddData("coordTransEtaX", emptyData);
		data.AddData("coordTransEtaY", emptyData);
		data.AddData("coordTransEtaZ", emptyData);
		data.AddData("coordTransEtaT", emptyData);
		data.AddData("coordTransZetaX", emptyData);
		data.AddData("coordTransZetaY", emptyData);
		data.AddData("coordTransZetaZ", emptyData);
		data.AddData("coordTransZetaT", emptyData);
		data.AddData("coordTransTauX", emptyData);
		data.AddData("coordTransTauY", emptyData);
		data.AddData("coordTransTauZ", emptyData);
		data.AddData("coordTransTauT", emptyData);
		data.AddData("coordTransJ", emptyData);
	}
	void NSSolver::RegisterFieldData()
	{
		auto& data = *GetFieldData();
		m_Primtive.push_back(&data.GetData("rho"));
		m_Primtive.push_back(&data.GetData("u"));
		m_Primtive.push_back(&data.GetData("v"));
		m_Primtive.push_back(&data.GetData("w"));
		m_Primtive.push_back(&data.GetData("p"));
		m_Conservative.push_back(&data.GetData("cons0"));
		m_Conservative.push_back(&data.GetData("cons1"));
		m_Conservative.push_back(&data.GetData("cons2"));
		m_Conservative.push_back(&data.GetData("cons3"));
		m_Conservative.push_back(&data.GetData("cons4"));
		m_Residual.push_back(&data.GetData("res0"));
		m_Residual.push_back(&data.GetData("res1"));
		m_Residual.push_back(&data.GetData("res2"));
		m_Residual.push_back(&data.GetData("res3"));
		m_Residual.push_back(&data.GetData("res4"));
		m_TimeStep = &data.GetData("dt");
		m_LimiterCoef.push_back(&data.GetData("limiterCoef0"));
		m_LimiterCoef.push_back(&data.GetData("limiterCoef1"));
		m_LimiterCoef.push_back(&data.GetData("limiterCoef2"));
		m_LimiterCoef.push_back(&data.GetData("limiterCoef3"));
		m_LimiterCoef.push_back(&data.GetData("limiterCoef4"));
		m_PrimGradX.push_back(&data.GetData("rhoGradX"));
		m_PrimGradX.push_back(&data.GetData("uGradX"));
		m_PrimGradX.push_back(&data.GetData("vGradX"));
		m_PrimGradX.push_back(&data.GetData("wGradX"));
		m_PrimGradX.push_back(&data.GetData("pGradX"));
		m_PrimGradY.push_back(&data.GetData("rhoGradY"));
		m_PrimGradY.push_back(&data.GetData("uGradY"));
		m_PrimGradY.push_back(&data.GetData("vGradY"));
		m_PrimGradY.push_back(&data.GetData("wGradY"));
		m_PrimGradY.push_back(&data.GetData("pGradY"));
		m_PrimGradZ.push_back(&data.GetData("rhoGradZ"));
		m_PrimGradZ.push_back(&data.GetData("uGradZ"));
		m_PrimGradZ.push_back(&data.GetData("vGradZ"));
		m_PrimGradZ.push_back(&data.GetData("wGradZ"));
		m_PrimGradZ.push_back(&data.GetData("pGradZ"));
		m_ConservativeRK.push_back(&data.GetData("consRK0"));
		m_ConservativeRK.push_back(&data.GetData("consRK1"));
		m_ConservativeRK.push_back(&data.GetData("consRK2"));
		m_ConservativeRK.push_back(&data.GetData("consRK3"));
		m_ConservativeRK.push_back(&data.GetData("consRK4"));
		m_CoordTrans.push_back(&data.GetData("coordTransXXi"));
		m_CoordTrans.push_back(&data.GetData("coordTransXEta"));
		m_CoordTrans.push_back(&data.GetData("coordTransXZeta"));
		m_CoordTrans.push_back(&data.GetData("coordTransXTau"));
		m_CoordTrans.push_back(&data.GetData("coordTransYXi"));
		m_CoordTrans.push_back(&data.GetData("coordTransYEta"));
		m_CoordTrans.push_back(&data.GetData("coordTransYZeta"));
		m_CoordTrans.push_back(&data.GetData("coordTransYTau"));
		m_CoordTrans.push_back(&data.GetData("coordTransZXi"));
		m_CoordTrans.push_back(&data.GetData("coordTransZEta"));
		m_CoordTrans.push_back(&data.GetData("coordTransZZeta"));
		m_CoordTrans.push_back(&data.GetData("coordTransZTau"));
		m_CoordTrans.push_back(&data.GetData("coordTransTXi"));
		m_CoordTrans.push_back(&data.GetData("coordTransTEta"));
		m_CoordTrans.push_back(&data.GetData("coordTransTZeta"));
		m_CoordTrans.push_back(&data.GetData("coordTransTTau"));
		m_CoordTrans.push_back(&data.GetData("coordTransXiX"));
		m_CoordTrans.push_back(&data.GetData("coordTransXiY"));
		m_CoordTrans.push_back(&data.GetData("coordTransXiZ"));
		m_CoordTrans.push_back(&data.GetData("coordTransXiT"));
		m_CoordTrans.push_back(&data.GetData("coordTransEtaX"));
		m_CoordTrans.push_back(&data.GetData("coordTransEtaY"));
		m_CoordTrans.push_back(&data.GetData("coordTransEtaZ"));
		m_CoordTrans.push_back(&data.GetData("coordTransEtaT"));
		m_CoordTrans.push_back(&data.GetData("coordTransZetaX"));
		m_CoordTrans.push_back(&data.GetData("coordTransZetaY"));
		m_CoordTrans.push_back(&data.GetData("coordTransZetaZ"));
		m_CoordTrans.push_back(&data.GetData("coordTransZetaT"));
		m_CoordTrans.push_back(&data.GetData("coordTransTauX"));
		m_CoordTrans.push_back(&data.GetData("coordTransTauY"));
		m_CoordTrans.push_back(&data.GetData("coordTransTauZ"));
		m_CoordTrans.push_back(&data.GetData("coordTransTauT"));
		m_CoordTrans.push_back(&data.GetData("coordTransJ"));
	}
	void NSSolver::Solve()
	{
		ComputeTimeStep();
		ComputePrimtiveGradient();
		ComputeLimiterCoef();
		TimeAdvance();
		UpdateField();
		BoundaryCondition();
	}
	double NSSolver::ComputeMaxResidual()
	{
		double maxRes = 0;
		auto& rhoRes = *m_Residual[0];
		for (int i = 0; i < rhoRes.size(); ++i)
		{
			double res = abs(rhoRes[i]);
			if (res > maxRes)
				maxRes = res;
		}
		return maxRes;
	}

	void NSSolver::ComputeTimeStep()
	{
		GlobalData::Update("dt", LARGE_NUMBER);
		ComputeTimeStepLocal();
		//是否使用全局时间步
		int useGlobalTimeStep = GlobalData::GetInt("useGlobalTimeStep");
		if (useGlobalTimeStep == 1)
		{
			double globlaTimeStep = GlobalData::GetDouble("dt");
			double globalTime = GlobalData::GetDouble("globalTime");
			globalTime += globlaTimeStep;
			GlobalData::Update("globalTime", globalTime);
			SnycTimeStepWithGlobal(globlaTimeStep);
		}
	}

	void NSSolver::SnycTimeStepWithGlobal(double& dt)
	{
		GridPtr grid = GetGrid();
		FlowSolverParaPtr para = GetPara();
		for (int iNode = 0; iNode < m_TimeStep->size(); ++iNode)
		{
			(*m_TimeStep)[iNode] = dt;
		}
	}

	void NSSolver::TimeAdvance()
	{
		RungeKutta();
	}

	void NSSolver::BoundaryCondition()
	{
		GridPtr& grid = GetGrid();
		BoundaryMapPtr& boundaryMapPtr = grid->GetBoundaryMap();
		auto& boundaryMap = boundaryMapPtr->GetBoundaryMap();
		auto& wallBound = boundaryMap["slipWall"];
		for (int iBound = 0; iBound < wallBound.size(); ++iBound)
			ComputeWallBC(wallBound[iBound]);
		auto& outletBound = boundaryMap["outlet"];
		for (int iBound = 0; iBound < outletBound.size(); ++iBound)
			ComputeOutletBC(outletBound[iBound]);
		auto& inletletBound = boundaryMap["inlet"];
		for (int iBound = 0; iBound < inletletBound.size(); ++iBound)
			ComputeInletBC(inletletBound[iBound]);
	}


	void NSSolver::ComputePrimtiveGradient()
	{
		FlowSolverParaPtr para = GetPara();
		if (para->GetGradScheme() == GradScheme::wls)
		{
			ComputeGradientWLS();
		}
		else if (para->GetGradScheme() == GradScheme::ufdm)
		{
			ComputeGradientUFDM();
		}
		else if (para->GetGradScheme() == GradScheme::noGrad)
		{
			NoGradient();
		}
		else
		{
			ZaranLog::warn("Unsupported Gradiend Scheme!");
		}
		ComputeBoundaryPrimtiveGradient();
	}
	void NSSolver::ComputeBoundaryPrimtiveGradient()
	{
		auto& grid = GetGrid();
		BoundaryMapPtr& boundaryMapPtr = grid->GetBoundaryMap();
		auto& boundaryMap = boundaryMapPtr->GetBoundaryMap();
		auto& primGradX = m_PrimGradX;
		auto& primGradY = m_PrimGradY;
		auto& primGradZ = m_PrimGradZ;
		for (auto& boundary : boundaryMap)
		{
			auto& boundName = boundary.first;
			auto& bound = boundary.second;
			for (int iBound = 0; iBound < bound.size(); ++iBound)
			{
				auto& boundIndex = bound[iBound].GetIndex();
				auto& innerIndex = bound[iBound].GetInnerNodeIndex();
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					(*primGradX[iVal])[boundIndex] = 0;
					(*primGradY[iVal])[boundIndex] = 0;
					(*primGradZ[iVal])[boundIndex] = 0;
					//(*primGradX[iVal])[boundIndex] = (*primGradX[iVal])[innerIndex];
					//(*primGradY[iVal])[boundIndex] = (*primGradY[iVal])[innerIndex];
					//(*primGradZ[iVal])[boundIndex] = (*primGradZ[iVal])[innerIndex];
				}
			}
		}
	}

	void NSSolver::RungeKutta()
	{
		GridPtr grid = GetGrid();
		FlowSolverParaPtr para = GetPara();
		const DArray& rkCoef = para->GetRKCoef();
		int rkStage = rkCoef.size();
		auto& cons = m_Conservative;
		auto& cons_RK = m_ConservativeRK;
		auto& dt = *m_TimeStep;
		auto& res = m_Residual;
		auto& coordTrans = m_CoordTrans;
		auto& nodeTopo = grid->GetNodeTopo();
		int nInnerNode = grid->GetInnerNodeNum();
		int nBoundNode = grid->GetBoundNodeNum();
		for (int iStage = 0; iStage < rkStage; ++iStage)
		{
			ComputeResidual();
			for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
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

	void NSSolver::Primitive2Conservative()
	{
		double gamma = 1.4;
		GridPtr grid = GetGrid();
		auto& rho = *m_Primtive[0];
		auto& u = *m_Primtive[1];
		auto& v = *m_Primtive[2];
		auto& w = *m_Primtive[3];
		auto& p = *m_Primtive[4];
		auto& cons0 = *m_Conservative[0];
		auto& cons1 = *m_Conservative[1];
		auto& cons2 = *m_Conservative[2];
		auto& cons3 = *m_Conservative[3];
		auto& cons4 = *m_Conservative[4];
		int nTotalNode = grid->GetTotalNodeNum();
		for (int iNode = 0; iNode < nTotalNode; ++iNode)
			Primitive2Conservative(rho[iNode], u[iNode], v[iNode], w[iNode], p[iNode], cons0[iNode], cons1[iNode], cons2[iNode], cons3[iNode], cons4[iNode]);
	}

	void NSSolver::Primitive2Conservative(double& rho, double& u, double& v, double& w, double& p, double& cons0, double& cons1, double& cons2, double& cons3, double& cons4)
	{
		double gamma = 1.4;
		double v2 = u * u + v * v + w * w;
		cons0 = rho;
		cons1 = rho * u;
		cons2 = rho * v;
		cons3 = rho * w;
		cons4 = p / (gamma - 1) + 0.5 * rho * v2;
	}

	void NSSolver::Conservative2Primitive()
	{

		GridPtr grid = GetGrid();
		auto& rho = *m_Primtive[0];
		auto& u = *m_Primtive[1];
		auto& v = *m_Primtive[2];
		auto& w = *m_Primtive[3];
		auto& p = *m_Primtive[4];
		auto& cons0 = *m_Conservative[0];
		auto& cons1 = *m_Conservative[1];
		auto& cons2 = *m_Conservative[2];
		auto& cons3 = *m_Conservative[3];
		auto& cons4 = *m_Conservative[4];
		auto& inflowPrim = GetPara()->GetPrimitiveInflow();
		int nTotalNode = grid->GetTotalNodeNum();
		for (int iNode = 0; iNode < nTotalNode; ++iNode)
		{
			Conservative2Primitive(cons0[iNode], cons1[iNode], cons2[iNode], cons3[iNode], cons4[iNode], rho[iNode], u[iNode], v[iNode], w[iNode], p[iNode]);
		}
	}

	void NSSolver::Conservative2Primitive(double& cons0, double& cons1, double& cons2, double& cons3, double& cons4, double& rho, double& u, double& v, double& w, double& p)
	{
		double gamma = 1.4;
		rho = cons0;
		u = cons1 / rho;
		v = cons2 / rho;
		w = cons3 / rho;
		double v2 = u * u + v * v + w * w;
		p = (gamma - 1) * (cons4 - 0.5 * rho * v2);
	}

	void NSSolver::ComputeResidual()
	{
		ZeroResidual();
		InviscidFlux();
		ViscousFlux();
		SourceFlux();
	}

	void NSSolver::ComputeGradientUFDM()
	{
		ZaranLog::warn("TO DO Gradient Function UFDM!");
	}

	void NSSolver::UpdateField()
	{
		Conservative2Primitive();
	}

	void NSSolver::ComputeInletBC(Boundary& bound)
	{
		FlowSolverParaPtr para = GetPara();
		int boundIndex = bound.GetIndex();
		auto& rho = *m_Primtive[0];
		auto& u = *m_Primtive[1];
		auto& v = *m_Primtive[2];
		auto& w = *m_Primtive[3];
		auto& p = *m_Primtive[4];
		auto& cons0 = *m_Conservative[0];
		auto& cons1 = *m_Conservative[1];
		auto& cons2 = *m_Conservative[2];
		auto& cons3 = *m_Conservative[3];
		auto& cons4 = *m_Conservative[4];
		auto& inletPara = para->GetPrimitiveInflow();
		rho[boundIndex] = inletPara[0];
		u[boundIndex] = inletPara[1];
		v[boundIndex] = inletPara[2];
		w[boundIndex] = inletPara[3];
		p[boundIndex] = inletPara[4];
		Primitive2Conservative(rho[boundIndex], u[boundIndex], v[boundIndex], w[boundIndex], p[boundIndex],
			cons0[boundIndex], cons1[boundIndex], cons2[boundIndex], cons3[boundIndex], cons4[boundIndex]);
	}

	void NSSolver::ComputeOutletBC(Boundary& bound)
	{
		auto& rho = *m_Primtive[0];
		auto& u = *m_Primtive[1];
		auto& v = *m_Primtive[2];
		auto& w = *m_Primtive[3];
		auto& p = *m_Primtive[4];
		auto& cons0 = *m_Conservative[0];
		auto& cons1 = *m_Conservative[1];
		auto& cons2 = *m_Conservative[2];
		auto& cons3 = *m_Conservative[3];
		auto& cons4 = *m_Conservative[4];
		int boundIndex = bound.GetIndex();
		int innerIndex = bound.GetInnerNodeIndex();
		rho[boundIndex] = rho[innerIndex];
		u[boundIndex] = u[innerIndex];
		v[boundIndex] = v[innerIndex];
		w[boundIndex] = w[innerIndex];
		p[boundIndex] = p[innerIndex];
		Primitive2Conservative(rho[boundIndex], u[boundIndex], v[boundIndex], w[boundIndex], p[boundIndex],
			cons0[boundIndex], cons1[boundIndex], cons2[boundIndex], cons3[boundIndex], cons4[boundIndex]);


	}

	void NSSolver::ComputeWallBC(Boundary& bound)
	{
		int& innerIndex = bound.GetInnerNodeIndex();
		int boundIndex = bound.GetIndex();
		auto& rho = *m_Primtive[0];
		auto& u = *m_Primtive[1];
		auto& v = *m_Primtive[2];
		auto& w = *m_Primtive[3];
		auto& p = *m_Primtive[4];
		auto& cons0 = *m_Conservative[0];
		auto& cons1 = *m_Conservative[1];
		auto& cons2 = *m_Conservative[2];
		auto& cons3 = *m_Conservative[3];
		auto& cons4 = *m_Conservative[4];
		auto& boundNorm = bound.GetNorm();
		rho[boundIndex] = rho[innerIndex];
		u[boundIndex] = u[innerIndex];
		v[boundIndex] = v[innerIndex];
		w[boundIndex] = w[innerIndex];
		p[boundIndex] = p[innerIndex];
		DVector3D innerVel(u[innerIndex], v[innerIndex], w[innerIndex]);
		DVector3D boundVel = innerVel - (innerVel.dot(boundNorm)) * boundNorm / (boundNorm.norm() * boundNorm.norm());
		u[boundIndex] = boundVel(0);
		v[boundIndex] = boundVel(1);
		w[boundIndex] = boundVel(2);
		Primitive2Conservative(rho[boundIndex], u[boundIndex], v[boundIndex], w[boundIndex], p[boundIndex],
			cons0[boundIndex], cons1[boundIndex], cons2[boundIndex], cons3[boundIndex], cons4[boundIndex]);
	}

	void NSSolver::ComputeLimiterCoef()
	{
		ComputeLimiterCoefVK();
		return;
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& nodeNeighbor = nodeTopo->GetNeighborCloud();
		auto& prim = m_Primtive;
		auto& limiterCoef = m_LimiterCoef;
		auto& primGradX = m_PrimGradX;
		auto& primGradY = m_PrimGradY;
		auto& primGradZ = m_PrimGradZ;
		int nTotalNodeNum = grid->GetTotalNodeNum();
		double maxVal, minVal;
		for (int iNode = 0; iNode < nTotalNodeNum; ++iNode)
		{
			if (nodeType[iNode] != NodeType::inner)
				continue;
			auto& neighborNode = nodeNeighbor[iNode];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				maxVal = (*prim[iVal])[iNode];
				minVal = (*prim[iVal])[iNode];
				for (int iNeighbor = 0; iNeighbor < neighborNode.size(); ++iNeighbor)
				{
					maxVal = Max(maxVal, (*prim[iVal])[neighborNode[iNeighbor]]);
					minVal = Min(minVal, (*prim[iVal])[neighborNode[iNeighbor]]);
				}
				double gradx = (*primGradX[iVal])[iNode];
				double grady = (*primGradY[iVal])[iNode];
				double gradz = (*primGradZ[iVal])[iNode];
				double deltaMax = maxVal - (*prim[iVal])[iNode];
				double deltaMin = minVal - (*prim[iVal])[iNode];
				double tempCoef = LARGE_NUMBER;
				(*limiterCoef[iVal])[iNode] = LARGE_NUMBER;
				for (int iNeighbor = 0; iNeighbor < neighborNode.size(); ++iNeighbor)
				{
					auto current2Neighbor = nodeCoord[neighborNode[iNeighbor]] - nodeCoord[iNode];
					double delta2 = current2Neighbor(0) * gradx + current2Neighbor(1) * grady + current2Neighbor(2) * gradz;
					delta2 *= 0.5;
					if (delta2 > 0)
					{
						tempCoef = limiter(maxVal - (*prim[iVal])[iNode], delta2);
					}
					else if (delta2 < 0)
					{
						tempCoef = limiter(minVal - (*prim[iVal])[iNode], delta2);
					}
					else
					{
						tempCoef = 1.0;
					}
					(*limiterCoef[iVal])[iNode] = Min((*limiterCoef[iVal])[iNode], tempCoef);
				}
			}
		}
		ComputeBoundaryLimiterCoef();
	}
	void NSSolver::ComputeLimiterCoefVK()
	{
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& nodeNeighbor = nodeTopo->GetNeighborCloud();
		auto& prim = m_Primtive;
		auto& limiterCoef = m_LimiterCoef;
		auto& primGradX = m_PrimGradX;
		auto& primGradY = m_PrimGradY;
		auto& primGradZ = m_PrimGradZ;
		int nTotalNodeNum = grid->GetTotalNodeNum();
		double maxVal, minVal;
		double eps = 1e-6;
		double venkatCoeff = 1.0e-5;
		for (int iNode = 0; iNode < nTotalNodeNum; ++iNode)
		{
			if (nodeType[iNode] != NodeType::inner)
				continue;
			auto& neighborNode = nodeNeighbor[iNode];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				maxVal = (*prim[iVal])[iNode];
				minVal = (*prim[iVal])[iNode];
				for (int iNeighbor = 0; iNeighbor < neighborNode.size(); ++iNeighbor)
				{
					maxVal = Max(maxVal, (*prim[iVal])[neighborNode[iNeighbor]]);
					minVal = Min(minVal, (*prim[iVal])[neighborNode[iNeighbor]]);
				}
				eps = venkatCoeff * (maxVal - minVal);
				eps = eps * eps + SMALL_NUMBER;
				//eps = venkatCoeff * (maxVal - minVal) + SMALL_NUMBER;
				double gradx = (*primGradX[iVal])[iNode];
				double grady = (*primGradY[iVal])[iNode];
				double gradz = (*primGradZ[iVal])[iNode];
				double deltaMax = maxVal - (*prim[iVal])[iNode];
				double deltaMin = minVal - (*prim[iVal])[iNode];
				double tempCoef = LARGE_NUMBER;
				(*limiterCoef[iVal])[iNode] = LARGE_NUMBER;
				for (int iNeighbor = 0; iNeighbor < neighborNode.size(); ++iNeighbor)
				{
					auto current2Neighbor = nodeCoord[neighborNode[iNeighbor]] - nodeCoord[iNode];
					double delta2 = current2Neighbor(0) * gradx + current2Neighbor(1) * grady + current2Neighbor(2) * gradz;
					delta2 *= 0.5;
					if (delta2 > 0)
					{
						tempCoef = VenFun(maxVal - (*prim[iVal])[iNode], delta2, eps);
					}
					else if (delta2 < 0)
					{
						tempCoef = VenFun(minVal - (*prim[iVal])[iNode], delta2, eps);
					}
					else
					{
						tempCoef = 1.0;
					}
					(*limiterCoef[iVal])[iNode] = Min((*limiterCoef[iVal])[iNode], tempCoef);
				}

			}
		}
		ComputeBoundaryLimiterCoef();
	}
	void NSSolver::ComputeBoundaryLimiterCoef()
	{
		auto& grid = GetGrid();
		BoundaryMapPtr& boundaryMapPtr = grid->GetBoundaryMap();
		auto& boundaryMap = boundaryMapPtr->GetBoundaryMap();
		auto& limiterCoef = m_LimiterCoef;
		for (auto& boundary : boundaryMap)
		{
			auto& boundName = boundary.first;
			auto& bound = boundary.second;
			for (int iBound = 0; iBound < bound.size(); ++iBound)
			{
				auto& boundIndex = bound[iBound].GetIndex();
				auto& innerIndex = bound[iBound].GetInnerNodeIndex();
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					(*limiterCoef[iVal])[boundIndex] = 0;
					//(*limiterCoef[iVal])[boundIndex] = (*limiterCoef[iVal])[innerIndex];
				}
			}
		}
	}
	void NSSolver::CreateLimiter()
	{
		FlowSolverParaPtr& para = GetPara();
		LimiterType limiterType = para->GetLimiterType();
		if (limiterType == LimiterType::minmod)
			limiter = MinMod;
		else if (limiterType == LimiterType::vanleer)
			limiter = VanLeer;
		else if (limiterType == LimiterType::vanalbada)
			limiter = VanAlbada;
		else if (limiterType == LimiterType::mixminmodvanleer)
			limiter = MixMinModVanLeer;
		else if (limiterType == LimiterType::nolimit)
			limiter = NoLimiter;
		else if (limiterType == LimiterType::oneorder)
		{
			limiter = OneOrder;
		}
		else if (limiterType == LimiterType::barth)
		{
			limiter = Barth;
		}
		else
		{
			ZaranLog::error("Limiter Type Error!");
		}
	}

	void NSSolver::NoGradient()
	{
		// do nothing
	}

}