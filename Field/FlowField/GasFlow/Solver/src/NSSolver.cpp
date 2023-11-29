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
		auto& rho = *m_Primitive[0];
		auto& u = *m_Primitive[1];
		auto& v = *m_Primitive[2];
		auto& w = *m_Primitive[3];
		auto& p = *m_Primitive[4];
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
			ZaranLog::warn("Initialize Failed!");
			ZaranLog::warn("Wrong Flow field init parameter:{}", initType);
			exit(0);
		}
		int n_data = rho.size();
		double x, y, z;
		for (int iNode = 0; iNode < n_data; ++iNode)
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
	}

	void NSSolver::CreateFieldData()
	{

		auto& grid = *GetGrid();
		int nTotalNodeNum = grid.GetTotalNodeNum();
		DArray emptyData(nTotalNodeNum);
		auto& dataPtr = GetFieldData();
		// dataPtr = std::make_shared<FieldData>();
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

		// Ԥ�����ڴ�
		m_Primitive.reserve(5);
		m_Conservative.reserve(5);
		m_Residual.reserve(5);
		m_LimiterCoef.reserve(5);
		m_PrimGradX.reserve(5);
		m_PrimGradY.reserve(5);
		m_PrimGradZ.reserve(5);
		m_ConservativeRK.reserve(5);
		m_CoordTrans.reserve(36);

		auto addDataToVector = [&data](Array<DArray*>& vec, const std::string& name) {
			vec.push_back(&data.GetData(name));
			};

		addDataToVector(m_Primitive, "rho");
		addDataToVector(m_Primitive, "u");
		addDataToVector(m_Primitive, "v");
		addDataToVector(m_Primitive, "w");
		addDataToVector(m_Primitive, "p");
		addDataToVector(m_Conservative, "cons0");
		addDataToVector(m_Conservative, "cons1");
		addDataToVector(m_Conservative, "cons2");
		addDataToVector(m_Conservative, "cons3");
		addDataToVector(m_Conservative, "cons4");
		addDataToVector(m_Residual, "res0");
		addDataToVector(m_Residual, "res1");
		addDataToVector(m_Residual, "res2");
		addDataToVector(m_Residual, "res3");
		addDataToVector(m_Residual, "res4");
		addDataToVector(m_LimiterCoef, "limiterCoef0");
		addDataToVector(m_LimiterCoef, "limiterCoef1");
		addDataToVector(m_LimiterCoef, "limiterCoef2");
		addDataToVector(m_LimiterCoef, "limiterCoef3");
		addDataToVector(m_LimiterCoef, "limiterCoef4");
		addDataToVector(m_PrimGradX, "rhoGradX");
		addDataToVector(m_PrimGradX, "uGradX");
		addDataToVector(m_PrimGradX, "vGradX");
		addDataToVector(m_PrimGradX, "wGradX");
		addDataToVector(m_PrimGradX, "pGradX");
		addDataToVector(m_PrimGradY, "rhoGradY");
		addDataToVector(m_PrimGradY, "uGradY");
		addDataToVector(m_PrimGradY, "vGradY");
		addDataToVector(m_PrimGradY, "wGradY");
		addDataToVector(m_PrimGradY, "pGradY");
		addDataToVector(m_PrimGradZ, "rhoGradZ");
		addDataToVector(m_PrimGradZ, "uGradZ");
		addDataToVector(m_PrimGradZ, "vGradZ");
		addDataToVector(m_PrimGradZ, "wGradZ");
		addDataToVector(m_PrimGradZ, "pGradZ");
		addDataToVector(m_ConservativeRK, "consRK0");
		addDataToVector(m_ConservativeRK, "consRK1");
		addDataToVector(m_ConservativeRK, "consRK2");
		addDataToVector(m_ConservativeRK, "consRK3");
		addDataToVector(m_ConservativeRK, "consRK4");
		addDataToVector(m_CoordTrans, "coordTransXXi");
		addDataToVector(m_CoordTrans, "coordTransXEta");
		addDataToVector(m_CoordTrans, "coordTransXZeta");
		addDataToVector(m_CoordTrans, "coordTransXTau");
		addDataToVector(m_CoordTrans, "coordTransYXi");
		addDataToVector(m_CoordTrans, "coordTransYEta");
		addDataToVector(m_CoordTrans, "coordTransYZeta");
		addDataToVector(m_CoordTrans, "coordTransYTau");
		addDataToVector(m_CoordTrans, "coordTransZXi");
		addDataToVector(m_CoordTrans, "coordTransZEta");
		addDataToVector(m_CoordTrans, "coordTransZZeta");
		addDataToVector(m_CoordTrans, "coordTransZTau");
		addDataToVector(m_CoordTrans, "coordTransTXi");
		addDataToVector(m_CoordTrans, "coordTransTEta");
		addDataToVector(m_CoordTrans, "coordTransTZeta");
		addDataToVector(m_CoordTrans, "coordTransTTau");
		addDataToVector(m_CoordTrans, "coordTransXiX");
		addDataToVector(m_CoordTrans, "coordTransXiY");
		addDataToVector(m_CoordTrans, "coordTransXiZ");
		addDataToVector(m_CoordTrans, "coordTransXiT");
		addDataToVector(m_CoordTrans, "coordTransEtaX");
		addDataToVector(m_CoordTrans, "coordTransEtaY");
		addDataToVector(m_CoordTrans, "coordTransEtaZ");
		addDataToVector(m_CoordTrans, "coordTransEtaT");
		addDataToVector(m_CoordTrans, "coordTransZetaX");
		addDataToVector(m_CoordTrans, "coordTransZetaY");
		addDataToVector(m_CoordTrans, "coordTransZetaZ");
		addDataToVector(m_CoordTrans, "coordTransZetaT");
		addDataToVector(m_CoordTrans, "coordTransTauX");
		addDataToVector(m_CoordTrans, "coordTransTauY");
		addDataToVector(m_CoordTrans, "coordTransTauZ");
		addDataToVector(m_CoordTrans, "coordTransTauT");
		addDataToVector(m_CoordTrans, "coordTransJ");
		m_TimeStep = &data.GetData("dt");

	}

	void NSSolver::Solve()
	{
		ComputeTimeStep();
		ComputePrimitiveGradient();
		ComputeLimiterCoef();
		BoundaryCondition();
		TimeAdvance();
		UpdateField();
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
		//�Ƿ�ʹ��ȫ��ʱ�䲽
		int useGlobalTimeStep = GlobalData::GetInt("useGlobalTimeStep");
		if (useGlobalTimeStep == 1)
		{
			double globlaTimeStep = GlobalData::GetDouble("dt");
			double globalTime = GlobalData::GetDouble("globalTime");
			double endTime = GlobalData::GetDouble("endTime");
			if (globalTime + globlaTimeStep > endTime)
			{
				globlaTimeStep = endTime - globalTime;
				globalTime = endTime;
			}
			else
				globalTime += globlaTimeStep;
			GlobalData::Update("globalTime", globalTime);
			SnycTimeStepWithGlobal(globlaTimeStep);
		}
	}

	void NSSolver::SnycTimeStepWithGlobal(double& dt)
	{
		GridPtr grid = GetGrid();
		FlowSolverParaPtr para = GetPara();
		for (auto var : *m_TimeStep)
		{
			(*m_TimeStep)[var] = dt;
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
		auto& inletBound = boundaryMap["inlet"];
		for (int iBound = 0; iBound < inletBound.size(); ++iBound)
			ComputeInletBC(inletBound[iBound]);

	}


	void NSSolver::ComputePrimitiveGradient()
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
			if (boundName == "hole")
				continue;
			for (int iBound = 0; iBound < bound.size(); ++iBound)
			{
				auto& boundIndex = bound[iBound].GetIndex();
				auto& innerIndex = bound[iBound].GetInnerIndex();
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
		int n_data = rho.size();
#pragma omp parallel for
		for (int iNode = 0; iNode < n_data; ++iNode)
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
		auto& inflowPrim = GetPara()->GetPrimitiveInflow();
#pragma omp parallel for
		for (int iNode = 0; iNode < rho.size(); ++iNode)
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
		int& innerIndex = bound.GetInnerIndex();
		int boundIndex = bound.GetIndex();
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
		string limiterType = GlobalData::GetString("limiterType");
		if (limiterType == "vk")
			ComputeLimiterCoefVK();
		else if (limiterType == "barth")
			ComputeLimiterCoefBJ();
		else if (limiterType == "noLimiter")
			ComputeLimiterCoefNoLimiter();
		else if (limiterType == "oneOrder")
			ComputeLimiterCoefOneOrder();
		else
			ZaranLog::warn("Unsupported Limiter Type: {}", limiterType);
		ComputeBoundaryLimiterCoef();
	}
	void NSSolver::ComputeLimiterCoefVK()
	{
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& nodeNeighbor = nodeTopo->GetNeighborCloud();
		auto& prim = m_Primitive;
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
			if (nodeType[iNode] != NodeType::inner && nodeType[iNode] != NodeType::hole)
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
	}
	void NSSolver::ComputeLimiterCoefBJ()
	{
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& nodeNeighbor = nodeTopo->GetNeighborCloud();
		auto& prim = m_Primitive;
		auto& limiterCoef = m_LimiterCoef;
		auto& primGradX = m_PrimGradX;
		auto& primGradY = m_PrimGradY;
		auto& primGradZ = m_PrimGradZ;
		int nTotalNodeNum = grid->GetTotalNodeNum();
		double maxVal, minVal;
		for (int iNode = 0; iNode < nTotalNodeNum; ++iNode)
		{
			if (nodeType[iNode] != NodeType::inner && nodeType[iNode] != NodeType::hole)
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
						tempCoef = Barth(maxVal - (*prim[iVal])[iNode], delta2);
					}
					else if (delta2 < 0)
					{
						tempCoef = Barth(minVal - (*prim[iVal])[iNode], delta2);
					}
					else
					{
						tempCoef = 1.0;
					}
					(*limiterCoef[iVal])[iNode] = Min((*limiterCoef[iVal])[iNode], tempCoef);
				}
			}
		}
	}
	void NSSolver::ComputeLimiterCoefNoLimiter()
	{
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		auto& limiterCoef = m_LimiterCoef;
		int nTotalNodeNum = grid->GetTotalNodeNum();
		for (int iNode = 0; iNode < nTotalNodeNum; ++iNode)
		{
			if (nodeType[iNode] != NodeType::inner && nodeType[iNode] != NodeType::hole)
				continue;
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				(*limiterCoef[iVal])[iNode] = 1.0;
			}
		}
	}
	void NSSolver::ComputeLimiterCoefOneOrder()
	{
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& nodeNeighbor = nodeTopo->GetNeighborCloud();
		auto& prim = m_Primitive;
		auto& limiterCoef = m_LimiterCoef;
		auto& primGradX = m_PrimGradX;
		auto& primGradY = m_PrimGradY;
		auto& primGradZ = m_PrimGradZ;
		int nTotalNodeNum = grid->GetTotalNodeNum();
		double maxVal, minVal;
		for (int iNode = 0; iNode < nTotalNodeNum; ++iNode)
		{
			if (nodeType[iNode] != NodeType::inner && nodeType[iNode] != NodeType::hole)
				continue;
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				(*limiterCoef[iVal])[iNode] = 0.0;
			}
		}
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
			if (boundName == "hole")
				continue;
			for (int iBound = 0; iBound < bound.size(); ++iBound)
			{
				auto& boundIndex = bound[iBound].GetIndex();
				auto& innerIndex = bound[iBound].GetInnerIndex();
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					(*limiterCoef[iVal])[boundIndex] = 0;
					//(*limiterCoef[iVal])[boundIndex] = (*limiterCoef[iVal])[innerIndex];
				}
			}
		}
	}
	void NSSolver::NoGradient()
	{
		// do nothing
	}

}