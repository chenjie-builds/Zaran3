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
		ZaranLog::info("Flow Field Initialize Finished!");
	}

	void NSSolver::InitSolver()
	{
		ZaranLog::info("Initialize NS Solver!");
		FlowSolver::InitSolver();
		SetNumberOfEquations(5);
		std::string riemannSolverType = GlobalData::GetString("riemannSolver");
		RiemannSolverFactory riemannSolverFactory;
		riemannSolverFactory.Create(riemannSolver_, riemannSolverType);
		ZaranLog::info("NS Solver Initialize Finished!");
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
		data.AddData("nonPhysical", emptyData);
	}
	void NSSolver::RegisterFieldData()
	{

		auto& data = *GetFieldData();

		// 预锟斤拷锟斤拷锟节达�?
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
		m_non_physical = &data.GetData("nonPhysical");

	}

	void NSSolver::Solve()
	{
		ComputeTimeStep();
		ComputePrimitiveGradient();
		ComputeLimiterCoef();
		BoundaryCondition();
		TimeAdvance();
		UpdateField();
		CheckPrimtive();
		FixPrimtive();
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
		//锟角凤拷使锟斤拷全锟斤拷时锟戒步
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

	void NSSolver::MidPointReconstruct(int index_left, int index_right, double* value_rec_left, double* value_rec_right)
	{
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& prim = m_Primitive;
		auto& cons = m_Conservative;
		auto& primGradX = m_PrimGradX;
		auto& primGradY = m_PrimGradY;
		auto& primGradZ = m_PrimGradZ;
		auto& limiterCoef = m_LimiterCoef;
		DVector3D r = nodeCoord[index_right] - nodeCoord[index_left];
		for (int iVal = 0;iVal < GetNumberOfEquations();++iVal)
		{
			value_rec_left[iVal] = (*prim[iVal])[index_left] +
				0.5 * (*limiterCoef[iVal])[index_left] *
				(r.x() * (*primGradX[iVal])[index_left] +
					r.y() * (*primGradY[iVal])[index_left] +
					r.z() * (*primGradZ[iVal])[index_left]);
			value_rec_right[iVal] = (*prim[iVal])[index_right] -
				0.5 * (*limiterCoef[iVal])[index_right] *
				(r.x() * (*primGradX[iVal])[index_right] +
					r.y() * (*primGradY[iVal])[index_right] +
					r.z() * (*primGradZ[iVal])[index_right]);
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
				value_rec_left[iVal] = (*prim[iVal])[index_left];
				value_rec_right[iVal] = (*prim[iVal])[index_right];
			}
			OutputError(index_left);
		}
		if (value_rec_right[0] < 0 || value_rec_right[4] < 0)
		{
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				value_rec_left[iVal] = (*prim[iVal])[index_left];
				value_rec_right[iVal] = (*prim[iVal])[index_right];
			}
			OutputError(index_right);
		}
#endif
	}

	void NSSolver::MidPointReconstructOneOrder(int index_left, int index_right, double* value_rec_left, double* value_rec_right)
	{
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& prim = m_Primitive;
		auto& cons = m_Conservative;
		auto& primGradX = m_PrimGradX;
		auto& primGradY = m_PrimGradY;
		auto& primGradZ = m_PrimGradZ;
		auto& limiterCoef = m_LimiterCoef;
		DVector3D r = nodeCoord[index_right] - nodeCoord[index_left];
		for (int iVal = 0;iVal < GetNumberOfEquations();++iVal)
		{
			value_rec_left[iVal] = (*prim[iVal])[index_left];
			value_rec_right[iVal] = (*prim[iVal])[index_right];
		}
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
		auto& nodeTopo = grid->GetNodeTopo();
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
				int boundNeighborNum = nodeTopo->GetNeighborCloud()[boundIndex].size();
				// if (boundNeighborNum <= 6)
				// 	ZaranLog::info("Boundary {} has {} neighbors", boundIndex, boundNeighborNum);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					(*primGradX[iVal])[boundIndex] = 0;
					(*primGradY[iVal])[boundIndex] = 0;
					(*primGradZ[iVal])[boundIndex] = 0;
					// (*primGradX[iVal])[boundIndex] = (*primGradX[iVal])[innerIndex];
					// (*primGradY[iVal])[boundIndex] = (*primGradY[iVal])[innerIndex];
					// (*primGradZ[iVal])[boundIndex] = (*primGradZ[iVal])[innerIndex];
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
		if (limiterType != "oneOrder")
		{
			int firstOrderSteps = GlobalData::GetInt("firstOrderSteps");
			int currentStep = GlobalData::GetInt("step");
			if (currentStep < firstOrderSteps)
			{
				limiterType = "oneOrder";
				ZaranLog::info("First {}/{} steps use one order scheme", currentStep, firstOrderSteps);
			}
		}
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
		//ComputeBoundaryLimiterCoef();
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
	void NSSolver::CheckPrimtive()
	{
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& nodeNeighbor = nodeTopo->GetNeighborCloud();
		auto& prim = m_Primitive;
		auto& res = m_Residual;
		auto& limiterCoef = m_LimiterCoef;
		int nTotalNodeNum = grid->GetTotalNodeNum();
		int equation_num = GetNumberOfEquations();
		DArray ave_prim(equation_num, 0.0);
		int nonphysical_node_num = 0;
		// #pragma omp parallel for private(ave_prim) reduction(+:nonphysical_node_num)
		for (int iNode = 0; iNode < nTotalNodeNum; ++iNode)
		{
			bool exist_nonphysical = false;
			(*m_non_physical)[iNode] = -1.0;
			if (nodeType[iNode] != NodeType::inner && nodeType[iNode] != NodeType::hole)
				continue;
			if ((*prim[0])[iNode] < 0 || (*prim[4])[iNode] < 0)
			{
				exist_nonphysical = true;
			}
			if (!exist_nonphysical)
			{
				for (int iVal = 0; iVal < equation_num; ++iVal)
				{
					if (isnan((*prim[iVal])[iNode]) || isinf((*prim[iVal])[iNode]))
					{
						exist_nonphysical = true;
						break;
					}
				}
			}
			if (exist_nonphysical)
			{
				(*m_non_physical)[iNode] = 1.0;
				nonphysical_node_num++;
				ZaranLog::info("Non-physical Node: {}, neighbor num: {}, prim: {},{},{},{},{}", iNode, nodeNeighbor[iNode].size(), (*prim[0])[iNode], (*prim[1])[iNode], (*prim[2])[iNode], (*prim[3])[iNode], (*prim[4])[iNode]);
				ZaranLog::info("Non-physical Node: {}, coord: {},{},{},", iNode, nodeCoord[iNode].x(), nodeCoord[iNode].y(), nodeCoord[iNode].z());
			}
		}
		if (nonphysical_node_num > 0)
		{
			ZaranLog::warn("Non-physical Node Num: {}", nonphysical_node_num);
			auto& para = GetPara();
			double cfl = para->GetCflNumber();
			cfl = cfl / 10.0;
			para->SetCflNumber(cfl);
			ZaranLog::warn("CFL Number is reduced to {}", cfl);
		}
		else
		{
			auto& para = GetPara();
			double cfl = para->GetCflNumber();
			double cfl_max = GlobalData::GetDouble("cflNumber");
			if (cfl < cfl_max)
			{
				cfl = Min(cfl_max, cfl * 1.15);
				para->SetCflNumber(cfl);
				ZaranLog::info("CFL Number is increased to {}", cfl);
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
		auto& prim = m_Primitive;
		auto& res = m_Residual;
		auto& limiterCoef = m_LimiterCoef;
		int total_node_num = grid->GetTotalNodeNum();
		int equation_num = GetNumberOfEquations();
		DArray weight, distance;
		IArray physical_neighbor;
		double sum = 0;
		// #pragma omp parallel for private(ave_prim)
		for (int iNode = 0; iNode < total_node_num; ++iNode)
		{
			if (node_type[iNode] != NodeType::inner && node_type[iNode] != NodeType::hole)
				continue;
			if ((*m_non_physical)[iNode] < 0)
				continue;
			physical_neighbor.clear();
			auto& neighborNode = node_neighbor[iNode];
			for (int iNeighbor = 0; iNeighbor < neighborNode.size(); ++iNeighbor)
			{
				if ((*m_non_physical)[neighborNode[iNeighbor]] > 0)
					continue;
				physical_neighbor.push_back(neighborNode[iNeighbor]);
			}
			weight.resize(physical_neighbor.size());
			distance.resize(physical_neighbor.size());
			sum = 0;
			for (int iNeighbor = 0; iNeighbor < physical_neighbor.size(); ++iNeighbor)
			{
				if ((*m_non_physical)[physical_neighbor[iNeighbor]] > 0)
					continue;
				distance[iNeighbor] = (node_coord[physical_neighbor[iNeighbor]] - node_coord[iNode]).norm();
				sum += 1.0 / distance[iNeighbor];
			}
			for (int iNeighbor = 0; iNeighbor < physical_neighbor.size(); ++iNeighbor)
			{
				if ((*m_non_physical)[physical_neighbor[iNeighbor]] > 0)
					continue;
				weight[iNeighbor] = 1.0 / (distance[iNeighbor] * sum);
			}
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				if (isnan((*prim[iVal])[iNode]) || isinf((*prim[iVal])[iNode]))
				{
					(*prim[iVal])[iNode] = 0;
					for (int iNeighbor = 0; iNeighbor < physical_neighbor.size(); ++iNeighbor)
					{
						if ((*m_non_physical)[physical_neighbor[iNeighbor]] > 0)
							continue;
						(*prim[iVal])[iNode] += (*prim[iVal])[physical_neighbor[iNeighbor]] * weight[iNeighbor];
					}
				}
			}
			if ((*prim[0])[iNode] < 0)
				(*prim[0])[iNode] = 0;
			for (int iNeighbor = 0; iNeighbor < physical_neighbor.size(); ++iNeighbor)
			{
				if ((*m_non_physical)[physical_neighbor[iNeighbor]] > 0)
					continue;
				(*prim[0])[iNode] += (*prim[0])[physical_neighbor[iNeighbor]] * weight[iNeighbor];
			}
			if ((*prim[4])[iNode] < 0)
				(*prim[4])[iNode] = 0;
			for (int iNeighbor = 0; iNeighbor < physical_neighbor.size(); ++iNeighbor)
			{
				if ((*m_non_physical)[physical_neighbor[iNeighbor]] > 0)
					continue;
				(*prim[4])[iNode] += (*prim[4])[physical_neighbor[iNeighbor]] * weight[iNeighbor];
			}
			Primitive2Conservative((*prim[0])[iNode], (*prim[1])[iNode], (*prim[2])[iNode], (*prim[3])[iNode], (*prim[4])[iNode],
				(*m_Conservative[0])[iNode], (*m_Conservative[1])[iNode], (*m_Conservative[2])[iNode], (*m_Conservative[3])[iNode], (*m_Conservative[4])[iNode]);
			(*res[0])[iNode] = 0;
		}

	}
	void NSSolver::NoGradient()
	{
		// do nothing
	}

}