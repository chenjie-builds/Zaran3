#include "NSSolverFNFDM.h"
#include"File.h"
#include <fstream>
namespace zaran
{
	NSSolverFNFDM::NSSolverFNFDM(int index, string name, FlowSolverPara* para, GridFN* grid, FieldData* fieldData, DataManagerNS_FNFDM* data_manager)
		:NSSolver(index, name, para, grid, fieldData), m_data_manager(data_manager)
	{
		m_node_metric = new NodeMetric(grid);
		double** residual = new double* [GetEquNum()];
		for (int i = 0; i < GetEquNum(); i++)
		{
			residual[i] = m_data_manager->GetResidual(i);
		}
		m_grad_wlsq = new GradWLSQ(grid);
	}

	NSSolverFNFDM::~NSSolverFNFDM()
	{
		delete[] m_data_manager;
		delete[] m_node_metric;
	}

	void NSSolverFNFDM::InitFieldFarFlow()
	{
		GridFN* grid = GetGrid();
		FlowSolverPara* para = GetPara();
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
				m_data_manager->SetPrimitive(iVal, iNode, inflow_prim[iVal]);
		}
	}
	void NSSolverFNFDM::InitFieldFarFieldNoVelocity()
	{
		auto grid = GetGrid();
		FlowSolverPara* para = GetPara();
		int n_node = grid->GetTotalNodeNum();
		for (int iNode = 0; iNode < n_node; ++iNode)
		{
			m_data_manager->SetPrimitive(0, iNode, para->GetInflowDensity());
			m_data_manager->SetPrimitive(1, iNode, 0.0);
			m_data_manager->SetPrimitive(2, iNode, 0.0);
			m_data_manager->SetPrimitive(3, iNode, 0.0);
			m_data_manager->SetPrimitive(4, iNode, para->GetInflowPressure());
		}
	}
	void NSSolverFNFDM::InitFieldBackup()
	{
		auto grid = GetGrid();
		FlowSolverPara* para = GetPara();
		std::string restartFileName = para->GetBackupFieldFileName();
		std::ifstream fin(restartFileName);
		if (!fin.is_open())
		{
			Log::warn("Backup file not found!");
			exit(0);
		}
		int n_node = grid->GetTotalNodeNum();
		double value;
		for (int iNode = 0; iNode < n_node; ++iNode)
		{
			for (int iVal = 0; iVal < GetEquNum(); ++iVal)
			{
				fin >> value;
				m_data_manager->SetPrimitive(iVal, iNode, value);
			}
		}
		fin.close();
	}
	void NSSolverFNFDM::CalcMetric()
	{
		Log::info("Compute NS 3D Coordination Transformation Coefficients");
		m_node_metric->CalcMetric();
		m_node_metric->CheckJacobian();
		m_node_metric->FixJacobian();
		Log::info("max jacobi: {}, index: {}", m_node_metric->GetMaxJacobian(), m_node_metric->GetMaxJacobianNode());
		Log::info("min jacobi: {}, index: {}", m_node_metric->GetMinJacobian(), m_node_metric->GetMinJacobianNode());
		Log::info("NS 3D Coordination Transformation Coefficients are computed");
	}


	void NSSolverFNFDM::BackupField(std::string& back_folder)
	{
		auto grid = GetGrid();
		FlowSolverPara* para = GetPara();
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
				fout << m_data_manager->GetPrimitive(iVal, iNode) << " ";
			}
			fout << std::endl;
		}
		fout.close();
	}

	void NSSolverFNFDM::RungeKutta()
	{
		auto grid = GetGrid();
		FlowSolverPara* para = GetPara();
		const DArray& rk_coef = para->GetRKCoef();
		int rkStage = rk_coef.size();
		auto nodeTopo = grid->GetNodeTopo();
		int nInnerNode = grid->GetInnerNodeNum();
		int nBoundNode = grid->GetBoundNodeNum();
		double dt, jacobi;
		for (int iStage = 0; iStage < rkStage; ++iStage)
		{
			CalcResidual();
			#pragma omp parallel for private(dt, jacobi)
			for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
			{
				if (nodeTopo->GetType(iNode) != NodeType::inner)
					continue;
				dt = m_data_manager->GetTimeStep(iNode);
				jacobi = m_node_metric->GetJacobian(iNode);
				for (int iVal = 0; iVal < 5; ++iVal)
				{
					m_data_manager->SetConservative(iVal, iNode, m_data_manager->GetConservative(iVal, iNode) + dt * rk_coef[iStage] * m_data_manager->GetResidual(iVal, iNode) * jacobi);
				}
			}
		}
	}

	void NSSolverFNFDM::Prim2Cons()
	{
		auto grid = GetGrid();
		int n_node = grid->GetTotalNodeNum();
#pragma omp parallel for
		for (int iNode = 0; iNode < n_node; ++iNode)
		{
			double prim[5] = { m_data_manager->GetPrimitive(0, iNode), m_data_manager->GetPrimitive(1, iNode), m_data_manager->GetPrimitive(2, iNode),
							m_data_manager->GetPrimitive(3, iNode), m_data_manager->GetPrimitive(4, iNode) };
			double cons[5];
			NSSolver::Prim2Cons(prim, cons);
			for (int iVal = 0; iVal < 5; ++iVal)
			{
				m_data_manager->SetConservative(iVal, iNode, cons[iVal]);
			}
		}
	}

	void NSSolverFNFDM::Cons2Prim()
	{
		auto grid = GetGrid();
		int n_node = grid->GetTotalNodeNum();
#pragma omp parallel for
		for (int iNode = 0; iNode < n_node; ++iNode)
		{
			double cons[5];
			for (int iVal = 0; iVal < 5; ++iVal)
			{
				cons[iVal] = m_data_manager->GetConservative(iVal, iNode);
			}
			double prim[5];
			NSSolver::Cons2Prim(cons, prim);
			for (int iVal = 0; iVal < 5; ++iVal)
			{
				m_data_manager->SetPrimitive(iVal, iNode, prim[iVal]);
			}
		}
	}

	void NSSolverFNFDM::MidPointReconstruct2ndOrder(int index_left, int index_right, double* value_rec_left, double* value_rec_right)
	{
		auto grid = GetGrid();
		auto nodeTopo = grid->GetNodeTopo();
		double r[3];
		r[0] = nodeTopo->GetCoord(index_right)[0] - nodeTopo->GetCoord(index_left)[0];
		r[1] = nodeTopo->GetCoord(index_right)[1] - nodeTopo->GetCoord(index_left)[1];
		r[2] = nodeTopo->GetCoord(index_right)[2] - nodeTopo->GetCoord(index_left)[2];
		for (int iVal = 0; iVal < GetEquNum(); ++iVal)
		{
			value_rec_left[iVal] = m_data_manager->GetPrimitive(iVal, index_left) + 0.5 * m_data_manager->GetLimiter(iVal, index_left) *
				(r[0] * m_data_manager->GetPrimitiveGrad(iVal, 0, index_left) +
					r[1] * m_data_manager->GetPrimitiveGrad(iVal, 1, index_left) +
					r[2] * m_data_manager->GetPrimitiveGrad(iVal, 2, index_left));
			value_rec_right[iVal] = m_data_manager->GetPrimitive(iVal, index_right) - 0.5 * m_data_manager->GetLimiter(iVal, index_right) *
				(r[0] * m_data_manager->GetPrimitiveGrad(iVal, 0, index_right) +
					r[1] * m_data_manager->GetPrimitiveGrad(iVal, 1, index_right) +
					r[2] * m_data_manager->GetPrimitiveGrad(iVal, 2, index_right));
		}
	}

	void NSSolverFNFDM::MidPointReconstruct1stOrder(int index_left, int index_right, double* value_rec_left, double* value_rec_right)
	{
		auto  grid = GetGrid();
		auto nodeTopo = grid->GetNodeTopo();
		for (int iVal = 0; iVal < GetEquNum(); ++iVal)
		{
			value_rec_left[iVal] = m_data_manager->GetPrimitive(iVal, index_left);
			value_rec_right[iVal] = m_data_manager->GetPrimitive(iVal, index_right);
		}
	}

	void NSSolverFNFDM::BoundaryCondition()
	{
		auto grid = GetGrid();
		BoundaryMap* boundaryMap = grid->GetBoundaryMap();
		auto& wallBound = boundaryMap->GetBoundary("slipWall");
#pragma omp parallel for
		for (int iBound = 0; iBound < wallBound.size(); ++iBound)
			WallBC(wallBound[iBound]);
		auto& outletBound = boundaryMap->GetBoundary("outlet");
#pragma omp parallel for
		for (int iBound = 0; iBound < outletBound.size(); ++iBound)
			OutletBC(outletBound[iBound]);
		auto& inletBound = boundaryMap->GetBoundary("inlet");
#pragma omp parallel for
		for (int iBound = 0; iBound < inletBound.size(); ++iBound)
			InletBC(inletBound[iBound]);
	}
	void NSSolverFNFDM::InletBC(Boundary& bound)
	{
		FlowSolverPara* para = GetPara();
		int boundIndex = bound.GetIndex();
		double inflow_prim[5];
		inflow_prim[0] = para->GetInflowDensity();
		inflow_prim[1] = para->GetInflowVelocityX();
		inflow_prim[2] = para->GetInflowVelocityY();
		inflow_prim[3] = para->GetInflowVelocityZ();
		inflow_prim[4] = para->GetInflowPressure();
		double inflow_cons[5];
		NSSolver::Prim2Cons(inflow_prim, inflow_cons);
		for (int iVal = 0; iVal < 5; ++iVal)
		{
			m_data_manager->SetPrimitive(iVal, boundIndex, inflow_prim[iVal]);
			m_data_manager->SetConservative(iVal, boundIndex, inflow_cons[iVal]);
		}
	}

	void NSSolverFNFDM::OutletBC(Boundary& bound)
	{
		int boundIndex = bound.GetIndex();
		int innerIndex = bound.GetInnerIndex();
		double prim[5], cons[5];
		for (int iVal = 0; iVal < 5; ++iVal)
		{
			prim[iVal] = m_data_manager->GetPrimitive(iVal, innerIndex);
		}
		NSSolver::Prim2Cons(prim, cons);
		for (int iVal = 0; iVal < 5; ++iVal)
		{
			m_data_manager->SetPrimitive(iVal, boundIndex, prim[iVal]);
			m_data_manager->SetConservative(iVal, boundIndex, cons[iVal]);
		}
	}

	void NSSolverFNFDM::WallBC(Boundary& bound)
	{
		int& innerIndex = bound.GetInnerIndex();
		int boundIndex = bound.GetIndex();
		auto boundNorm = bound.GetNorm();
		double prim_bound[5];
		for (int iVal = 0; iVal < 5; ++iVal)
		{
			prim_bound[iVal] = m_data_manager->GetPrimitive(iVal, innerIndex);
		}
		double boundVel[3] = { prim_bound[1], prim_bound[2], prim_bound[3] };
		double vn = boundVel[0] * boundNorm[0] + boundVel[1] * boundNorm[1] + boundVel[2] * boundNorm[2];
		prim_bound[1] = prim_bound[1] - 1.0 * vn * boundNorm[0];
		prim_bound[2] = prim_bound[2] - 1.0 * vn * boundNorm[1];
		prim_bound[3] = prim_bound[3] - 1.0 * vn * boundNorm[2];
		double cons_bound[5];
		NSSolver::Prim2Cons(prim_bound, cons_bound);
		for (int iVal = 0; iVal < 5; ++iVal)
		{
			m_data_manager->SetPrimitive(iVal, boundIndex, prim_bound[iVal]);
			m_data_manager->SetConservative(iVal, boundIndex, cons_bound[iVal]);
		}
	}

	void NSSolverFNFDM::RiemannBC(Boundary& bound)
	{
		int id_in = bound.GetInnerIndex();
		int id_bound = bound.GetIndex();
		auto norm_bound = bound.GetNorm();
		double prim_in[5] = { m_data_manager->GetPrimitive(0, id_in), m_data_manager->GetPrimitive(1, id_in), m_data_manager->GetPrimitive(2, id_in),
							m_data_manager->GetPrimitive(3, id_in), m_data_manager->GetPrimitive(4, id_in) };
		double prim_bound[5] = { m_data_manager->GetPrimitive(0, id_bound), m_data_manager->GetPrimitive(1, id_bound), m_data_manager->GetPrimitive(2, id_bound),
								m_data_manager->GetPrimitive(3, id_bound), m_data_manager->GetPrimitive(4, id_bound) };
		FlowSolverPara* para = GetPara();
		Gas* gas = para->GetGas();
		double gamma = gas->GetGamma();
		double prim_far[5] = { para->GetInflowDensity(), para->GetInflowVelocityX(), para->GetInflowVelocityY(),
							  para->GetInflowVelocityZ(), para->GetInflowPressure() };
		double vel_in[3] = { prim_in[1], prim_in[2], prim_in[3] };
		// 速度在边界法向上的投影
		double vn_in = vel_in[0] * norm_bound[0] + vel_in[1] * norm_bound[1] + vel_in[2] * norm_bound[2];
		double vel_far[3] = { para->GetInflowVelocityX(), para->GetInflowVelocityY(), para->GetInflowVelocityZ() };
		double vn_far = vel_far[0] * norm_bound[0] + vel_far[1] * norm_bound[1] + vel_far[2] * norm_bound[2];
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
					m_data_manager->SetPrimitive(iVal, id_bound, prim_in[iVal]);
				}
			}
			else // 超声速入口
			{
				m_data_manager->SetPrimitive(0, id_bound, para->GetInflowDensity());
				m_data_manager->SetPrimitive(1, id_bound, para->GetInflowVelocityX());
				m_data_manager->SetPrimitive(2, id_bound, para->GetInflowVelocityY());
				m_data_manager->SetPrimitive(3, id_bound, para->GetInflowVelocityZ());
				m_data_manager->SetPrimitive(4, id_bound, para->GetInflowPressure());
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
				m_data_manager->SetPrimitive(0, id_bound, pow(c_bound * c_bound / (entropy * gamma), 1.0 / gamma1));
				m_data_manager->SetPrimitive(1, id_bound, prim_far[1] + norm_bound[0] * (vn_bound - vn_far));
				m_data_manager->SetPrimitive(2, id_bound, prim_far[2] + norm_bound[1] * (vn_bound - vn_far));
				m_data_manager->SetPrimitive(3, id_bound, prim_far[3] + norm_bound[2] * (vn_bound - vn_far));
				m_data_manager->SetPrimitive(4, id_bound, c_bound * c_bound * m_data_manager->GetPrimitive(0, id_bound) / gamma);
			}
			else // 亚声速出口
			{
				double entropy = prim_in[4] / pow(prim_in[0], gamma);
				m_data_manager->SetPrimitive(0, id_bound, pow(c_bound * c_bound / (entropy * gamma), 1.0 / gamma1));
				m_data_manager->SetPrimitive(1, id_bound, prim_in[1] + norm_bound[0] * (vn_bound - vn_in));
				m_data_manager->SetPrimitive(2, id_bound, prim_in[2] + norm_bound[1] * (vn_bound - vn_in));
				m_data_manager->SetPrimitive(3, id_bound, prim_in[3] + norm_bound[2] * (vn_bound - vn_in));
				m_data_manager->SetPrimitive(4, id_bound, c_bound * c_bound * m_data_manager->GetPrimitive(0, id_bound) / gamma);
			}
		}
		for (int iVal = 0; iVal < 5; ++iVal)
			prim_bound[iVal] = m_data_manager->GetPrimitive(iVal, id_bound);
		double cons_bound[5];
		NSSolver::Prim2Cons(prim_bound, cons_bound);
		for (int iVal = 0; iVal < 5; ++iVal)
		{
			m_data_manager->SetConservative(iVal, id_bound, cons_bound[iVal]);
		}
	}
	void NSSolverFNFDM::CalcPrimGradBound()
	{
		auto grid = GetGrid();
		auto nodeTopo = grid->GetNodeTopo();
		BoundaryMap* boundaryMapPtr = grid->GetBoundaryMap();
		auto& boundaryMap = boundaryMapPtr->GetBoundaryMap();
		for (auto& boundary : boundaryMap)
		{
			auto& boundName = boundary.first;
			auto& bound = boundary.second;
			if (boundName == "hole")
				continue;
#pragma omp parallel for
			for (int iVal = 0; iVal < GetEquNum(); ++iVal)
			{
				for (int iBound = 0; iBound < bound.size(); ++iBound)
				{
					auto& boundIndex = bound[iBound].GetIndex();
					auto& innerIndex = bound[iBound].GetInnerIndex();
					for (int iDim = 0; iDim < 3; ++iDim)
					{
						m_data_manager->SetPrimitiveGrad(iVal, iDim, boundIndex, 0.0);
					}
				}
			}
		}
	}

	void NSSolverFNFDM::CalcGradWLS()
	{
		auto grid = GetGrid();
		auto nodeTopo = grid->GetNodeTopo();
		int nInnerNode = grid->GetInnerNodeNum();
		int nBoundNode = grid->GetBoundNodeNum();
		Matrix3d A, A_inv;
		DVector3D b, grad;
		double omega = 0;
		double deltaVal;
		double deltaX, deltaY, deltaZ;
#pragma omp parallel for private(A,A_inv, b, grad, omega, deltaVal, deltaX, deltaY, deltaZ)
		for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
		{
			if (nodeTopo->GetType(iNode) != NodeType::inner)
				continue;
			int neighbor_num = nodeTopo->GetNeighborNodeNum(iNode);
			auto neighborNodeVec = nodeTopo->GetNeighborNode(iNode);
			auto nodeCoord = nodeTopo->GetCoord(iNode);
			A.setZero();
			for (size_t iNeib = 0; iNeib < neighbor_num; ++iNeib)
			{
				auto neighborCoord = nodeTopo->GetCoord(neighborNodeVec[iNeib]);
				omega = DistanceOfTwoPoints(nodeCoord, neighborCoord);
				if (abs(omega) < SMALL_NUMBER)
					continue;
				omega = 1.0 / omega;
				deltaX = neighborCoord[0] - nodeCoord[0];
				deltaY = neighborCoord[1] - nodeCoord[1];
				deltaZ = neighborCoord[2] - nodeCoord[2];
				A(0, 0) += omega * deltaX * deltaX;
				A(0, 1) += omega * deltaX * deltaY;
				A(0, 2) += omega * deltaX * deltaZ;
				A(1, 0) += omega * deltaY * deltaX;
				A(1, 1) += omega * deltaY * deltaY;
				A(1, 2) += omega * deltaY * deltaZ;
				A(2, 0) += omega * deltaZ * deltaX;
				A(2, 1) += omega * deltaZ * deltaY;
				A(2, 2) += omega * deltaZ * deltaZ;
			}
			A_inv = A.inverse();
			for (size_t iVal = 0; iVal < GetEquNum(); ++iVal)
			{
				b.setZero();
				for (size_t iNeib = 0; iNeib < neighbor_num; ++iNeib)
				{
					auto neighborCoord = nodeTopo->GetCoord(neighborNodeVec[iNeib]);
					omega = DistanceOfTwoPoints(nodeCoord, neighborCoord);
					if (abs(omega) < SMALL_NUMBER)
						continue;
					omega = 1.0 / omega;
					deltaVal = m_data_manager->GetPrimitive(iVal, neighborNodeVec[iNeib]) - m_data_manager->GetPrimitive(iVal, iNode);
					deltaX = neighborCoord[0] - nodeCoord[0];
					deltaY = neighborCoord[1] - nodeCoord[1];
					deltaZ = neighborCoord[2] - nodeCoord[2];
					b(0) += omega * deltaVal * deltaX;
					b(1) += omega * deltaVal * deltaY;
					b(2) += omega * deltaVal * deltaZ;
				}
				grad = A_inv * b;
				for (size_t iDim = 0; iDim < 3; ++iDim)
				{
					m_data_manager->SetPrimitiveGrad(iVal, iDim, iNode, grad(iDim));
				}
			}
		}
	}

	void NSSolverFNFDM::CalcLimiter()
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

	void NSSolverFNFDM::CalcLimiterVK()
	{
		auto grid = GetGrid();
		auto nodeTopo = grid->GetNodeTopo();
		int nTotalNodeNum = grid->GetTotalNodeNum();
		double maxVal, minVal;
		double eps = 1e-6;
		double venkatCoeff = 1.0e-5;
		for (int iVal = 0; iVal < GetEquNum(); ++iVal)
		{
#pragma omp parallel for private(maxVal, minVal, eps)
			for (int iNode = 0; iNode < nTotalNodeNum; ++iNode)
			{
				if (nodeTopo->GetType(iNode) != NodeType::inner && nodeTopo->GetType(iNode) != NodeType::hole)
					continue;
				auto currentNodeCoord = nodeTopo->GetCoord(iNode);
				auto neighborNode = nodeTopo->GetNeighborNode(iNode);
				int nNeighbor = nodeTopo->GetNeighborNodeNum(iNode);
				maxVal = m_data_manager->GetPrimitive(iVal, iNode);
				minVal = m_data_manager->GetPrimitive(iVal, iNode);
				for (int iNeighbor = 0; iNeighbor < nNeighbor; ++iNeighbor)
				{
					maxVal = Max(maxVal, m_data_manager->GetPrimitive(iVal, neighborNode[iNeighbor]));
					minVal = Min(minVal, m_data_manager->GetPrimitive(iVal, neighborNode[iNeighbor]));
				}
				eps = venkatCoeff * (maxVal - minVal);
				eps = eps * eps + SMALL_NUMBER;
				// eps = venkatCoeff * (maxVal - minVal) + SMALL_NUMBER;
				double deltaMax = maxVal - m_data_manager->GetPrimitive(iVal, iNode);
				double deltaMin = minVal - m_data_manager->GetPrimitive(iVal, iNode);
				double tempCoef = LARGE_NUMBER;
				m_data_manager->SetLimiter(iVal, iNode, LARGE_NUMBER);
				for (int iNeighbor = 0; iNeighbor < nNeighbor; ++iNeighbor)
				{
					auto neighborNodeCoord = nodeTopo->GetCoord(neighborNode[iNeighbor]);
					double delta2 = 0.0;
					for (int iDim = 0; iDim < 3; ++iDim)
						delta2 += (neighborNodeCoord[iDim] - currentNodeCoord[iDim]) * m_data_manager->GetPrimitiveGrad(iVal, iDim, iNode);
					delta2 *= 0.5;
					if (delta2 > 0)
					{
						tempCoef = LimiterVK(maxVal - m_data_manager->GetPrimitive(iVal, iNode), delta2, eps);
					}
					else if (delta2 < 0)
					{
						tempCoef = LimiterVK(minVal - m_data_manager->GetPrimitive(iVal, iNode), delta2, eps);
					}
					else
					{
						tempCoef = 1.0;
					}
					m_data_manager->SetLimiter(iVal, iNode, Min(m_data_manager->GetLimiter(iNode, iVal), tempCoef));
				}
			}
		}
	}

	void NSSolverFNFDM::CalcLimiterBJ()
	{
		auto grid = GetGrid();
		auto nodeTopo = grid->GetNodeTopo();
		int nTotalNodeNum = grid->GetTotalNodeNum();
		double maxVal, minVal;
#pragma omp parallel for private(maxVal, minVal)
		for (int iVal = 0; iVal < GetEquNum(); ++iVal)
		{
			for (int iNode = 0; iNode < nTotalNodeNum; ++iNode)
			{
				if (nodeTopo->GetType(iNode) != NodeType::inner && nodeTopo->GetType(iNode) != NodeType::hole)
					continue;
				auto currentNodeCoord = nodeTopo->GetCoord(iNode);
				auto neighborNode = nodeTopo->GetNeighborNode(iNode);
				auto nNeighbor = nodeTopo->GetNeighborNodeNum(iNode);
				maxVal = m_data_manager->GetPrimitive(iVal, iNode);
				minVal = m_data_manager->GetPrimitive(iVal, iNode);
				for (int iNeighbor = 0; iNeighbor < nNeighbor; ++iNeighbor)
				{
					maxVal = Max(maxVal, m_data_manager->GetPrimitive(iVal, neighborNode[iNeighbor]));
					minVal = Min(minVal, m_data_manager->GetPrimitive(iVal, neighborNode[iNeighbor]));
				}
				double deltaMax = maxVal - m_data_manager->GetPrimitive(iVal, iNode);
				double deltaMin = minVal - m_data_manager->GetPrimitive(iVal, iNode);
				double tempCoef = LARGE_NUMBER;
				m_data_manager->SetLimiter(iVal, iNode, LARGE_NUMBER);
				for (int iNeighbor = 0; iNeighbor < nNeighbor; ++iNeighbor)
				{
					auto neighborNodeCoord = nodeTopo->GetCoord(neighborNode[iNeighbor]);
					double delta2 = 0.0;
					for (int iDim = 0; iDim < 3; ++iDim)
						delta2 += (neighborNodeCoord[iDim] - currentNodeCoord[iDim]) * m_data_manager->GetPrimitiveGrad(iVal, iDim, iNode);
					delta2 *= 0.5;
					if (delta2 > 0)
					{
						tempCoef = LimiterBarth(maxVal - m_data_manager->GetPrimitive(iVal, iNode), delta2);
					}
					else if (delta2 < 0)
					{
						tempCoef = LimiterBarth(minVal - m_data_manager->GetPrimitive(iVal, iNode), delta2);
					}
					else
					{
						tempCoef = 1.0;
					}
					m_data_manager->SetLimiter(iVal, iNode, Min(m_data_manager->GetLimiter(iNode, iVal), tempCoef));
				}
			}
		}
	}

	void NSSolverFNFDM::CalcLimiterNone()
	{
		auto grid = GetGrid();
		auto nodeTopo = grid->GetNodeTopo();
		int nTotalNodeNum = grid->GetTotalNodeNum();
		for (int iVal = 0; iVal < GetEquNum(); ++iVal)
		{
#pragma omp parallel for
			for (int iNode = 0; iNode < nTotalNodeNum; ++iNode)
			{
				if (nodeTopo->GetType(iNode) != NodeType::inner && nodeTopo->GetType(iNode) != NodeType::hole)
					continue;
				m_data_manager->SetLimiter(iVal, iNode, 1.0);
			}
		}
	}

	void NSSolverFNFDM::CalcLimiterFirstOrder()
	{
		auto grid = GetGrid();
		auto nodeTopo = grid->GetNodeTopo();
		int nTotalNodeNum = grid->GetTotalNodeNum();
		for (int iVal = 0; iVal < GetEquNum(); ++iVal)
		{
#pragma omp parallel for
			for (int iNode = 0; iNode < nTotalNodeNum; ++iNode)
			{
				m_data_manager->SetLimiter(iVal, iNode, 0.0);
			}
		}
	}

	void NSSolverFNFDM::CalcLimiterBound()
	{
		auto grid = GetGrid();
		auto boundaryMapPtr = grid->GetBoundaryMap();
		auto& boundaryMap = boundaryMapPtr->GetBoundaryMap();
		for (auto& boundary : boundaryMap)
		{
			auto& boundName = boundary.first;
			auto& bound = boundary.second;
			if (boundName == "hole")
				continue;
			for (int iVal = 0; iVal < GetEquNum(); ++iVal)
			{
#pragma omp parallel for
				for (int iBound = 0; iBound < bound.size(); ++iBound)
				{
					auto& boundIndex = bound[iBound].GetIndex();
					auto& innerIndex = bound[iBound].GetInnerIndex();
					m_data_manager->SetLimiter(iVal, boundIndex, 0.0);
				}
			}
		}
	}

	void NSSolverFNFDM::CheckPrimtive()
	{
		auto grid = GetGrid();
		auto nodeTopo = grid->GetNodeTopo();
		int nTotalNodeNum = grid->GetTotalNodeNum();
		int equation_num = GetEquNum();
		int nonphysical_node_num = 0;
#pragma omp parallel for  reduction(+ : nonphysical_node_num)
		for (int iNode = 0; iNode < nTotalNodeNum; ++iNode)
		{
			bool exist_nonphysical = false;
			m_data_manager->SetNonPhysical(iNode, -1);
			if (nodeTopo->GetType(iNode) != NodeType::inner && nodeTopo->GetType(iNode) != NodeType::hole)
				continue;
			if (m_data_manager->GetDensity(iNode) < 0 || m_data_manager->GetPressure(iNode) < 0)
			{
				exist_nonphysical = true;
			}
			if (!exist_nonphysical)
			{
				for (int iVal = 0; iVal < equation_num; ++iVal)
				{
					if (isnan(m_data_manager->GetPrimitive(iVal, iNode)) || isinf(m_data_manager->GetPrimitive(iVal, iNode)))
					{
						exist_nonphysical = true;
						break;
					}
				}
			}
			if (exist_nonphysical)
			{
				m_data_manager->SetNonPhysical(iNode, 1);
				nonphysical_node_num++;
				Log::info("Non-physical Node: {}, neighbor num: {}, prim: {},{},{},{},{}", iNode,
					nodeTopo->GetNeighborNodeNum(iNode), m_data_manager->GetPrimitive(0, iNode), m_data_manager->GetPrimitive(1, iNode),
					m_data_manager->GetPrimitive(2, iNode), m_data_manager->GetPrimitive(3, iNode), m_data_manager->GetPrimitive(4, iNode));
				Log::info("Non-physical Node: {}, coord: {},{},{},", iNode, nodeTopo->GetCoord(iNode)[0],
					nodeTopo->GetCoord(iNode)[1], nodeTopo->GetCoord(iNode)[2]);
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

	void NSSolverFNFDM::CheckResidual()
	{

	}

	void NSSolverFNFDM::FixPrimtive()
	{
		auto grid = GetGrid();
		auto node_topo = grid->GetNodeTopo();
		int total_node_num = grid->GetTotalNodeNum();
		int equation_num = GetEquNum();
		DArray weight, distance;
		IArray physical_neighbor;
		double sum = 0;
#pragma omp parallel for private(physical_neighbor, weight, distance, sum)
		for (int iNode = 0; iNode < total_node_num; ++iNode)
		{

			if (node_topo->GetType(iNode) != NodeType::inner && node_topo->GetType(iNode) != NodeType::hole)
				continue;
			if (m_data_manager->GetNonPhysical(iNode) == -1)
				continue;
			physical_neighbor.clear();
			auto current_node_coord = node_topo->GetCoord(iNode);
			auto neighborNode = node_topo->GetNeighborNode(iNode);
			int nNeighbor = node_topo->GetNeighborNodeNum(iNode);
			for (int iNeighbor = 0; iNeighbor < nNeighbor; ++iNeighbor)
			{
				if (m_data_manager->GetNonPhysical(neighborNode[iNeighbor]) == 1)
					continue;
				physical_neighbor.push_back(neighborNode[iNeighbor]);
			}
			weight.resize(physical_neighbor.size());
			distance.resize(physical_neighbor.size());
			sum = 0;
			for (int iNeighbor = 0; iNeighbor < physical_neighbor.size(); ++iNeighbor)
			{
				if (m_data_manager->GetNonPhysical(physical_neighbor[iNeighbor]) == 1)
					continue;
				auto neighbor_node_coord = node_topo->GetCoord(physical_neighbor[iNeighbor]);
				distance[iNeighbor] = DistanceOfTwoPoints(current_node_coord, neighbor_node_coord);
				sum += 1.0 / distance[iNeighbor];
			}
			for (int iNeighbor = 0; iNeighbor < physical_neighbor.size(); ++iNeighbor)
			{
				if (m_data_manager->GetNonPhysical(physical_neighbor[iNeighbor]) == 1)
					continue;
				weight[iNeighbor] = 1.0 / (distance[iNeighbor] * sum);
			}
			for (int iVal = 0; iVal < GetEquNum(); ++iVal)
			{
				if (isnan(m_data_manager->GetPrimitive(iVal, iNode)) || isinf(m_data_manager->GetPrimitive(iVal, iNode)))
				{
					m_data_manager->SetPrimitive(iVal, iNode, 0);
					for (int iNeighbor = 0; iNeighbor < physical_neighbor.size(); ++iNeighbor)
					{
						if (m_data_manager->GetNonPhysical(physical_neighbor[iNeighbor]) == 1)
							continue;
						m_data_manager->SetPrimitive(iVal, iNode, m_data_manager->GetPrimitive(iVal, physical_neighbor[iNeighbor]) * weight[iNeighbor]);
					}
				}
			}
			if (m_data_manager->GetPrimitive(0, iNode) < 0)
				m_data_manager->SetPrimitive(0, iNode, 0);
			for (int iNeighbor = 0; iNeighbor < physical_neighbor.size(); ++iNeighbor)
			{
				if (m_data_manager->GetNonPhysical(physical_neighbor[iNeighbor]) == 1)
					continue;
				m_data_manager->SetPrimitive(0, iNode, m_data_manager->GetPrimitive(0, iNode) + m_data_manager->GetPrimitive(0, physical_neighbor[iNeighbor]) * weight[iNeighbor]);
			}
			if (m_data_manager->GetPrimitive(0, iNode) < 0)
				m_data_manager->SetPrimitive(0, iNode, 0);
			for (int iNeighbor = 0; iNeighbor < physical_neighbor.size(); ++iNeighbor)
			{
				if (m_data_manager->GetNonPhysical(physical_neighbor[iNeighbor]) == 1)
					continue;
				m_data_manager->SetPrimitive(4, iNode, m_data_manager->GetPrimitive(1, iNode) + m_data_manager->GetPrimitive(4, physical_neighbor[iNeighbor]) * weight[iNeighbor]);
			}
			double prim[5] = { m_data_manager->GetPrimitive(0, iNode),m_data_manager->GetPrimitive(1, iNode),m_data_manager->GetPrimitive(2, iNode),m_data_manager->GetPrimitive(3, iNode),m_data_manager->GetPrimitive(4, iNode) };
			double cons[5];
			NSSolver::Prim2Cons(prim, cons);
			for (int iVal = 0; iVal < 5; ++iVal)
			{
				m_data_manager->SetConservative(iVal, iNode, cons[iVal]);
				m_data_manager->SetResidual(iVal, iNode, 0.0);
			}
		}
	}

	void NSSolverFNFDM::CalcTimeStepLocal()
	{
		auto grid = GetGrid();
		auto nodeTopo = grid->GetNodeTopo();
		FlowSolverPara* para = GetPara();
		double gamma = para->GetGas()->GetGamma();
		double cfl = para->GetCflNumber();
		int nInnerNode = grid->GetInnerNodeNum();
		double min_dt = LARGE_NUMBER;
#pragma omp parallel for reduction(min:min_dt)
		for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
		{
			if (nodeTopo->GetType(iNode) != NodeType::inner)
				continue;
			auto xi = m_node_metric->GetMetricXi(iNode);
			auto eta = m_node_metric->GetMetricEta(iNode);
			auto zeta = m_node_metric->GetMetricZeta(iNode);
			auto jacobi = m_node_metric->GetJacobian(iNode);
			double c = sqrt(gamma * m_data_manager->GetPressure(iNode) / m_data_manager->GetDensity(iNode));
			double norm_xi = sqrt(xi[0] * xi[0] + xi[1] * xi[1] + xi[2] * xi[2]);
			double norm_eta = sqrt(eta[0] * eta[0] + eta[1] * eta[1] + eta[2] * eta[2]);
			double norm_zeta = sqrt(zeta[0] * zeta[0] + zeta[1] * zeta[1] + zeta[2] * zeta[2]);
			double u_xi = m_data_manager->GetVelocity(0, iNode) * xi[0] + m_data_manager->GetVelocity(1, iNode) * xi[1] + m_data_manager->GetVelocity(2, iNode) * xi[2];
			double u_eta = m_data_manager->GetVelocity(0, iNode) * eta[0] + m_data_manager->GetVelocity(1, iNode) * eta[1] + m_data_manager->GetVelocity(2, iNode) * eta[2];
			double u_zeta = m_data_manager->GetVelocity(0, iNode) * zeta[0] + m_data_manager->GetVelocity(1, iNode) * zeta[1] + m_data_manager->GetVelocity(2, iNode) * zeta[2];
			double lamda = abs(u_xi) + abs(u_eta) + abs(u_zeta) + c * (norm_xi + norm_eta + norm_zeta);
			m_data_manager->SetTimeStep(iNode, cfl / lamda);
			if (min_dt > m_data_manager->GetTimeStep(iNode))
			{
				min_dt = m_data_manager->GetTimeStep(iNode);
			}
		}
		GlobalData::Update("dt", min_dt);
	}

	void NSSolverFNFDM::SnycTimeStepWithGlobal(double& dt)
	{
		auto grid = GetGrid();
		auto para = GetPara();
#pragma omp parallel for
		for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
		{
			m_data_manager->SetTimeStep(iNode, dt);
		}
	}

	void NSSolverFNFDM::ZeroResidual()
	{
		auto grid = GetGrid();
		auto nodeTopo = grid->GetNodeTopo();
		int nNode = grid->GetTotalNodeNum();
		for (int iVar = 0; iVar < GetEquNum(); ++iVar)
		{
#pragma omp parallel for
			for (int iNode = 0; iNode < nNode; ++iNode)
			{
				m_data_manager->SetResidual(iVar, iNode, 0.0);
			}
		}
	}

	void NSSolverFNFDM::ConvectiveResidual()
	{
		auto grid = GetGrid();
		auto nodeTopo = grid->GetNodeTopo();
		RiemannSolverPara riemann_para[6];
		for (int i = 0; i < 6; ++i)
		{
			riemann_para[i].gamma_left = riemann_para[i].gamma_right = GetPara()->GetGas()->GetGamma();
		}
		bool exist_negative = false;
#pragma omp parallel for private(riemann_para,exist_negative)
		for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
		{
			if (nodeTopo->GetType(iNode) != NodeType::inner)
				continue;
			exist_negative = false;
			double jacobi = m_node_metric->GetJacobian(iNode);
			auto neighbor = nodeTopo->GetNeighborNode(iNode);
			// i direction
			riemann_para[0].norm(0) = m_node_metric->GetMetricXi(iNode)[0];
			riemann_para[0].norm(1) = m_node_metric->GetMetricXi(iNode)[1];
			riemann_para[0].norm(2) = m_node_metric->GetMetricXi(iNode)[2];
			riemann_para[0].nt = m_node_metric->GetMetricXi(iNode)[3];
			riemann_para[1].norm(0) = m_node_metric->GetMetricXi(iNode)[0];
			riemann_para[1].norm(1) = m_node_metric->GetMetricXi(iNode)[1];
			riemann_para[1].norm(2) = m_node_metric->GetMetricXi(iNode)[2];
			riemann_para[1].nt = m_node_metric->GetMetricXi(iNode)[3];
			MidPointReconstruct2ndOrder(iNode, neighbor[1], &riemann_para[0].prim_left(0), &riemann_para[0].prim_right(0));
			MidPointReconstruct2ndOrder(neighbor[0], iNode, &riemann_para[1].prim_left(0), &riemann_para[1].prim_right(0));
			// j direction
			riemann_para[2].norm(0) = m_node_metric->GetMetricEta(iNode)[0];
			riemann_para[2].norm(1) = m_node_metric->GetMetricEta(iNode)[1];
			riemann_para[2].norm(2) = m_node_metric->GetMetricEta(iNode)[2];
			riemann_para[2].nt = m_node_metric->GetMetricEta(iNode)[3];
			riemann_para[3].norm(0) = m_node_metric->GetMetricEta(iNode)[0];
			riemann_para[3].norm(1) = m_node_metric->GetMetricEta(iNode)[1];
			riemann_para[3].norm(2) = m_node_metric->GetMetricEta(iNode)[2];
			riemann_para[3].nt = m_node_metric->GetMetricEta(iNode)[3];
			MidPointReconstruct2ndOrder(iNode, neighbor[3], &riemann_para[2].prim_left(0), &riemann_para[2].prim_right(0));
			MidPointReconstruct2ndOrder(neighbor[2], iNode, &riemann_para[3].prim_left(0), &riemann_para[3].prim_right(0));

			// k direction
			riemann_para[4].norm(0) = m_node_metric->GetMetricZeta(iNode)[0];
			riemann_para[4].norm(1) = m_node_metric->GetMetricZeta(iNode)[1];
			riemann_para[4].norm(2) = m_node_metric->GetMetricZeta(iNode)[2];
			riemann_para[4].nt = m_node_metric->GetMetricZeta(iNode)[3];
			riemann_para[5].norm(0) = m_node_metric->GetMetricZeta(iNode)[0];
			riemann_para[5].norm(1) = m_node_metric->GetMetricZeta(iNode)[1];
			riemann_para[5].norm(2) = m_node_metric->GetMetricZeta(iNode)[2];
			riemann_para[5].nt = m_node_metric->GetMetricZeta(iNode)[3];
			MidPointReconstruct2ndOrder(iNode, neighbor[5], &riemann_para[4].prim_left(0), &riemann_para[4].prim_right(0));
			MidPointReconstruct2ndOrder(neighbor[4], iNode, &riemann_para[5].prim_left(0), &riemann_para[5].prim_right(0));
			// check negative density and pressure
			for (int i = 0; i < 6; ++i)
			{
				if (riemann_para[i].prim_left(0) < 0 || riemann_para[i].prim_right(0) < 0 || riemann_para[i].prim_left(4) < 0 || riemann_para[i].prim_right(4) < 0)
				{
					exist_negative = true;
					break;
				}
			}
			if (exist_negative || m_data_manager->GetNonPhysical(iNode) > 0)
			{
				MidPointReconstruct1stOrder(iNode, neighbor[1], &riemann_para[0].prim_left(0), &riemann_para[0].prim_right(0));
				MidPointReconstruct1stOrder(neighbor[0], iNode, &riemann_para[1].prim_left(0), &riemann_para[1].prim_right(0));
				MidPointReconstruct1stOrder(iNode, neighbor[3], &riemann_para[2].prim_left(0), &riemann_para[2].prim_right(0));
				MidPointReconstruct1stOrder(neighbor[2], iNode, &riemann_para[3].prim_left(0), &riemann_para[3].prim_right(0));
				MidPointReconstruct1stOrder(iNode, neighbor[5], &riemann_para[4].prim_left(0), &riemann_para[4].prim_right(0));
				MidPointReconstruct1stOrder(neighbor[4], iNode, &riemann_para[5].prim_left(0), &riemann_para[5].prim_right(0));
			}
			// calculate flux
			for (int i = 0; i < 6; ++i)
			{
				m_riemann_solver->Solver(riemann_para[i]);
			}
			for (int iVar = 0; iVar < GetEquNum(); ++iVar)
			{
				double flux = (riemann_para[0].flux[iVar] - riemann_para[1].flux[iVar] + riemann_para[2].flux[iVar] - riemann_para[3].flux[iVar] + riemann_para[4].flux[iVar] - riemann_para[5].flux[iVar]) / jacobi;
				m_data_manager->SetResidual(iVar, iNode, m_data_manager->GetResidual(iVar, iNode) - flux);
			}
		}
	}

	void NSSolverFNFDM::ViscousResidual()
	{
		return;
		CalcViscousFlux();
		CalcViscousFluxGrad();
		auto grid = GetGrid();
		auto nodeTopo = grid->GetNodeTopo();
		int nNode = grid->GetTotalNodeNum();
		for (int iNode = 0; iNode < nNode; ++iNode)
		{
			if (nodeTopo->GetType(iNode) != NodeType::inner)
				continue;
			for (int iVar = 0; iVar < GetEquNum(); ++iVar)
			{
				m_data_manager->SetResidual(iVar, iNode, m_data_manager->GetResidual(iVar, iNode) - m_data_manager->GetViscousFluxGrad(iVar, 0, 0, iNode) - m_data_manager->GetViscousFluxGrad(iVar, 1, 1, iNode) - m_data_manager->GetViscousFluxGrad(iVar, 2, 2, iNode));
			}
		}
	}
	void NSSolverFNFDM::CalcViscousFlux()
	{
		return;
		FlowSolverPara* para = GetPara();
		Gas* gas = para->GetGas();
		auto grid = GetGrid();
		auto nodeTopo = grid->GetNodeTopo();
		int nNode = grid->GetTotalNodeNum();
		double* viscous_flux_x, * viscous_flux_y, * viscous_flux_z;
		double tau_xx, tau_yy, tau_zz, tau_xy, tau_xz, tau_yz;
		double grad_rho[3], grad_u[3], grad_v[3], grad_w[3], grad_p[3];
		double vel[3];
		double tempeture;
		double vis_coef;
		double therm_coef;
		double lamda;
		for (int iNode = 0; iNode < nNode; ++iNode)
		{
			vel[0] = m_data_manager->GetPrimitive(1, iNode);
			vel[1] = m_data_manager->GetPrimitive(2, iNode);
			vel[2] = m_data_manager->GetPrimitive(3, iNode);
			for (int iDim = 0; iDim < 3; ++iDim)
			{
				grad_rho[iDim] = m_data_manager->GetPrimitiveGrad(0, iDim, iNode);
				grad_u[iDim] = m_data_manager->GetPrimitiveGrad(1, iDim, iNode);
				grad_v[iDim] = m_data_manager->GetPrimitiveGrad(2, iDim, iNode);
				grad_w[iDim] = m_data_manager->GetPrimitiveGrad(3, iDim, iNode);
				grad_p[iDim] = m_data_manager->GetPrimitiveGrad(4, iDim, iNode);
			}
			tempeture = gas->CalcTemperature(m_data_manager->GetDensity(iNode), m_data_manager->GetPressure(iNode));
			vis_coef = gas->CalcMu(tempeture);
			lamda = 2.0 / 3.0 * vis_coef;
			therm_coef = gas->CalcK(tempeture);
			tau_xx = 2 * vis_coef * grad_u[0] - lamda * (grad_u[0] + grad_v[1] + grad_w[2]);
			tau_yy = 2 * vis_coef * grad_v[1] - lamda * (grad_u[0] + grad_v[1] + grad_w[2]);
			tau_zz = 2 * vis_coef * grad_w[2] - lamda * (grad_u[0] + grad_v[1] + grad_w[2]);
			tau_xy = vis_coef * (grad_u[1] + grad_v[0]);
			tau_xz = vis_coef * (grad_u[2] + grad_w[0]);
			tau_yz = vis_coef * (grad_v[2] + grad_w[1]);
			m_data_manager->SetViscousFlux(0, 0, iNode, 0);
			m_data_manager->SetViscousFlux(1, 0, iNode, tau_xx);
			m_data_manager->SetViscousFlux(2, 0, iNode, tau_xy);
			m_data_manager->SetViscousFlux(3, 0, iNode, tau_xz);
			m_data_manager->SetViscousFlux(4, 0, iNode, vel[0] * tau_xx + vel[1] * tau_xy + vel[2] * tau_xz + therm_coef * grad_rho[0]);
			m_data_manager->SetViscousFlux(0, 1, iNode, 0);
			m_data_manager->SetViscousFlux(1, 1, iNode, tau_xy);
			m_data_manager->SetViscousFlux(2, 1, iNode, tau_yy);
			m_data_manager->SetViscousFlux(3, 1, iNode, tau_yz);
			m_data_manager->SetViscousFlux(4, 1, iNode, vel[0] * tau_xy + vel[1] * tau_yy + vel[2] * tau_yz + therm_coef * grad_rho[1]);
			m_data_manager->SetViscousFlux(0, 2, iNode, 0);
			m_data_manager->SetViscousFlux(1, 2, iNode, tau_xz);
			m_data_manager->SetViscousFlux(2, 2, iNode, tau_yz);
			m_data_manager->SetViscousFlux(3, 2, iNode, tau_zz);
			m_data_manager->SetViscousFlux(4, 2, iNode, vel[0] * tau_xz + vel[1] * tau_yz + vel[2] * tau_zz + therm_coef * grad_rho[2]);
		}
	}
	void NSSolverFNFDM::CalcViscousFluxGrad()
	{
		auto grid = GetGrid();
		auto nodeTopo = grid->GetNodeTopo();
		Matrix3d A, A_inv;
		DVector3D b, grad;
		double omega = 0;
		double deltaVal;
		double deltaX, deltaY, deltaZ;


	}
	void NSSolverFNFDM::SourceTermResidual()
	{
	}

	void NSSolverFNFDM::CalcForce()
	{
		return;
		auto grid = GetGrid();
		auto faceTopo = grid->GetFaceTopo();
		FlowSolverPara* para = GetPara();
		const Dimensionless& dimensionless = para->GetDimensionless();
		int nFace = faceTopo->GetFaceNum();
		double face_pressure;
		double force[3] = { 0 };
		for (int iFace = 0; iFace < nFace; ++iFace)
		{
			face_pressure = 0;
			int* face2node = faceTopo->GetFace2Node(iFace);
			for (int iNode = 0; iNode < faceTopo->GetFaceNodeNum(iFace); ++iNode)
			{
				face_pressure += m_data_manager->GetPressure(face2node[iNode]);
			}
			face_pressure /= faceTopo->GetFaceNodeNum(iFace);
			face_pressure -= dimensionless.GetPressureDL(0);
			force[0] -= face_pressure * faceTopo->GetArea(iFace) * faceTopo->GetNormal(iFace)[0];
			force[1] -= face_pressure * faceTopo->GetArea(iFace) * faceTopo->GetNormal(iFace)[1];
			force[2] -= face_pressure * faceTopo->GetArea(iFace) * faceTopo->GetNormal(iFace)[2];
		}
		force[0] = dimensionless.GetForceReal(force[0]);
		force[1] = dimensionless.GetForceReal(force[1]);
		force[2] = dimensionless.GetForceReal(force[2]);
		Log::info("Force: {}, {}, {}", force[0], force[1], force[2]);
	}


	GridFN* NSSolverFNFDM::GetGrid()
	{
		return static_cast<GridFN*>(FlowSolver::GetGrid());
	}
}