#include "NSSolverFNFDM.h"
#include"File.h"
#include <fstream>
namespace zaran
{
	NSSolverFNFDM::NSSolverFNFDM(int index, string name, FlowSolverPara* para, GridFN* grid, FieldData* fieldData, DataManagerNS_FNFDM* data_manager)
		:NSSolver(index, name, para, grid, fieldData), m_data_manager(data_manager)
	{
		m_node_metric = new NodeMetric(grid);
		m_grad_wlsq = new GradWLSQ(grid);
	}

	NSSolverFNFDM::~NSSolverFNFDM()
	{
		delete[] m_node_metric;
		delete[] m_grad_wlsq;
	}

	void NSSolverFNFDM::InitFieldFarFlow()
	{
		GridFN* grid = GetGrid();
		FlowSolverPara* para = GetPara();
		double prim_far[5];
		prim_far[0] = para->GetInflowDensity();
		prim_far[1] = para->GetInflowVelocityX();
		prim_far[2] = para->GetInflowVelocityY();
		prim_far[3] = para->GetInflowVelocityZ();
		prim_far[4] = para->GetInflowPressure();
		int n_node = grid->GetTotalNodeNum();
		for (int iVal = 0; iVal < para->GetEquNum(); ++iVal)
		{
			for (int iNode = 0; iNode < n_node; ++iNode)
				m_data_manager->SetPrimitive(iVal, iNode, prim_far[iVal]);
		}
	}
	void NSSolverFNFDM::InitFieldFarFieldNoVelocity()
	{
		auto grid = GetGrid();
		auto para = GetPara();
		int node_num = grid->GetTotalNodeNum();
		for (int iNode = 0; iNode < node_num; ++iNode)
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
		auto para = GetPara();
		std::string restart_file_name = para->GetBackupFieldFileName();
		std::ifstream fin(restart_file_name);
		if (!fin.is_open())
		{
			Log::warn("Backup file not found!");
			exit(0);
		}
		int node_num = grid->GetTotalNodeNum();
		double value;
		for (int iNode = 0; iNode < node_num; ++iNode)
		{
			for (int iVal = 0; iVal < para->GetEquNum(); ++iVal)
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
		auto para = GetPara();
		std::string backup_file_name = para->GetBackupFieldFileName();
		backup_file_name = back_folder + "/" + backup_file_name;
		if (IsFileExist(backup_file_name))
		{
			DeleteFile(backup_file_name);
		}
		std::ofstream fout(backup_file_name);
		int n_node = grid->GetTotalNodeNum();
		for (int iNode = 0; iNode < n_node; ++iNode)
		{
			for (int iVal = 0; iVal < para->GetEquNum(); ++iVal)
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
		auto node = grid->GetNode();
		int nInnerNode = grid->GetInnerNodeNum();
		int nBoundNode = grid->GetBoundNodeNum();
		double dt, jacobi;
		for (int iStage = 0; iStage < rkStage; ++iStage)
		{
			CalcResidual();
#pragma omp parallel for private(dt, jacobi)
			for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
			{
				if (node->GetType(iNode) != NodeType::inner)
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
		auto gas = GetGas();
		int node_num = grid->GetTotalNodeNum();
#pragma omp parallel for
		for (int iNode = 0; iNode < node_num; ++iNode)
		{
			double prim[5] = { m_data_manager->GetPrimitive(0, iNode), m_data_manager->GetPrimitive(1, iNode), m_data_manager->GetPrimitive(2, iNode),
							m_data_manager->GetPrimitive(3, iNode), m_data_manager->GetPrimitive(4, iNode) };
			double cons[5];
			gas->Prim2Cons(prim, cons);
			for (int iVal = 0; iVal < 5; ++iVal)
			{
				m_data_manager->SetConservative(iVal, iNode, cons[iVal]);
			}
		}
	}

	void NSSolverFNFDM::Cons2Prim()
	{
		auto grid = GetGrid();
		int node_num = grid->GetTotalNodeNum();
#pragma omp parallel for
		for (int iNode = 0; iNode < node_num; ++iNode)
		{
			double cons[5];
			for (int iVal = 0; iVal < 5; ++iVal)
			{
				cons[iVal] = m_data_manager->GetConservative(iVal, iNode);
			}
			double prim[5];
			GetGas()->Cons2Prim(cons, prim);
			for (int iVal = 0; iVal < 5; ++iVal)
			{
				m_data_manager->SetPrimitive(iVal, iNode, prim[iVal]);
			}
		}
	}

	void NSSolverFNFDM::MidPointReconstruct2ndOrder(int index_left, int index_right, double* value_rec_left, double* value_rec_right)
	{
		int equ_num = GetPara()->GetEquNum();
		auto grid = GetGrid();
		auto node = grid->GetNode();
		double vec_node2neighbor[3];
		vec_node2neighbor[0] = node->GetCoord(index_right)[0] - node->GetCoord(index_left)[0];
		vec_node2neighbor[1] = node->GetCoord(index_right)[1] - node->GetCoord(index_left)[1];
		vec_node2neighbor[2] = node->GetCoord(index_right)[2] - node->GetCoord(index_left)[2];
		for (int iVal = 0; iVal < equ_num; ++iVal)
		{
			value_rec_left[iVal] = m_data_manager->GetPrimitive(iVal, index_left) + 0.5 * m_data_manager->GetLimiter(iVal, index_left) *
				(vec_node2neighbor[0] * m_data_manager->GetPrimitiveGrad(iVal, 0, index_left) +
					vec_node2neighbor[1] * m_data_manager->GetPrimitiveGrad(iVal, 1, index_left) +
					vec_node2neighbor[2] * m_data_manager->GetPrimitiveGrad(iVal, 2, index_left));
			value_rec_right[iVal] = m_data_manager->GetPrimitive(iVal, index_right) - 0.5 * m_data_manager->GetLimiter(iVal, index_right) *
				(vec_node2neighbor[0] * m_data_manager->GetPrimitiveGrad(iVal, 0, index_right) +
					vec_node2neighbor[1] * m_data_manager->GetPrimitiveGrad(iVal, 1, index_right) +
					vec_node2neighbor[2] * m_data_manager->GetPrimitiveGrad(iVal, 2, index_right));
		}
	}

	void NSSolverFNFDM::MidPointReconstruct1stOrder(int index_left, int index_right, double* value_rec_left, double* value_rec_right)
	{
		int equ_num = GetPara()->GetEquNum();
		auto  grid = GetGrid();
		auto node = grid->GetNode();
		for (int iVal = 0; iVal < equ_num; ++iVal)
		{
			value_rec_left[iVal] = m_data_manager->GetPrimitive(iVal, index_left);
			value_rec_right[iVal] = m_data_manager->GetPrimitive(iVal, index_right);
		}
	}

	void NSSolverFNFDM::BoundaryCondition()
	{
		auto grid = GetGrid();
		BoundaryMap* bound_map = grid->GetBoundaryMap();
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
	void NSSolverFNFDM::InletBC(Boundary& bound)
	{
		FlowSolverPara* para = GetPara();
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
			m_data_manager->SetPrimitive(iVal, bound_index, prim_far[iVal]);
			m_data_manager->SetConservative(iVal, bound_index, cons_far[iVal]);
		}
	}

	void NSSolverFNFDM::OutletBC(Boundary& bound)
	{
		int bound_index = bound.GetIdxBound();
		int inner_index = bound.GetIdxRef();
		double prim[5], cons[5];
		for (int iVal = 0; iVal < 5; ++iVal)
		{
			prim[iVal] = m_data_manager->GetPrimitive(iVal, inner_index);
		}
		GetGas()->Prim2Cons(prim, cons);
		for (int iVal = 0; iVal < 5; ++iVal)
		{
			m_data_manager->SetPrimitive(iVal, bound_index, prim[iVal]);
			m_data_manager->SetConservative(iVal, bound_index, cons[iVal]);
		}
	}

	void NSSolverFNFDM::WallBC(Boundary& bound)
	{
		int idx_ref = bound.GetIdxRef();
		int idx_bound = bound.GetIdxBound();
		auto norm_bound = bound.GetNormBound();
		double prim_bound[5];
		for (int iVal = 0; iVal < 5; ++iVal)
		{
			prim_bound[iVal] = m_data_manager->GetPrimitive(iVal, idx_ref);
		}
		double bound_Vel[3] = { prim_bound[1], prim_bound[2], prim_bound[3] };
		double vn = bound_Vel[0] * norm_bound[0] + bound_Vel[1] * norm_bound[1] + bound_Vel[2] * norm_bound[2];
		prim_bound[1] = prim_bound[1] - 1.0 * vn * norm_bound[0];
		prim_bound[2] = prim_bound[2] - 1.0 * vn * norm_bound[1];
		prim_bound[3] = prim_bound[3] - 1.0 * vn * norm_bound[2];
		double cons_bound[5];
		GetGas()->Prim2Cons(prim_bound, cons_bound);
		for (int iVal = 0; iVal < 5; ++iVal)
		{
			m_data_manager->SetPrimitive(iVal, idx_bound, prim_bound[iVal]);
			m_data_manager->SetConservative(iVal, idx_bound, cons_bound[iVal]);
		}
	}

	void NSSolverFNFDM::RiemannBC(Boundary& bound)
	{
		int idx_ref = bound.GetIdxRef();
		int idx_bound = bound.GetIdxBound();
		auto norm_bound = bound.GetNormBound();
		double prim_in[5] = { m_data_manager->GetPrimitive(0, idx_ref), m_data_manager->GetPrimitive(1, idx_ref), m_data_manager->GetPrimitive(2, idx_ref),
							m_data_manager->GetPrimitive(3, idx_ref), m_data_manager->GetPrimitive(4, idx_ref) };
		double prim_bound[5] = { m_data_manager->GetPrimitive(0, idx_bound), m_data_manager->GetPrimitive(1, idx_bound), m_data_manager->GetPrimitive(2, idx_bound),
								m_data_manager->GetPrimitive(3, idx_bound), m_data_manager->GetPrimitive(4, idx_bound) };
		FlowSolverPara* para = GetPara();
		auto gas = GetGas();
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
					m_data_manager->SetPrimitive(iVal, idx_bound, prim_in[iVal]);
				}
			}
			else // 超声速入口
			{
				m_data_manager->SetPrimitive(0, idx_bound, para->GetInflowDensity());
				m_data_manager->SetPrimitive(1, idx_bound, para->GetInflowVelocityX());
				m_data_manager->SetPrimitive(2, idx_bound, para->GetInflowVelocityY());
				m_data_manager->SetPrimitive(3, idx_bound, para->GetInflowVelocityZ());
				m_data_manager->SetPrimitive(4, idx_bound, para->GetInflowPressure());
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
				m_data_manager->SetPrimitive(0, idx_bound, pow(c_bound * c_bound / (entropy * gamma), 1.0 / gamma1));
				m_data_manager->SetPrimitive(1, idx_bound, prim_far[1] + norm_bound[0] * (vn_bound - vn_far));
				m_data_manager->SetPrimitive(2, idx_bound, prim_far[2] + norm_bound[1] * (vn_bound - vn_far));
				m_data_manager->SetPrimitive(3, idx_bound, prim_far[3] + norm_bound[2] * (vn_bound - vn_far));
				m_data_manager->SetPrimitive(4, idx_bound, c_bound * c_bound * m_data_manager->GetPrimitive(0, idx_bound) / gamma);
			}
			else // 亚声速出口
			{
				double entropy = prim_in[4] / pow(prim_in[0], gamma);
				m_data_manager->SetPrimitive(0, idx_bound, pow(c_bound * c_bound / (entropy * gamma), 1.0 / gamma1));
				m_data_manager->SetPrimitive(1, idx_bound, prim_in[1] + norm_bound[0] * (vn_bound - vn_in));
				m_data_manager->SetPrimitive(2, idx_bound, prim_in[2] + norm_bound[1] * (vn_bound - vn_in));
				m_data_manager->SetPrimitive(3, idx_bound, prim_in[3] + norm_bound[2] * (vn_bound - vn_in));
				m_data_manager->SetPrimitive(4, idx_bound, c_bound * c_bound * m_data_manager->GetPrimitive(0, idx_bound) / gamma);
			}
		}
		for (int iVal = 0; iVal < 5; ++iVal)
			prim_bound[iVal] = m_data_manager->GetPrimitive(iVal, idx_bound);
		double cons_bound[5];
		GetGas()->Prim2Cons(prim_bound, cons_bound);
		for (int iVal = 0; iVal < 5; ++iVal)
		{
			m_data_manager->SetConservative(iVal, idx_bound, cons_bound[iVal]);
		}
	}
	void NSSolverFNFDM::CalcPrimGradBound()
	{
		int equ_num = GetPara()->GetEquNum();
		auto grid = GetGrid();
		auto node = grid->GetNode();
		BoundaryMap* bound_map = grid->GetBoundaryMap();
		auto& boundaryMap = bound_map->GetBoundaryMap();
		for (auto& boundary : boundaryMap)
		{
			auto& boundName = boundary.first;
			auto& bound = boundary.second;
			if (boundName == "hole")
				continue;
#pragma omp parallel for
			for (int iVal = 0; iVal < equ_num; ++iVal)
			{
				for (int iBound = 0; iBound < bound.size(); ++iBound)
				{
					auto& bound_index = bound[iBound].GetIdxBound();
					auto& inner_index = bound[iBound].GetIdxRef();
					for (int iDim = 0; iDim < 3; ++iDim)
					{
						m_data_manager->SetPrimitiveGrad(iVal, iDim, bound_index, 0.0);
					}
				}
			}
		}
	}

	void NSSolverFNFDM::CalcGradWLS()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		int equ_num = GetPara()->GetEquNum();
		for (int iEqu = 0; iEqu < equ_num; ++iEqu)
		{
			m_grad_wlsq->CalcGradient(grid, m_data_manager->GetPrimitive(iEqu), m_data_manager->GetPrimitiveGrad(iEqu, 0), m_data_manager->GetPrimitiveGrad(iEqu, 1), m_data_manager->GetPrimitiveGrad(iEqu, 2));
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
		auto node = grid->GetNode();
		int node_num = grid->GetTotalNodeNum();
		int equ_num = GetPara()->GetEquNum();
		double maxVal, minVal;
		double eps = 1e-6;
		double venkatCoeff = 1.0e-5;
		for (int iVal = 0; iVal < equ_num; ++iVal)
		{
#pragma omp parallel for private(maxVal, minVal, eps)
			for (int iNode = 0; iNode < node_num; ++iNode)
			{
				if (node->GetType(iNode) != NodeType::inner && node->GetType(iNode) != NodeType::hole)
					continue;
				auto currentNodeCoord = node->GetCoord(iNode);
				auto neighborNode = node->GetNeighborNode(iNode);
				int nNeighbor = node->GetNeighborNodeNum(iNode);
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
					auto neighborNodeCoord = node->GetCoord(neighborNode[iNeighbor]);
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
					m_data_manager->SetLimiter(iVal, iNode, Min(m_data_manager->GetLimiter(iVal, iNode), tempCoef));
				}
			}
		}
	}

	void NSSolverFNFDM::CalcLimiterBJ()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		int node_num = grid->GetTotalNodeNum();
		int equ_num = GetPara()->GetEquNum();
		double maxVal, minVal;
#pragma omp parallel for private(maxVal, minVal)
		for (int iVal = 0; iVal < equ_num; ++iVal)
		{
			for (int iNode = 0; iNode < node_num; ++iNode)
			{
				if (node->GetType(iNode) != NodeType::inner && node->GetType(iNode) != NodeType::hole)
					continue;
				auto currentNodeCoord = node->GetCoord(iNode);
				auto neighborNode = node->GetNeighborNode(iNode);
				auto nNeighbor = node->GetNeighborNodeNum(iNode);
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
					auto neighborNodeCoord = node->GetCoord(neighborNode[iNeighbor]);
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
					m_data_manager->SetLimiter(iVal, iNode, Min(m_data_manager->GetLimiter(iVal, iNode), tempCoef));
				}
			}
		}
	}

	void NSSolverFNFDM::CalcLimiterNone()
	{
		auto para = GetPara();
		auto grid = GetGrid();
		auto node = grid->GetNode();
		int node_num = grid->GetTotalNodeNum();
		for (int iVal = 0; iVal < para->GetEquNum(); ++iVal)
		{
#pragma omp parallel for
			for (int iNode = 0; iNode < node_num; ++iNode)
			{
				if (node->GetType(iNode) != NodeType::inner && node->GetType(iNode) != NodeType::hole)
					continue;
				m_data_manager->SetLimiter(iVal, iNode, 1.0);
			}
		}
	}

	void NSSolverFNFDM::CalcLimiterFirstOrder()
	{
		auto para = GetPara();
		auto grid = GetGrid();
		auto node = grid->GetNode();
		int node_num = grid->GetTotalNodeNum();
		for (int iVal = 0; iVal < para->GetEquNum(); ++iVal)
		{
#pragma omp parallel for
			for (int iNode = 0; iNode < node_num; ++iNode)
			{
				m_data_manager->SetLimiter(iVal, iNode, 0.0);
			}
		}
	}

	void NSSolverFNFDM::CalcLimiterBound()
	{
		int equ_num = GetPara()->GetEquNum();
		auto grid = GetGrid();
		auto boundaryMapPtr = grid->GetBoundaryMap();
		auto& boundaryMap = boundaryMapPtr->GetBoundaryMap();
		for (auto& boundary : boundaryMap)
		{
			auto& boundName = boundary.first;
			auto& bound = boundary.second;
			if (boundName == "hole")
				continue;
			for (int iVal = 0; iVal < equ_num; ++iVal)
			{
#pragma omp parallel for
				for (int iBound = 0; iBound < bound.size(); ++iBound)
				{
					auto& bound_index = bound[iBound].GetIdxBound();
					auto& inner_index = bound[iBound].GetIdxRef();
					m_data_manager->SetLimiter(iVal, bound_index, 0.0);
				}
			}
		}
	}

	void NSSolverFNFDM::CheckPrimtive()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		int node_num = grid->GetTotalNodeNum();
		int equ_num = GetPara()->GetEquNum();
		int nonphysical_node_num = 0;
#pragma omp parallel for  reduction(+ : nonphysical_node_num)
		for (int iNode = 0; iNode < node_num; ++iNode)
		{
			bool exist_nonphysical = false;
			m_data_manager->SetNonPhysical(iNode, -1);
			if (node->GetType(iNode) != NodeType::inner && node->GetType(iNode) != NodeType::hole)
				continue;
			if (m_data_manager->GetDensity(iNode) < 0 || m_data_manager->GetPressure(iNode) < 0)
			{
				exist_nonphysical = true;
			}
			if (!exist_nonphysical)
			{
				for (int iVal = 0; iVal < equ_num; ++iVal)
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
					node->GetNeighborNodeNum(iNode), m_data_manager->GetPrimitive(0, iNode), m_data_manager->GetPrimitive(1, iNode),
					m_data_manager->GetPrimitive(2, iNode), m_data_manager->GetPrimitive(3, iNode), m_data_manager->GetPrimitive(4, iNode));
				Log::info("Non-physical Node: {}, coord: {},{},{},", iNode, node->GetCoord(iNode)[0],
					node->GetCoord(iNode)[1], node->GetCoord(iNode)[2]);
			}
		}
		if (nonphysical_node_num > 0)
		{
			Log::warn("Non-physical Node Num: {}", nonphysical_node_num);
			auto para = GetPara();
			double cfl = para->GetCflNumber();
			cfl = cfl / 5.0;
			para->SetCflNumber(cfl);
			Log::warn("CFL Number is reduced to {}", cfl);
		}
		else
		{
			auto para = GetPara();
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
		auto node_topo = grid->GetNode();
		int total_node_num = grid->GetTotalNodeNum();
		int equ_num = GetPara()->GetEquNum();
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
			for (int iVal = 0; iVal < equ_num; ++iVal)
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
			GetGas()->Prim2Cons(prim, cons);
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
		auto node = grid->GetNode();
		auto para = GetPara();
		double gamma = para->GetGas()->GetGamma();
		double cfl = para->GetCflNumber();
		int inner_node_num = grid->GetInnerNodeNum();
		int* inner_node = grid->GetInnerNode();
		double min_dt = LARGE_NUMBER;
#pragma omp parallel for reduction(min:min_dt)
		for (int iNode = 0; iNode < inner_node_num; ++iNode)
		{
			int idx_node = inner_node[iNode];
			auto xi = m_node_metric->GetMetricXi(idx_node);
			auto eta = m_node_metric->GetMetricEta(idx_node);
			auto zeta = m_node_metric->GetMetricZeta(idx_node);
			auto jacobi = m_node_metric->GetJacobian(idx_node);
			double c = sqrt(gamma * m_data_manager->GetPressure(idx_node) / m_data_manager->GetDensity(idx_node));
			double norm_xi = sqrt(xi[0] * xi[0] + xi[1] * xi[1] + xi[2] * xi[2]);
			double norm_eta = sqrt(eta[0] * eta[0] + eta[1] * eta[1] + eta[2] * eta[2]);
			double norm_zeta = sqrt(zeta[0] * zeta[0] + zeta[1] * zeta[1] + zeta[2] * zeta[2]);
			double u_xi = m_data_manager->GetVelocity(0, idx_node) * xi[0] + m_data_manager->GetVelocity(1, idx_node) * xi[1] + m_data_manager->GetVelocity(2, idx_node) * xi[2];
			double u_eta = m_data_manager->GetVelocity(0, idx_node) * eta[0] + m_data_manager->GetVelocity(1, idx_node) * eta[1] + m_data_manager->GetVelocity(2, idx_node) * eta[2];
			double u_zeta = m_data_manager->GetVelocity(0, idx_node) * zeta[0] + m_data_manager->GetVelocity(1, idx_node) * zeta[1] + m_data_manager->GetVelocity(2, idx_node) * zeta[2];
			double lamda = abs(u_xi) + abs(u_eta) + abs(u_zeta) + c * (norm_xi + norm_eta + norm_zeta);
			m_data_manager->SetTimeStep(idx_node, cfl / lamda);
			if (min_dt > m_data_manager->GetTimeStep(idx_node))
			{
				min_dt = m_data_manager->GetTimeStep(idx_node);
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
		int equ_num = GetPara()->GetEquNum();
		auto grid = GetGrid();
		auto node = grid->GetNode();
		int node_num = grid->GetTotalNodeNum();
		for (int iVar = 0; iVar < equ_num; ++iVar)
		{
#pragma omp parallel for
			for (int iNode = 0; iNode < node_num; ++iNode)
			{
				m_data_manager->SetResidual(iVar, iNode, 0.0);
			}
		}
	}

	void NSSolverFNFDM::ConvectiveResidual()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		int equ_num = GetPara()->GetEquNum();
		int inner_node_num = grid->GetInnerNodeNum();
		int* inner_node = grid->GetInnerNode();
		RiemannSolverPara riemann_para[6];
		for (int i = 0; i < 6; ++i)
		{
			riemann_para[i].gamma_left = riemann_para[i].gamma_right = GetPara()->GetGas()->GetGamma();
		}
		bool exist_negative = false;
#pragma omp parallel for private(riemann_para,exist_negative)
		for (int iNode = 0; iNode < inner_node_num; ++iNode)
		{
			int idx = inner_node[iNode];
			exist_negative = false;
			double jacobi = m_node_metric->GetJacobian(idx);
			auto neighbor = node->GetNeighborNode(idx);
			// i direction
			riemann_para[0].norm(0) = m_node_metric->GetMetricXi(idx)[0];
			riemann_para[0].norm(1) = m_node_metric->GetMetricXi(idx)[1];
			riemann_para[0].norm(2) = m_node_metric->GetMetricXi(idx)[2];
			riemann_para[0].nt = m_node_metric->GetMetricXi(idx)[3];
			riemann_para[1].norm(0) = m_node_metric->GetMetricXi(idx)[0];
			riemann_para[1].norm(1) = m_node_metric->GetMetricXi(idx)[1];
			riemann_para[1].norm(2) = m_node_metric->GetMetricXi(idx)[2];
			riemann_para[1].nt = m_node_metric->GetMetricXi(idx)[3];
			MidPointReconstruct2ndOrder(idx, neighbor[1], &riemann_para[0].prim_left(0), &riemann_para[0].prim_right(0));
			MidPointReconstruct2ndOrder(neighbor[0], idx, &riemann_para[1].prim_left(0), &riemann_para[1].prim_right(0));
			// j direction
			riemann_para[2].norm(0) = m_node_metric->GetMetricEta(idx)[0];
			riemann_para[2].norm(1) = m_node_metric->GetMetricEta(idx)[1];
			riemann_para[2].norm(2) = m_node_metric->GetMetricEta(idx)[2];
			riemann_para[2].nt = m_node_metric->GetMetricEta(idx)[3];
			riemann_para[3].norm(0) = m_node_metric->GetMetricEta(idx)[0];
			riemann_para[3].norm(1) = m_node_metric->GetMetricEta(idx)[1];
			riemann_para[3].norm(2) = m_node_metric->GetMetricEta(idx)[2];
			riemann_para[3].nt = m_node_metric->GetMetricEta(idx)[3];
			MidPointReconstruct2ndOrder(idx, neighbor[3], &riemann_para[2].prim_left(0), &riemann_para[2].prim_right(0));
			MidPointReconstruct2ndOrder(neighbor[2], idx, &riemann_para[3].prim_left(0), &riemann_para[3].prim_right(0));

			// k direction
			riemann_para[4].norm(0) = m_node_metric->GetMetricZeta(idx)[0];
			riemann_para[4].norm(1) = m_node_metric->GetMetricZeta(idx)[1];
			riemann_para[4].norm(2) = m_node_metric->GetMetricZeta(idx)[2];
			riemann_para[4].nt = m_node_metric->GetMetricZeta(idx)[3];
			riemann_para[5].norm(0) = m_node_metric->GetMetricZeta(idx)[0];
			riemann_para[5].norm(1) = m_node_metric->GetMetricZeta(idx)[1];
			riemann_para[5].norm(2) = m_node_metric->GetMetricZeta(idx)[2];
			riemann_para[5].nt = m_node_metric->GetMetricZeta(idx)[3];
			MidPointReconstruct2ndOrder(idx, neighbor[5], &riemann_para[4].prim_left(0), &riemann_para[4].prim_right(0));
			MidPointReconstruct2ndOrder(neighbor[4], idx, &riemann_para[5].prim_left(0), &riemann_para[5].prim_right(0));
			// check negative density and pressure
			for (int i = 0; i < 6; ++i)
			{
				if (riemann_para[i].prim_left(0) < 0 || riemann_para[i].prim_right(0) < 0 || riemann_para[i].prim_left(4) < 0 || riemann_para[i].prim_right(4) < 0)
				{
					exist_negative = true;
					break;
				}
			}
			if (exist_negative || m_data_manager->GetNonPhysical(idx) > 0)
			{
				MidPointReconstruct1stOrder(idx, neighbor[1], &riemann_para[0].prim_left(0), &riemann_para[0].prim_right(0));
				MidPointReconstruct1stOrder(neighbor[0], idx, &riemann_para[1].prim_left(0), &riemann_para[1].prim_right(0));
				MidPointReconstruct1stOrder(idx, neighbor[3], &riemann_para[2].prim_left(0), &riemann_para[2].prim_right(0));
				MidPointReconstruct1stOrder(neighbor[2], idx, &riemann_para[3].prim_left(0), &riemann_para[3].prim_right(0));
				MidPointReconstruct1stOrder(idx, neighbor[5], &riemann_para[4].prim_left(0), &riemann_para[4].prim_right(0));
				MidPointReconstruct1stOrder(neighbor[4], idx, &riemann_para[5].prim_left(0), &riemann_para[5].prim_right(0));
			}
			// calculate flux
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

	void NSSolverFNFDM::ViscousResidual()
	{
		return;
		CalcViscousFlux();
		CalcViscousFluxGrad();
		auto grid = GetGrid();
		auto node = grid->GetNode();
		int node_num = grid->GetTotalNodeNum();
		int inner_node_num = grid->GetInnerNodeNum();
		int* inner_node = grid->GetInnerNode();
		auto para = GetPara();
		double res_vis;
		for (int iNode = 0; iNode < inner_node_num; ++iNode)
		{
			int idx = inner_node[iNode];
			for (int iVar = 0; iVar < para->GetEquNum(); ++iVar)
			{
				res_vis = m_data_manager->GetViscousFluxGrad(iVar, 0, 0, idx) + m_data_manager->GetViscousFluxGrad(iVar, 1, 1, idx) + m_data_manager->GetViscousFluxGrad(iVar, 2, 2, idx);
				res_vis /= para->GetDimensionless().GetRe();
				m_data_manager->SetResidual(iVar, idx, m_data_manager->GetResidual(iVar, idx) - res_vis);
			}
		}
	}
	void NSSolverFNFDM::CalcViscousFlux()
	{
		auto para = GetPara();
		auto gas = GetGas();
		auto grid = GetGrid();
		auto node = grid->GetNode();
		int node_num = grid->GetTotalNodeNum();
		double* viscous_flux_x, * viscous_flux_y, * viscous_flux_z;
		double tau_xx, tau_yy, tau_zz, tau_xy, tau_xz, tau_yz;
		double grad_rho[3], grad_u[3], grad_v[3], grad_w[3], grad_p[3];
		double vel[3];
		double tempeture;
		double vis_coef;
		double therm_coef;
		double lamda;
		for (int iNode = 0; iNode < node_num; ++iNode)
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
		int equ_num = GetPara()->GetEquNum();
		auto grid = GetGrid();
		for (int iEqu = 0; iEqu < equ_num; ++iEqu)
		{
			m_grad_wlsq->CalcGradient(grid, m_data_manager->GetViscousFlux(iEqu, 0), m_data_manager->GetViscousFluxGrad(iEqu, 0, 0), m_data_manager->GetViscousFluxGrad(iEqu, 0, 1), m_data_manager->GetViscousFluxGrad(iEqu, 0, 2));
			m_grad_wlsq->CalcGradient(grid, m_data_manager->GetViscousFlux(iEqu, 1), m_data_manager->GetViscousFluxGrad(iEqu, 1, 0), m_data_manager->GetViscousFluxGrad(iEqu, 1, 1), m_data_manager->GetViscousFluxGrad(iEqu, 1, 2));
			m_grad_wlsq->CalcGradient(grid, m_data_manager->GetViscousFlux(iEqu, 2), m_data_manager->GetViscousFluxGrad(iEqu, 2, 0), m_data_manager->GetViscousFluxGrad(iEqu, 2, 1), m_data_manager->GetViscousFluxGrad(iEqu, 2, 2));
		}
	}
	void NSSolverFNFDM::SourceTermResidual()
	{
	}

	void NSSolverFNFDM::CalcForce()
	{
		return;
		auto grid = GetGrid();
		auto face = grid->GetFace();
		auto para = GetPara();
		auto dimensionless = para->GetDimensionless();
		int face_num = face->GetFaceNum();
		double face_pressure;
		double force[3] = { 0 };
		for (int iFace = 0; iFace < face_num; ++iFace)
		{
			face_pressure = 0;
			int* face2node = face->GetFace2Node(iFace);
			for (int iNode = 0; iNode < face->GetFaceNodeNum(iFace); ++iNode)
			{
				face_pressure += m_data_manager->GetPressure(face2node[iNode]);
			}
			face_pressure /= face->GetFaceNodeNum(iFace);
			face_pressure -= dimensionless.GetPressureDL(0);
			force[0] -= face_pressure * face->GetArea(iFace) * face->GetNormal(iFace)[0];
			force[1] -= face_pressure * face->GetArea(iFace) * face->GetNormal(iFace)[1];
			force[2] -= face_pressure * face->GetArea(iFace) * face->GetNormal(iFace)[2];
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