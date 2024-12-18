#include "NSSolverStruct.h"
#include "FlowSolverStructPara.h"
#include "MathBasic.h"
#include"Log.h"
#include <omp.h>
#include "Metric.h"

namespace zaran
{
	NSSolverStruct::NSSolverStruct(index_type index, string name, shared_ptr<FlowSolverParamStruct> para, shared_ptr < GridStruct> grid, shared_ptr < DataManagerNSStruct>data_manager)
		: NSSolver(index, name, para, grid, data_manager)
	{
		int ni, nj, nk;
		ni = grid->GetNi();
		nj = grid->GetNj();
		nk = grid->GetNk();
		m_idx_proxy = make_shared<IdProxyStruct>(ni, nj, nk);
		int node_num = ni * nj * nk;
		m_node_metrics = make_shared<Metric>(node_num);
		m_metrics_mid_i = new Metric(node_num);
		m_metrics_mid_j = new Metric(node_num);
		m_metrics_mid_k = new Metric(node_num);
		m_grid = grid;
		m_data_manager = data_manager;
		m_para = para;
	}
	NSSolverStruct::~NSSolverStruct()
	{

		if (m_metrics_mid_i)
		{
			delete m_metrics_mid_i;
			m_metrics_mid_i = nullptr;
		}
		if (m_metrics_mid_j)
		{
			delete m_metrics_mid_j;
			m_metrics_mid_j = nullptr;
		}
		if (m_metrics_mid_k)
		{
			delete m_metrics_mid_k;
			m_metrics_mid_k = nullptr;
		}
	}
	DataManagerNSStruct* NSSolverStruct::GetDataManager()
	{
		return m_data_manager.get();
	}
	Metric* NSSolverStruct::GetNodeMetrics()
	{
		return m_node_metrics.get();
	}
	IdProxyStruct* NSSolverStruct::GetIdxProxy()
	{
		return m_idx_proxy.get();
	}
	FlowSolverParamStruct* NSSolverStruct::GetPara()
	{
		return m_para.get();
	}
	void NSSolverStruct::InitField()
	{
		auto para = GetPara();
		auto& init_type = para->GetInitFieldType();
		if (init_type == InitFieldType::Vortex)
		{
			InitFieldVortex();
			Prim2Cons();
		}
		else
		{
			NSSolver::InitField();
		}

	}
	void NSSolverStruct::InitFieldFarfield()
	{
		// InitFieldVortex();
		// return;
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		auto ni = grid->GetNi();
		auto nj = grid->GetNj();
		auto nk = grid->GetNk();
		auto para = GetPara();
		double prim_far[5];
		prim_far[0] = para->GetInflowDensity();
		prim_far[1] = para->GetInflowVelocityX();
		prim_far[2] = para->GetInflowVelocityY();
		prim_far[3] = para->GetInflowVelocityZ();
		prim_far[4] = para->GetInflowPressure();
		// prim_far[1] = 3.0;
		// prim_far[2] = 0.0;
		// prim_far[3] = 0.0;
		// prim_far[4] = 1.4;
		// prim_far[0] = 3.87330;
		// prim_far[1] = 0;
		// prim_far[2] = 0;
		// prim_far[3] = 0.0;
		// prim_far[4] = 4.82466;
		for (int k = 0; k < nk; ++k)
		{
			for (int j = 0; j < nj; ++j)
			{
				for (int i = 0; i < ni; ++i)
				{
					m_idx_proxy->SetIdx(i, j, k);
					int idx = m_idx_proxy->GetIdx();
					//double x = grid->GetNode()->GetCoord(i, j, k)[0];
					//double y = grid->GetNode()->GetCoord(i, j, k)[1];
					//prim_far[0] = pow(x, 3.5) + pow(y, 1.1);
					//prim_far[1] = pow(x, 2.5) + pow(y, 0.1);
					//prim_far[2] = pow(x, 1.0) + pow(y, 0.1);
					//prim_far[3] = 0.0;
					//prim_far[4] = pow(x, 0.2) + pow(y, 1.2);

					data_manager->SetPrim(idx, prim_far);
				}
			}
		}
	}
	void NSSolverStruct::InitFieldFarFieldZeroVel()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto data_manager = GetDataManager();
		auto ni = grid->GetNi();
		auto nj = grid->GetNj();
		auto nk = grid->GetNk();
		auto para = GetPara();
		double prim_far[5];
		prim_far[0] = para->GetInflowDensity();
		prim_far[1] = 0.0;
		prim_far[2] = 0.0;
		prim_far[3] = 0.0;
		prim_far[4] = para->GetInflowPressure();
		// prim_far[1] = 0.0;
		// prim_far[2] = 0.0;
		// prim_far[3] = 0.0;
		// prim_far[4] = 1.4;
		//double prim_left[5] = { 8.0,8.25 * cos(30.0 / 180.0 * PI),-8.25 * sin(30.0 / 180.0 * PI),0.0,116.5 };
		//double prim_left[5] = { 6.4,3.125,0,0,18.5 };
		//double prim_right[5] = { 1.4,0.0,0.0,0.0,1.0 };
		for (int k = 0; k < nk; ++k)
		{
			for (int j = 0; j < nj; ++j)
			{
				for (int i = 0; i < ni; ++i)
				{
					m_idx_proxy->SetIdx(i, j, k);
					int idx = m_idx_proxy->GetIdx();
					// double y = grid->GetNode()->GetCoord(i, j, k)[1];
					// prim_far[0] = 1.0 + 0.01 * y;
					data_manager->SetPrim(idx, prim_far);
					//auto coord = node->GetCoord(idx);
					//if (coord[0] <= 1.0 / 6.0 + coord[1] / tan(60.0 / 180.0 * PI))
					//	//if (coord[0] <= 0.1)
					//{
					//	data_manager->SetPrim(idx, prim_left);
					//}
					//else
					//{
					//	data_manager->SetPrim(idx, prim_right);
					//}
				}
			}
		}
	}
	void NSSolverStruct::InitFieldBackup()
	{
	}
	void NSSolverStruct::InitFieldVortex()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto data_manager = GetDataManager();
		auto ni = grid->GetNi();
		auto nj = grid->GetNj();
		auto nk = grid->GetNk();
		auto para = GetPara();
		double beta = 5.0;
		double x, y, z;
		double r2;
		double prim[5];
		double gamma = 1.4;
		for (int k = 0; k < nk; ++k)
		{
			for (int j = 0; j < nj; ++j)
			{
				for (int i = 0; i < ni; ++i)
				{
					m_idx_proxy->SetIdx(i, j, k);
					int idx = m_idx_proxy->GetIdx();
					x = node->GetCoord(i, j, k)[0];
					y = node->GetCoord(i, j, k)[1];
					z = node->GetCoord(i, j, k)[2];
					r2 = x * x + y * y;
					prim[0] = pow(1.0 - (gamma - 1.0) * beta * beta * exp(1.0 - r2) / (8.0 * gamma * PI * PI),
						1.0 / (gamma - 1.0));
					prim[4] = pow(prim[0], gamma);
					prim[1] = beta * exp(0.5 * (1.0 - r2)) / (2.0 * PI) * (-y);
					prim[2] = beta * exp(0.5 * (1.0 - r2)) / (2.0 * PI) * (x);
					prim[3] = 0.0;
					data_manager->SetPrim(idx, prim);
				}
			}
		}
	}
	void NSSolverStruct::InitFieldExplosion()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto data_manager = GetDataManager();
		auto ni = grid->GetNi();
		auto nj = grid->GetNj();
		auto nk = grid->GetNk();
		auto para = GetPara();
		double prim_far[5];
		prim_far[0] = para->GetInflowDensity();
		prim_far[1] = 0.0;
		prim_far[2] = 0.0;
		prim_far[3] = 0.0;
		prim_far[4] = para->GetInflowPressure();
		double explosion_pressure = GlobalData::GetDouble("explosion_pressure");
		double prim_explosion[5];
		prim_explosion[0] = prim_far[0];
		prim_explosion[1] = 0.0;
		prim_explosion[2] = 0.0;
		prim_explosion[3] = 0.0;
		prim_explosion[4] = explosion_pressure * prim_far[4];
		double explosion_center[3];
		explosion_center[0] = GlobalData::GetDouble("explosion_center_x");
		explosion_center[1] = GlobalData::GetDouble("explosion_center_y");
		explosion_center[2] = GlobalData::GetDouble("explosion_center_z");
		double explosion_radius = GlobalData::GetDouble("explosion_radius");
		double x, y, z;
		for (int k = 0; k < nk; ++k)
		{
			for (int j = 0; j < nj; ++j)
			{
				for (int i = 0; i < ni; ++i)
				{
					auto idx = m_idx_proxy->GetIdx();
					auto coord = node->GetCoord(idx);
					double dist = sqrt(pow(coord[0] - explosion_center[0], 2)
						+ pow(coord[1] - explosion_center[1], 2)
						+ pow(coord[2] - explosion_center[2], 2));
					if (dist < explosion_radius)
					{
						data_manager->SetPrim(idx, prim_explosion);
					}
					else
					{
						data_manager->SetPrim(idx, prim_far);
					}
				}
			}
		}
	}
	void NSSolverStruct::CalcMetricsS0()
	{
		auto para = GetPara();
		auto metric_diff_scheme = para->GetMetricDifferenceScheme();
		if (metric_diff_scheme == MetricDifferenceScheme::SecondOrder)
		{
			CalcMetricsS0_2nd();
		}
		else if (metric_diff_scheme == MetricDifferenceScheme::SixthOrder)
		{
			CalcMetricsS0_6th();
		}
		else
		{
			Log::warn("FluxDifferenceScheme is not defined, use SecondOrder as default");
			CalcMetricsS0_2nd();
		}
	}

	void NSSolverStruct::CalcMetricsS0_2nd()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto coef = GetNodeMetrics();
		Metric* coef_mid[3] = { GetMidMetricsI(), GetMidMetricsJ(), GetMidMetricsK() };
		auto idx_proxy = GetIdxProxy();
		auto Idx = [&](int i, int j, int k) { return idx_proxy->GetIdx(i, j, k); };
		int ni = grid->GetNi();
		int nj = grid->GetNj();
		int nk = grid->GetNk();
		// i+1/2,j+1/2,k+1/2处的坐标
		dynamic_array<dynamic_array<double>> coord_i(ni * nj * nk, dynamic_array<double>(3)),
			coord_j(ni * nj * nk, dynamic_array<double>(3)), coord_k(ni * nj * nk, dynamic_array<double>(3));
		// 计算之前，先把度量系数赋值为0
		for (int i = 0; i < ni; ++i)
		{
			for (int j = 0; j < nj; ++j)
			{
				for (int k = 0; k < nk; ++k)
				{
					for (int iDim = 0; iDim < 4; ++iDim)
					{
						coef->GetX(Idx(i, j, k))[iDim] = 0.0;
						coef->GetY(Idx(i, j, k))[iDim] = 0.0;
						coef->GetZ(Idx(i, j, k))[iDim] = 0.0;
						coef->GetXi(Idx(i, j, k))[iDim] = 0.0;
						coef->GetEta(Idx(i, j, k))[iDim] = 0.0;
						coef->GetZeta(Idx(i, j, k))[iDim] = 0.0;
						for (int jDim = 0; jDim < 3; ++jDim)
						{
							coef_mid[jDim]->GetX(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetY(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetZ(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetXi(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetEta(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetZeta(Idx(i, j, k))[iDim] = 0.0;
						}
					}
				}
			}
		}
		// 第一步：计算半点坐标
		for (int k = 0; k < nk - 1; ++k)
		{
			for (int j = 0; j < nj - 1; ++j)
			{
				for (int i = 0; i < ni - 1; ++i)
				{
					for (int iDim = 0; iDim < 3; ++iDim)
					{
						int idx = idx_proxy->GetIdx(i, j, k);
						coord_i[idx][iDim] = 0.5 * (node->GetCoord(i, j, k)[iDim] + node->GetCoord(i + 1, j, k)[iDim]);
						coord_j[idx][iDim] = 0.5 * (node->GetCoord(i, j, k)[iDim] + node->GetCoord(i, j + 1, k)[iDim]);
						coord_k[idx][iDim] = 0.5 * (node->GetCoord(i, j, k)[iDim] + node->GetCoord(i, j, k + 1)[iDim]);
					}
				}
			}
		}
		// 第二步：根据半点坐标计算整点度量系数
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					int idx = idx_proxy->GetIdx(i, j, k);
					auto coef_x = coef->GetX(idx);
					auto coef_y = coef->GetY(idx);
					auto coef_z = coef->GetZ(idx);
					coef_x[0] = coord_i[Idx(i, j, k)][0] - coord_i[Idx(i - 1, j, k)][0];
					coef_x[1] = coord_j[Idx(i, j, k)][0] - coord_j[Idx(i, j - 1, k)][0];
					coef_x[2] = coord_k[Idx(i, j, k)][0] - coord_k[Idx(i, j, k - 1)][0];
					coef_y[0] = coord_i[Idx(i, j, k)][1] - coord_i[Idx(i - 1, j, k)][1];
					coef_y[1] = coord_j[Idx(i, j, k)][1] - coord_j[Idx(i, j - 1, k)][1];
					coef_y[2] = coord_k[Idx(i, j, k)][1] - coord_k[Idx(i, j, k - 1)][1];
					coef_z[0] = coord_i[Idx(i, j, k)][2] - coord_i[Idx(i - 1, j, k)][2];
					coef_z[1] = coord_j[Idx(i, j, k)][2] - coord_j[Idx(i, j - 1, k)][2];
					coef_z[2] = coord_k[Idx(i, j, k)][2] - coord_k[Idx(i, j, k - 1)][2];
					coef_x[3] = 0.0;
					coef_y[3] = 0.0;
					if (grid->GetDim() == TWO_DIM)
					{
						coef_x[2] = 0.0;
						coef_y[2] = 0.0;
						coef_z[0] = 0.0;
						coef_z[1] = 0.0;
						coef_z[2] = 1.0;
					}
					// 计算xi,eta,zeta
					auto coef_xi = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, k));
					coef_xi[0] = coef_y[1] * coef_z[2] - coef_z[1] * coef_y[2];
					coef_xi[1] = coef_z[1] * coef_x[2] - coef_x[1] * coef_z[2];
					coef_xi[2] = coef_x[1] * coef_y[2] - coef_y[1] * coef_x[2];
					coef_xi[3] = 0.0;
					auto coef_eta = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, k));
					coef_eta[0] = coef_y[2] * coef_z[0] - coef_z[2] * coef_y[0];
					coef_eta[1] = coef_z[2] * coef_x[0] - coef_x[2] * coef_z[0];
					coef_eta[2] = coef_x[2] * coef_y[0] - coef_y[2] * coef_x[0];
					coef_eta[3] = 0.0;
					auto coef_zeta = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, k));
					coef_zeta[0] = coef_y[0] * coef_z[1] - coef_z[0] * coef_y[1];
					coef_zeta[1] = coef_z[0] * coef_x[1] - coef_x[0] * coef_z[1];
					coef_zeta[2] = coef_x[0] * coef_y[1] - coef_y[0] * coef_x[1];
					coef_zeta[3] = 0.0;
				}
			}
		}
	}

	void NSSolverStruct::CalcMetricsS0_6th()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto coef = GetNodeMetrics();
		Metric* coef_mid[3] = { GetMidMetricsI(), GetMidMetricsJ(), GetMidMetricsK() };
		auto idx_proxy = GetIdxProxy();
		auto Idx = [&](int i, int j, int k) { return idx_proxy->GetIdx(i, j, k); };
		int ni = grid->GetNi();
		int nj = grid->GetNj();
		int nk = grid->GetNk();
		double inter_temp[6];
		// 用于度量系数计算的临时变量
		dynamic_array<dynamic_array<double>> coord_i(ni * nj * nk, dynamic_array<double>(3)),
			coord_j(ni * nj * nk, dynamic_array<double>(3)), coord_k(ni * nj * nk, dynamic_array<double>(3));
		// 计算之前，先把度量系数赋值为0
		for (int i = 0; i < ni; ++i)
		{
			for (int j = 0; j < nj; ++j)
			{
				for (int k = 0; k < nk; ++k)
				{
					for (int iDim = 0; iDim < 4; ++iDim)
					{
						coef->GetX(Idx(i, j, k))[iDim] = 0.0;
						coef->GetY(Idx(i, j, k))[iDim] = 0.0;
						coef->GetZ(Idx(i, j, k))[iDim] = 0.0;
						coef->GetXi(Idx(i, j, k))[iDim] = 0.0;
						coef->GetEta(Idx(i, j, k))[iDim] = 0.0;
						coef->GetZeta(Idx(i, j, k))[iDim] = 0.0;
						for (int jDim = 0; jDim < 3; ++jDim)
						{
							coef_mid[jDim]->GetX(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetY(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetZ(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetXi(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetEta(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetZeta(Idx(i, j, k))[iDim] = 0.0;
						}
					}
				}
			}
		}
		// 第一步：计算半点坐标
		// i+1/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 3; i < ni - 4; ++i)
				{
					for (int iDim = 0; iDim < 3; ++iDim)
					{
						for (int iTemp = 0; iTemp < 6; iTemp++)
						{
							inter_temp[iTemp] = node->GetCoord(i + iTemp - 2, j, k)[iDim];
						}
						coord_i[Idx(i, j, k)][iDim] = MidNodeInter6th(inter_temp);
					}
				}
			}
		}
		// i=1/2,3/2,5/2;ni-3/2,ni-5/2,ni-7/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int iDim = 0; iDim < 3; ++iDim)
				{
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = node->GetCoord(1 + iTemp, j, k)[iDim];
					}
					coord_i[Idx(0, j, k)][iDim] = MidNodeInter4thRight1(inter_temp);
					coord_i[Idx(1, j, k)][iDim] = MidNodeInter4thRight2(inter_temp);
					coord_i[Idx(2, j, k)][iDim] = MidNodeInter4th(inter_temp);
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = node->GetCoord(ni - 5 + iTemp, j, k)[iDim];
					}
					coord_i[Idx(ni - 2, j, k)][iDim] = MidNodeInter4thLeft1(inter_temp);
					coord_i[Idx(ni - 3, j, k)][iDim] = MidNodeInter4thLeft2(inter_temp);
					coord_i[Idx(ni - 4, j, k)][iDim] = MidNodeInter4th(inter_temp);
				}
			}
		}
		// j+1/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 3; j < nj - 4; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					for (int iDim = 0; iDim < 3; ++iDim)
					{
						for (int iTemp = 0; iTemp < 6; iTemp++)
						{
							inter_temp[iTemp] = node->GetCoord(i, j + iTemp - 2, k)[iDim];
						}
						coord_j[Idx(i, j, k)][iDim] = MidNodeInter6th(inter_temp);
					}
				}
			}
		}
		// j=1/2,3/2,5/2;nj-3/2,nj-5/2,nj-7/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				for (int iDim = 0; iDim < 3; ++iDim)
				{
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = node->GetCoord(i, 1 + iTemp, k)[iDim];
					}
					coord_j[Idx(i, 0, k)][iDim] = MidNodeInter4thRight1(inter_temp);
					coord_j[Idx(i, 1, k)][iDim] = MidNodeInter4thRight2(inter_temp);
					coord_j[Idx(i, 2, k)][iDim] = MidNodeInter4th(inter_temp);
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = node->GetCoord(i, nj - 5 + iTemp, k)[iDim];
					}
					coord_j[Idx(i, nj - 2, k)][iDim] = MidNodeInter4thLeft1(inter_temp);
					coord_j[Idx(i, nj - 3, k)][iDim] = MidNodeInter4thLeft2(inter_temp);
					coord_j[Idx(i, nj - 4, k)][iDim] = MidNodeInter4th(inter_temp);
				}
			}
		}
		// k+1/2
		for (int k = 3; k < nk - 4; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					for (int iDim = 0; iDim < 3; ++iDim)
					{
						for (int iTemp = 0; iTemp < 6; iTemp++)
						{
							inter_temp[iTemp] = node->GetCoord(i, j, k + iTemp - 2)[iDim];
						}
						coord_k[Idx(i, j, k)][iDim] = MidNodeInter6th(inter_temp);
					}
				}
			}
		}
		// k=1/2,3/2,5/2;nk-3/2,nk-5/2,nk-7/2
		for (int j = 1; j < nj - 1; ++j)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				for (int iDim = 0; iDim < 3; ++iDim)
				{
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = node->GetCoord(i, j, 1 + iTemp)[iDim];
					}
					coord_k[Idx(i, j, 0)][iDim] = MidNodeInter4thRight1(inter_temp);
					coord_k[Idx(i, j, 1)][iDim] = MidNodeInter4thRight2(inter_temp);
					coord_k[Idx(i, j, 2)][iDim] = MidNodeInter4th(inter_temp);
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = node->GetCoord(i, j, nk - 5 + iTemp)[iDim];
					}
					coord_k[Idx(i, j, nk - 2)][iDim] = MidNodeInter4thLeft1(inter_temp);
					coord_k[Idx(i, j, nk - 3)][iDim] = MidNodeInter4thLeft2(inter_temp);
					coord_k[Idx(i, j, nk - 4)][iDim] = MidNodeInter4th(inter_temp);
				}
			}
		}
		// 第二步：根据半点坐标计算整点度量系数
		double diff_temp[6];
		// i direction
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 3; i < ni - 3; ++i)
				{
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_i[Idx(i + iTemp - 3, j, k)][0];
					}
					auto coef_x = coef->GetX(Idx(i, j, k));
					coef_x[0] = NodeDifferece6th(diff_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_i[Idx(i + iTemp - 3, j, k)][1];
					}
					auto coef_y = coef->GetY(Idx(i, j, k));
					coef_y[0] = NodeDifferece6th(diff_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_i[Idx(i + iTemp - 3, j, k)][2];
					}
					auto coef_z = coef->GetZ(Idx(i, j, k));
					coef_z[0] = NodeDifferece6th(diff_temp);
					if (grid->GetDim() == TWO_DIM)
					{
						coef_z[0] = 0.0;
					}
				}
			}
		}
		// j direction
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 3; j < nj - 3; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					auto coef_x = coef->GetX(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_j[Idx(i, j + iTemp - 3, k)][0];
					}
					coef_x[1] = NodeDifferece6th(diff_temp);
					auto coef_y = coef->GetY(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_j[Idx(i, j + iTemp - 3, k)][1];
					}
					coef_y[1] = NodeDifferece6th(diff_temp);
					auto coef_z = coef->GetZ(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_j[Idx(i, j + iTemp - 3, k)][2];
					}
					coef_z[1] = NodeDifferece6th(diff_temp);
					if (grid->GetDim() == TWO_DIM)
					{
						coef_z[1] = 0.0;
					}
				}
			}
		}
		// k direction
		for (int k = 3; k < nk - 3; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					auto coef_x = coef->GetX(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_k[Idx(i, j, k + iTemp - 3)][0];
					}
					coef_x[2] = NodeDifferece6th(diff_temp);
					auto coef_y = coef->GetY(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_k[Idx(i, j, k + iTemp - 3)][1];
					}
					coef_y[2] = NodeDifferece6th(diff_temp);
					auto coef_z = coef->GetZ(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_k[Idx(i, j, k + iTemp - 3)][2];
					}
					coef_z[2] = NodeDifferece6th(diff_temp);
					if (grid->GetDim() == TWO_DIM)
					{
						coef_x[2] = 0.0;
						coef_y[2] = 0.0;
						coef_z[2] = 1.0;
					}
				}
			}
		}
		// 计算边界整点的度量系数
		// i=1,2;ni-3,ni-2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_i[Idx(0 + iTemp, j, k)][0];
				}
				coef->GetX(Idx(1, j, k))[0] = NodeDifferece4thRight(diff_temp);
				coef->GetX(Idx(2, j, k))[0] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_i[Idx(0 + iTemp, j, k)][1];
				}
				coef->GetY(Idx(1, j, k))[0] = NodeDifferece4thRight(diff_temp);
				coef->GetY(Idx(2, j, k))[0] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_i[Idx(0 + iTemp, j, k)][2];
				}
				coef->GetZ(Idx(1, j, k))[0] = NodeDifferece4thRight(diff_temp);
				coef->GetZ(Idx(2, j, k))[0] = NodeDifferece4th(diff_temp);

				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_i[Idx(ni - 2 - iTemp, j, k)][0];
				}
				coef->GetX(Idx(ni - 2, j, k))[0] = NodeDifferece4thLeft(diff_temp);
				coef->GetX(Idx(ni - 3, j, k))[0] = -NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_i[Idx(ni - 2 - iTemp, j, k)][1];
				}
				coef->GetY(Idx(ni - 2, j, k))[0] = NodeDifferece4thLeft(diff_temp);
				coef->GetY(Idx(ni - 3, j, k))[0] = -NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_i[Idx(ni - 2 - iTemp, j, k)][2];
				}
				coef->GetZ(Idx(ni - 2, j, k))[0] = NodeDifferece4thLeft(diff_temp);
				coef->GetZ(Idx(ni - 3, j, k))[0] = -NodeDifferece4th(diff_temp);
			}
		}
		// j=1,2;nj-3,nj-2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_j[Idx(i, 0 + iTemp, k)][0];
				}
				coef->GetX(Idx(i, 1, k))[1] = NodeDifferece4thRight(diff_temp);
				coef->GetX(Idx(i, 2, k))[1] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_j[Idx(i, 0 + iTemp, k)][1];
				}
				coef->GetY(Idx(i, 1, k))[1] = NodeDifferece4thRight(diff_temp);
				coef->GetY(Idx(i, 2, k))[1] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_j[Idx(i, 0 + iTemp, k)][2];
				}
				coef->GetZ(Idx(i, 1, k))[1] = NodeDifferece4thRight(diff_temp);
				coef->GetZ(Idx(i, 2, k))[1] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_j[Idx(i, nj - 2 - iTemp, k)][0];
				}
				coef->GetX(Idx(i, nj - 2, k))[1] = NodeDifferece4thLeft(diff_temp);
				coef->GetX(Idx(i, nj - 3, k))[1] = -NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_j[Idx(i, nj - 2 - iTemp, k)][1];
				}
				coef->GetY(Idx(i, nj - 2, k))[1] = NodeDifferece4thLeft(diff_temp);
				coef->GetY(Idx(i, nj - 3, k))[1] = -NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_j[Idx(i, nj - 2 - iTemp, k)][2];
				}
				coef->GetZ(Idx(i, nj - 2, k))[1] = NodeDifferece4thLeft(diff_temp);
				coef->GetZ(Idx(i, nj - 3, k))[1] = -NodeDifferece4th(diff_temp);
			}
		}
		// k=1,2;nk-3,nk-2
		for (int j = 1; j < nj - 1; ++j)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_k[Idx(i, j, 0 + iTemp)][0];
				}
				coef->GetX(Idx(i, j, 1))[2] = NodeDifferece4thRight(diff_temp);
				coef->GetX(Idx(i, j, 2))[2] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_k[Idx(i, j, 0 + iTemp)][1];
				}
				coef->GetY(Idx(i, j, 1))[2] = NodeDifferece4thRight(diff_temp);
				coef->GetY(Idx(i, j, 2))[2] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_k[Idx(i, j, 0 + iTemp)][2];
				}
				coef->GetZ(Idx(i, j, 1))[2] = NodeDifferece4thRight(diff_temp);
				coef->GetZ(Idx(i, j, 2))[2] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_k[Idx(i, j, nk - 2 - iTemp)][0];
				}
				coef->GetX(Idx(i, j, nk - 2))[2] = NodeDifferece4thLeft(diff_temp);
				coef->GetX(Idx(i, j, nk - 3))[2] = -NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_k[Idx(i, j, nk - 2 - iTemp)][1];
				}
				coef->GetY(Idx(i, j, nk - 2))[2] = NodeDifferece4thLeft(diff_temp);
				coef->GetY(Idx(i, j, nk - 3))[2] = -NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_k[Idx(i, j, nk - 2 - iTemp)][2];
				}
				coef->GetZ(Idx(i, j, nk - 2))[2] = NodeDifferece4thLeft(diff_temp);
				coef->GetZ(Idx(i, j, nk - 3))[2] = -NodeDifferece4th(diff_temp);
				if (grid->GetDim() == TWO_DIM)
				{
					coef->GetX(Idx(i, j, 1))[2] = 0.0;
					coef->GetX(Idx(i, j, 2))[2] = 0.0;
					coef->GetX(Idx(i, j, nk - 2))[2] = 0.0;
					coef->GetX(Idx(i, j, nk - 3))[2] = 0.0;
					coef->GetY(Idx(i, j, 1))[2] = 0.0;
					coef->GetY(Idx(i, j, 2))[2] = 0.0;
					coef->GetY(Idx(i, j, nk - 2))[2] = 0.0;
					coef->GetY(Idx(i, j, nk - 3))[2] = 0.0;
					coef->GetZ(Idx(i, j, 1))[0] = coef->GetZ(Idx(i, j, 1))[1] = 0.0;
					coef->GetZ(Idx(i, j, 2))[0] = coef->GetZ(Idx(i, j, 2))[1] = 0.0;
					coef->GetZ(Idx(i, j, nk - 2))[0] = coef->GetZ(Idx(i, j, nk - 2))[1] = 0.0;
					coef->GetZ(Idx(i, j, nk - 3))[0] = coef->GetZ(Idx(i, j, nk - 3))[1] = 0.0;
					coef->GetZ(Idx(i, j, 1))[2] = 1.0;
					coef->GetZ(Idx(i, j, 2))[2] = 1.0;
					coef->GetZ(Idx(i, j, nk - 2))[2] = 1.0;
					coef->GetZ(Idx(i, j, nk - 3))[2] = 1.0;
				}
			}
		}
		// 计算整点的度量系数
		for (int k = 1; k <= nk - 2; ++k)
		{
			for (int j = 1; j <= nj - 2; ++j)
			{
				for (int i = 1; i <= ni - 2; ++i)
				{
					// 计算xi,eta,zeta
					auto coef_x = GetNodeMetrics()->GetX(idx_proxy->GetIdx(i, j, k));
					auto coef_y = GetNodeMetrics()->GetY(idx_proxy->GetIdx(i, j, k));
					auto coef_z = GetNodeMetrics()->GetZ(idx_proxy->GetIdx(i, j, k));
					auto coef_xi = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, k));
					coef_xi[0] = coef_y[1] * coef_z[2] - coef_z[1] * coef_y[2];
					coef_xi[1] = coef_z[1] * coef_x[2] - coef_x[1] * coef_z[2];
					coef_xi[2] = coef_x[1] * coef_y[2] - coef_y[1] * coef_x[2];
					coef_xi[3] = 0.0;
					auto coef_eta = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, k));
					coef_eta[0] = coef_y[2] * coef_z[0] - coef_z[2] * coef_y[0];
					coef_eta[1] = coef_z[2] * coef_x[0] - coef_x[2] * coef_z[0];
					coef_eta[2] = coef_x[2] * coef_y[0] - coef_y[2] * coef_x[0];
					coef_eta[3] = 0.0;
					auto coef_zeta = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, k));
					coef_zeta[0] = coef_y[0] * coef_z[1] - coef_z[0] * coef_y[1];
					coef_zeta[1] = coef_z[0] * coef_x[1] - coef_x[0] * coef_z[1];
					coef_zeta[2] = coef_x[0] * coef_y[1] - coef_y[0] * coef_x[1];
					coef_zeta[3] = 0.0;
				}
			}
		}
	}

	void NSSolverStruct::CalcCoordTransCoef()
	{
		CalcMetrics();
		CalcJacobian();
	}
	void NSSolverStruct::CalcMetrics()
	{
		if (GetGrid()->GetDim() == TWO_DIM)
		{
			Log::warn("2D case, use standard format metrics");
			CalcMetricsS0();
			return;
		}
		auto para = GetPara();
		auto metrics_type = para->GetMetricsType();
		if (metrics_type == MetricType::S0)
		{
			CalcMetricsS0();
		}
		else if (metrics_type == MetricType::S1)
		{
			CalcMetricsS1();
		}
		else if (metrics_type == MetricType::S2)
		{
			CalcMetricsS2();
		}
		else if (metrics_type == MetricType::S3)
		{
			CalcMetricsS3();
		}
		else
		{
			Log::warn("Metrics Type is not defined, use standard format as default");
			CalcMetricsS0();
		}
	}
	void NSSolverStruct::CalcJacobian()
	{
		if (GetGrid()->GetDim() == TWO_DIM)
		{
			Log::warn("2D case, use standard format Jacobian");
			CalcJacobianV1();
			return;
		}
		auto para = GetPara();
		auto jacobian_scheme = para->GetJacobianType();
		if (jacobian_scheme == JacobianType::V1)
		{
			CalcJacobianV1();
		}
		else if (jacobian_scheme == JacobianType::V2)
		{
			CalcJacobianV2();
		}
		else if (jacobian_scheme == JacobianType::V3)
		{
			CalcJacobianV3();
		}
		else
		{
			Log::warn("Jacobian Scheme is not defined, use standard format as default");
			CalcJacobianV1();
		}
	}
	void NSSolverStruct::CalcJacobianV1()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		int ni = grid->GetNi();
		int nj = grid->GetNj();
		int nk = grid->GetNk();
		auto idx_proxy = GetIdxProxy();
		auto node_metrics = GetNodeMetrics();
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					auto idx = idx_proxy->GetIdx(i, j, k);
					auto& jacobian = node_metrics->GetJacobian(idx);
					auto coef_x = node_metrics->GetX(idx);
					auto coef_y = node_metrics->GetY(idx);
					auto coef_z = node_metrics->GetZ(idx);
					jacobian = coef_x[0] * coef_y[1] * coef_z[2] - coef_x[0] * coef_y[2] * coef_z[1] +
						coef_x[1] * coef_y[2] * coef_z[0] - coef_x[1] * coef_y[0] * coef_z[2] +
						coef_x[2] * coef_y[0] * coef_z[1] - coef_x[2] * coef_y[1] * coef_z[0];
					jacobian = 1.0 / jacobian;
				}
			}
		}
	}
	void NSSolverStruct::CalcJacobianV2()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		int ni = grid->GetNi();
		int nj = grid->GetNj();
		int nk = grid->GetNk();
		auto idx_proxy = GetIdxProxy();
		auto node_metrics = GetNodeMetrics();
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					auto idx = idx_proxy->GetIdx(i, j, k);
					auto& jacobian = node_metrics->GetJacobian(idx);
					auto coef_x = node_metrics->GetX(idx);
					auto coef_y = node_metrics->GetY(idx);
					auto coef_z = node_metrics->GetZ(idx);
					auto coef_xi = node_metrics->GetXi(idx);
					auto coef_eta = node_metrics->GetEta(idx);
					auto coef_zeta = node_metrics->GetZeta(idx);
					jacobian = coef_x[0] * coef_xi[0] + coef_y[0] * coef_xi[1] + coef_z[0] * coef_xi[2] +
						coef_x[1] * coef_eta[0] + coef_y[1] * coef_eta[1] + coef_z[1] * coef_eta[2] +
						coef_x[2] * coef_zeta[0] + coef_y[2] * coef_zeta[1] + coef_z[2] * coef_zeta[2];
					jacobian = jacobian / 3.0;
					jacobian = 1.0 / jacobian;
				}
			}
		}
	}

	void NSSolverStruct::CalcJacobianV3()
	{
		auto para = GetPara();
		auto flux_diff_scheme = para->GetFluxDifferenceScheme();
		if (flux_diff_scheme == FluxDifferenceScheme::SecondOrder)
		{
			CalcJacobianV3_2nd();
		}
		else if (flux_diff_scheme == FluxDifferenceScheme::SixthOrder)
		{
			CalcJacobianV3_6th();
		}
		else
		{
			Log::warn("FluxDifferenceScheme is not defined, use SecondOrder as default");
			CalcJacobianV3_2nd();
		}
	}
	void NSSolverStruct::CalcJacobianV3_2nd()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto coef = GetNodeMetrics();
		auto idx_proxy = GetIdxProxy();
		auto Idx = [&](int i, int j, int k) { return idx_proxy->GetIdx(i, j, k); };
		int ni = grid->GetNi();
		int nj = grid->GetNj();
		int nk = grid->GetNk();
		double inter_temp[2];
		// i+1/2,j+1/2,k+1/2处的临时变量
		dynamic_array<double> temp_i(ni * nj * nk, 0.0), temp_j(ni * nj * nk, 0.0), temp_k(ni * nj * nk, 0.0);
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					int idx0 = Idx(i, j, k);
					int idx1 = Idx(i + 1, j, k);
					auto coord0 = node->GetCoord(i, j, k);
					auto coord1 = node->GetCoord(i + 1, j, k);
					inter_temp[0] = coord0[0] * coef->GetXi(idx0)[0] + coord0[1] * coef->GetXi(idx0)[1] +
						coord0[2] * coef->GetXi(idx0)[2];
					inter_temp[1] = coord1[0] * coef->GetXi(idx1)[0] + coord1[1] * coef->GetXi(idx1)[1] +
						coord1[2] * coef->GetXi(idx1)[2];
					temp_i[idx0] = 0.5 * (inter_temp[0] + inter_temp[1]);
					idx0 = Idx(i, j, k);
					idx1 = Idx(i, j + 1, k);
					coord0 = node->GetCoord(i, j, k);
					coord1 = node->GetCoord(i, j + 1, k);
					inter_temp[0] = coord0[0] * coef->GetEta(idx0)[0] + coord0[1] * coef->GetEta(idx0)[1] +
						coord0[2] * coef->GetEta(idx0)[2];
					inter_temp[1] = coord1[0] * coef->GetEta(idx1)[0] + coord1[1] * coef->GetEta(idx1)[1] +
						coord1[2] * coef->GetEta(idx1)[2];
					temp_j[idx0] = 0.5 * (inter_temp[0] + inter_temp[1]);
					idx0 = Idx(i, j, k);
					idx1 = Idx(i, j, k + 1);
					coord0 = node->GetCoord(i, j, k);
					coord1 = node->GetCoord(i, j, k + 1);
					inter_temp[0] = coord0[0] * coef->GetZeta(idx0)[0] + coord0[1] * coef->GetZeta(idx0)[1] +
						coord0[2] * coef->GetZeta(idx0)[2];
					inter_temp[1] = coord1[0] * coef->GetZeta(idx1)[0] + coord1[1] * coef->GetZeta(idx1)[1] +
						coord1[2] * coef->GetZeta(idx1)[2];
					temp_k[idx0] = 0.5 * (inter_temp[0] + inter_temp[1]);
				}
			}
		}
		// i=1/2,(ni-1)-1/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				int idx0 = Idx(1, j, k);
				int idx1 = Idx(2, j, k);
				auto coord0 = node->GetCoord(1, j, k);
				auto coord1 = node->GetCoord(2, j, k);
				inter_temp[0] =
					coord0[0] * coef->GetXi(idx0)[0] + coord0[1] * coef->GetXi(idx0)[1] + coord0[2] * coef->GetXi(idx0)[2];
				inter_temp[1] =
					coord1[0] * coef->GetXi(idx1)[0] + coord1[1] * coef->GetXi(idx1)[1] + coord1[2] * coef->GetXi(idx1)[2];
				temp_i[idx0] = 0.5 * (3.0 * inter_temp[0] - 1.0 * inter_temp[1]);
				idx0 = Idx(ni - 2, j, k);
				idx1 = Idx(ni - 3, j, k);
				coord0 = node->GetCoord(ni - 2, j, k);
				coord1 = node->GetCoord(ni - 3, j, k);
				inter_temp[0] =
					coord0[0] * coef->GetXi(idx0)[0] + coord0[1] * coef->GetXi(idx0)[1] + coord0[2] * coef->GetXi(idx0)[2];
				inter_temp[1] =
					coord1[0] * coef->GetXi(idx1)[0] + coord1[1] * coef->GetXi(idx1)[1] + coord1[2] * coef->GetXi(idx1)[2];
				temp_i[idx0] = 0.5 * (3.0 * inter_temp[0] - 1.0 * inter_temp[1]);
			}
		}
		// j=1/2,(nj-1)-1/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				int idx0 = Idx(i, 1, k);
				int idx1 = Idx(i, 2, k);
				auto coord0 = node->GetCoord(i, 1, k);
				auto coord1 = node->GetCoord(i, 2, k);
				inter_temp[0] = coord0[0] * coef->GetEta(idx0)[0] + coord0[1] * coef->GetEta(idx0)[1] +
					coord0[2] * coef->GetEta(idx0)[2];
				inter_temp[1] = coord1[0] * coef->GetEta(idx1)[0] + coord1[1] * coef->GetEta(idx1)[1] +
					coord1[2] * coef->GetEta(idx1)[2];
				temp_j[idx0] = 0.5 * (3.0 * inter_temp[0] - 1.0 * inter_temp[1]);
				idx0 = Idx(i, nj - 2, k);
				idx1 = Idx(i, nj - 3, k);
				coord0 = node->GetCoord(i, nj - 2, k);
				coord1 = node->GetCoord(i, nj - 3, k);
				inter_temp[0] = coord0[0] * coef->GetEta(idx0)[0] + coord0[1] * coef->GetEta(idx0)[1] +
					coord0[2] * coef->GetEta(idx0)[2];
				inter_temp[1] = coord1[0] * coef->GetEta(idx1)[0] + coord1[1] * coef->GetEta(idx1)[1] +
					coord1[2] * coef->GetEta(idx1)[2];
				temp_j[idx0] = 0.5 * (3.0 * inter_temp[0] - 1.0 * inter_temp[1]);
			}
		}
		// k=1/2,(nk-1)-1/2
		for (int j = 1; j < nj - 1; ++j)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				int idx0 = Idx(i, j, 1);
				int idx1 = Idx(i, j, 2);
				auto coord0 = node->GetCoord(i, j, 1);
				auto coord1 = node->GetCoord(i, j, 2);
				inter_temp[0] = coord0[0] * coef->GetZeta(idx0)[0] + coord0[1] * coef->GetZeta(idx0)[1] +
					coord0[2] * coef->GetZeta(idx0)[2];
				inter_temp[1] = coord1[0] * coef->GetZeta(idx1)[0] + coord1[1] * coef->GetZeta(idx1)[1] +
					coord1[2] * coef->GetZeta(idx1)[2];
				temp_k[idx0] = 0.5 * (3.0 * inter_temp[0] - 1.0 * inter_temp[1]);
				idx0 = Idx(i, j, nk - 2);
				idx1 = Idx(i, j, nk - 3);
				coord0 = node->GetCoord(i, j, nk - 2);
				coord1 = node->GetCoord(i, j, nk - 3);
				inter_temp[0] = coord0[0] * coef->GetZeta(idx0)[0] + coord0[1] * coef->GetZeta(idx0)[1] +
					coord0[2] * coef->GetZeta(idx0)[2];
				inter_temp[1] = coord1[0] * coef->GetZeta(idx1)[0] + coord1[1] * coef->GetZeta(idx1)[1] +
					coord1[2] * coef->GetZeta(idx1)[2];
				temp_k[idx0] = 0.5 * (3.0 * inter_temp[0] - 1.0 * inter_temp[1]);
			}
		}
		// 计算整点的Jacobian
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					int idx = Idx(i, j, k);
					auto& jacobian = coef->GetJacobian(idx);
					jacobian = temp_i[Idx(i, j, k)] - temp_i[Idx(i - 1, j, k)] + temp_j[Idx(i, j, k)] -
						temp_j[Idx(i, j - 1, k)] + temp_k[Idx(i, j, k)] - temp_k[Idx(i, j, k - 1)];
					jacobian = jacobian / 3.0;
					jacobian = 1.0 / jacobian;
				}
			}
		}
	}
	void NSSolverStruct::CalcJacobianV3_6th()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto coef = GetNodeMetrics();
		auto idx_proxy = GetIdxProxy();
		auto Idx = [&](int i, int j, int k) { return idx_proxy->GetIdx(i, j, k); };
		int ni = grid->GetNi();
		int nj = grid->GetNj();
		int nk = grid->GetNk();
		double inter_temp[6];
		// i+1/2,j+1/2,k+1/2处的临时变量
		dynamic_array<double> temp_i(ni * nj * nk, 0.0), temp_j(ni * nj * nk, 0.0), temp_k(ni * nj * nk, 0.0);
		// i+1/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 3; i < ni - 4; ++i)
				{
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							node->GetCoord(i + iTemp - 2, j, k)[0] * coef->GetXi(Idx(i + iTemp - 2, j, k))[0] +
							node->GetCoord(i + iTemp - 2, j, k)[1] * coef->GetXi(Idx(i + iTemp - 2, j, k))[1] +
							node->GetCoord(i + iTemp - 2, j, k)[2] * coef->GetXi(Idx(i + iTemp - 2, j, k))[2];
					}
					temp_i[Idx(i, j, k)] = MidNodeInter6th(inter_temp);
				}
			}
		}
		// i=1/2,3/2,5/2;ni-3/2,ni-5/2,ni-7/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = node->GetCoord(1 + iTemp, j, k)[0] * coef->GetXi(Idx(1 + iTemp, j, k))[0] +
						node->GetCoord(1 + iTemp, j, k)[1] * coef->GetXi(Idx(1 + iTemp, j, k))[1] +
						node->GetCoord(1 + iTemp, j, k)[2] * coef->GetXi(Idx(1 + iTemp, j, k))[2];
				}
				temp_i[Idx(0, j, k)] = MidNodeInter4thRight1(inter_temp);
				temp_i[Idx(1, j, k)] = MidNodeInter4thRight2(inter_temp);
				temp_i[Idx(2, j, k)] = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] =
						node->GetCoord(ni - 5 + iTemp, j, k)[0] * coef->GetXi(Idx(ni - 5 + iTemp, j, k))[0] +
						node->GetCoord(ni - 5 + iTemp, j, k)[1] * coef->GetXi(Idx(ni - 5 + iTemp, j, k))[1] +
						node->GetCoord(ni - 5 + iTemp, j, k)[2] * coef->GetXi(Idx(ni - 5 + iTemp, j, k))[2];
				}
				temp_i[Idx(ni - 2, j, k)] = MidNodeInter4thLeft1(inter_temp);
				temp_i[Idx(ni - 3, j, k)] = MidNodeInter4thLeft2(inter_temp);
				temp_i[Idx(ni - 4, j, k)] = MidNodeInter4th(inter_temp);
			}
		}
		// j+1/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 3; j < nj - 4; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							node->GetCoord(i, j + iTemp - 2, k)[0] * coef->GetEta(Idx(i, j + iTemp - 2, k))[0] +
							node->GetCoord(i, j + iTemp - 2, k)[1] * coef->GetEta(Idx(i, j + iTemp - 2, k))[1] +
							node->GetCoord(i, j + iTemp - 2, k)[2] * coef->GetEta(Idx(i, j + iTemp - 2, k))[2];
					}
					temp_j[Idx(i, j, k)] = MidNodeInter6th(inter_temp);
				}
			}
		}
		// j=1/2,3/2,5/2;nj-3/2,nj-5/2,nj-7/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = node->GetCoord(i, 1 + iTemp, k)[0] * coef->GetEta(Idx(i, 1 + iTemp, k))[0] +
						node->GetCoord(i, 1 + iTemp, k)[1] * coef->GetEta(Idx(i, 1 + iTemp, k))[1] +
						node->GetCoord(i, 1 + iTemp, k)[2] * coef->GetEta(Idx(i, 1 + iTemp, k))[2];
				}
				temp_j[Idx(i, 0, k)] = MidNodeInter4thRight1(inter_temp);
				temp_j[Idx(i, 1, k)] = MidNodeInter4thRight2(inter_temp);
				temp_j[Idx(i, 2, k)] = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] =
						node->GetCoord(i, nj - 5 + iTemp, k)[0] * coef->GetEta(Idx(i, nj - 5 + iTemp, k))[0] +
						node->GetCoord(i, nj - 5 + iTemp, k)[1] * coef->GetEta(Idx(i, nj - 5 + iTemp, k))[1] +
						node->GetCoord(i, nj - 5 + iTemp, k)[2] * coef->GetEta(Idx(i, nj - 5 + iTemp, k))[2];
				}
				temp_j[Idx(i, nj - 2, k)] = MidNodeInter4thLeft1(inter_temp);
				temp_j[Idx(i, nj - 3, k)] = MidNodeInter4thLeft2(inter_temp);
				temp_j[Idx(i, nj - 4, k)] = MidNodeInter4th(inter_temp);
			}
		}
		// k+1/2
		for (int k = 3; k < nk - 4; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							node->GetCoord(i, j, k + iTemp - 2)[0] * coef->GetZeta(Idx(i, j, k + iTemp - 2))[0] +
							node->GetCoord(i, j, k + iTemp - 2)[1] * coef->GetZeta(Idx(i, j, k + iTemp - 2))[1] +
							node->GetCoord(i, j, k + iTemp - 2)[2] * coef->GetZeta(Idx(i, j, k + iTemp - 2))[2];
					}
					temp_k[Idx(i, j, k)] = MidNodeInter6th(inter_temp);
				}
			}
		}
		// k=1/2,3/2,5/2;nk-3/2,nk-5/2,nk-7/2
		for (int j = 1; j < nj - 1; ++j)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = node->GetCoord(i, j, 1 + iTemp)[0] * coef->GetZeta(Idx(i, j, 1 + iTemp))[0] +
						node->GetCoord(i, j, 1 + iTemp)[1] * coef->GetZeta(Idx(i, j, 1 + iTemp))[1] +
						node->GetCoord(i, j, 1 + iTemp)[2] * coef->GetZeta(Idx(i, j, 1 + iTemp))[2];
				}
				temp_k[Idx(i, j, 0)] = MidNodeInter4thRight1(inter_temp);
				temp_k[Idx(i, j, 1)] = MidNodeInter4thRight2(inter_temp);
				temp_k[Idx(i, j, 2)] = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] =
						node->GetCoord(i, j, nk - 5 + iTemp)[0] * coef->GetZeta(Idx(i, j, nk - 5 + iTemp))[0] +
						node->GetCoord(i, j, nk - 5 + iTemp)[1] * coef->GetZeta(Idx(i, j, nk - 5 + iTemp))[1] +
						node->GetCoord(i, j, nk - 5 + iTemp)[2] * coef->GetZeta(Idx(i, j, nk - 5 + iTemp))[2];
				}
				temp_k[Idx(i, j, nk - 2)] = MidNodeInter4thLeft1(inter_temp);
				temp_k[Idx(i, j, nk - 3)] = MidNodeInter4thLeft2(inter_temp);
				temp_k[Idx(i, j, nk - 4)] = MidNodeInter4th(inter_temp);
			}
		}
		// 计算整点的Jacobian
		// 全部置零
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					coef->GetJacobian(Idx(i, j, k)) = 0.0;
				}
			}
		}
		double diff_temp[6];
		// i direction
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 3; i < ni - 3; ++i)
				{
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = temp_i[Idx(i + iTemp - 3, j, k)];
					}
					coef->GetJacobian(Idx(i, j, k)) += NodeDifferece6th(diff_temp);
				}
			}
		}
		// j direction
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 3; j < nj - 3; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = temp_j[Idx(i, j + iTemp - 3, k)];
					}
					coef->GetJacobian(Idx(i, j, k)) += NodeDifferece6th(diff_temp);
				}
			}
		}
		// k direction
		for (int k = 3; k < nk - 3; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = temp_k[Idx(i, j, k + iTemp - 3)];
					}
					coef->GetJacobian(Idx(i, j, k)) += NodeDifferece6th(diff_temp);
				}
			}
		}
		// 计算边界整点的Jacobian
		// i=1,2,ni-3,ni-2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = temp_i[Idx(0 + iTemp, j, k)];
				}
				coef->GetJacobian(Idx(1, j, k)) += NodeDifferece4thRight(diff_temp);
				coef->GetJacobian(Idx(2, j, k)) += NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = temp_i[Idx(ni - 2 - iTemp, j, k)];
				}
				coef->GetJacobian(Idx(ni - 2, j, k)) += NodeDifferece4thLeft(diff_temp);
				coef->GetJacobian(Idx(ni - 3, j, k)) += -NodeDifferece4th(diff_temp);
			}
		}
		// j=1,2,nj-3,nj-2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = temp_j[Idx(i, 0 + iTemp, k)];
				}
				coef->GetJacobian(Idx(i, 1, k)) += NodeDifferece4thRight(diff_temp);
				coef->GetJacobian(Idx(i, 2, k)) += NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = temp_j[Idx(i, nj - 2 - iTemp, k)];
				}
				coef->GetJacobian(Idx(i, nj - 2, k)) += NodeDifferece4thLeft(diff_temp);
				coef->GetJacobian(Idx(i, nj - 3, k)) += -NodeDifferece4th(diff_temp);
			}
		}
		// k=1,2,nk-3,nk-2
		for (int j = 1; j < nj - 1; ++j)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = temp_k[Idx(i, j, 0 + iTemp)];
				}
				coef->GetJacobian(Idx(i, j, 1)) += NodeDifferece4thRight(diff_temp);
				coef->GetJacobian(Idx(i, j, 2)) += NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = temp_k[Idx(i, j, nk - 2 - iTemp)];
				}
				coef->GetJacobian(Idx(i, j, nk - 2)) += NodeDifferece4thLeft(diff_temp);
				coef->GetJacobian(Idx(i, j, nk - 3)) += -NodeDifferece4th(diff_temp);
			}
		}
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					coef->GetJacobian(Idx(i, j, k)) = coef->GetJacobian(Idx(i, j, k)) / 3.0;
					coef->GetJacobian(Idx(i, j, k)) = 1.0 / coef->GetJacobian(Idx(i, j, k));
				}
			}
		}
	}
	void NSSolverStruct::CalcMetricsS1()
	{
		auto para = GetPara();
		auto metric_diff_scheme = para->GetMetricDifferenceScheme();
		if (metric_diff_scheme == MetricDifferenceScheme::SecondOrder)
		{
			CalcMetricsS1_2nd();
		}
		else if (metric_diff_scheme == MetricDifferenceScheme::SixthOrder)
		{
			CalcMetricsS1_6th();
		}
		else
		{
			Log::warn("Flux Difference Scheme is not defined, use SecondOrder as default");
			CalcMetricsS1_2nd();
		}
	}
	void NSSolverStruct::CalcMetricsS1_2nd()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto coef = GetNodeMetrics();
		Metric* coef_mid[3] = { GetMidMetricsI(), GetMidMetricsJ(), GetMidMetricsK() };
		auto idx_proxy = GetIdxProxy();
		auto Idx = [&](int i, int j, int k) { return idx_proxy->GetIdx(i, j, k); };
		int ni = grid->GetNi();
		int nj = grid->GetNj();
		int nk = grid->GetNk();
		int idx_temp[3];
		// i+1/2,j+1/2,k+1/2处的坐标
		dynamic_array<dynamic_array<double>> coord_i(ni * nj * nk, dynamic_array<double>(3)),
			coord_j(ni * nj * nk, dynamic_array<double>(3)), coord_k(ni * nj * nk, dynamic_array<double>(3));
		// 计算之前，先把度量系数赋值为0
		for (int i = 0; i < ni; ++i)
		{
			for (int j = 0; j < nj; ++j)
			{
				for (int k = 0; k < nk; ++k)
				{
					for (int iDim = 0; iDim < 4; ++iDim)
					{
						coef->GetX(Idx(i, j, k))[iDim] = 0.0;
						coef->GetY(Idx(i, j, k))[iDim] = 0.0;
						coef->GetZ(Idx(i, j, k))[iDim] = 0.0;
						coef->GetXi(Idx(i, j, k))[iDim] = 0.0;
						coef->GetEta(Idx(i, j, k))[iDim] = 0.0;
						coef->GetZeta(Idx(i, j, k))[iDim] = 0.0;
						for (int jDim = 0; jDim < 3; ++jDim)
						{
							coef_mid[jDim]->GetX(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetY(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetZ(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetXi(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetEta(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetZeta(Idx(i, j, k))[iDim] = 0.0;
						}
					}
				}
			}
		}
		// 第一步：计算半点坐标
		for (int k = 0; k < nk - 1; ++k)
		{
			for (int j = 0; j < nj - 1; ++j)
			{
				for (int i = 0; i < ni - 1; ++i)
				{
					for (int iDim = 0; iDim < 3; ++iDim)
					{
						coord_i[Idx(i, j, k)][iDim] =
							0.5 * (node->GetCoord(i, j, k)[iDim] + node->GetCoord(i + 1, j, k)[iDim]);
						coord_j[Idx(i, j, k)][iDim] =
							0.5 * (node->GetCoord(i, j, k)[iDim] + node->GetCoord(i, j + 1, k)[iDim]);
						coord_k[Idx(i, j, k)][iDim] =
							0.5 * (node->GetCoord(i, j, k)[iDim] + node->GetCoord(i, j, k + 1)[iDim]);
					}
				}
			}
		}
		// 第二步：根据半点坐标计算整点逆变换度量系数
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					int idx = idx_proxy->GetIdx(i, j, k);
					auto coef_x = coef->GetX(idx);
					auto coef_y = coef->GetY(idx);
					auto coef_z = coef->GetZ(idx);
					coef_x[0] = coord_i[Idx(i, j, k)][0] - coord_i[Idx(i - 1, j, k)][0];
					coef_x[1] = coord_j[Idx(i, j, k)][0] - coord_j[Idx(i, j - 1, k)][0];
					coef_x[2] = coord_k[Idx(i, j, k)][0] - coord_k[Idx(i, j, k - 1)][0];
					coef_y[0] = coord_i[Idx(i, j, k)][1] - coord_i[Idx(i - 1, j, k)][1];
					coef_y[1] = coord_j[Idx(i, j, k)][1] - coord_j[Idx(i, j - 1, k)][1];
					coef_y[2] = coord_k[Idx(i, j, k)][1] - coord_k[Idx(i, j, k - 1)][1];
					coef_z[0] = coord_i[Idx(i, j, k)][2] - coord_i[Idx(i - 1, j, k)][2];
					coef_z[1] = coord_j[Idx(i, j, k)][2] - coord_j[Idx(i, j - 1, k)][2];
					coef_z[2] = coord_k[Idx(i, j, k)][2] - coord_k[Idx(i, j, k - 1)][2];
					coef_x[3] = 0.0;
					coef_y[3] = 0.0;
					coef_z[3] = 0.0;
					if (grid->GetDim() == TWO_DIM)
					{
						coef_z[0] = 0.0;
						coef_z[1] = 0.0;
						coef_z[2] = 1.0;
					}
				}
			}
		}
		// 第三步：计算半点处逆变换度量系数和坐标乘积
		//  求出xi方向的临时变量：逆变换度量系数乘以坐标
		struct TempXi
		{
			double y_zeta_z, z_zeta_x, x_zeta_y, y_eta_z, z_eta_x, x_eta_y;
		};
		dynamic_array<TempXi> temp_xi(ni * nj * nk);
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					int idx = Idx(i, j, k);
					auto& temp = temp_xi[idx];
					temp.y_zeta_z = 0.5 * (coef->GetY(Idx(i + 1, j, k))[2] * node->GetCoord(i + 1, j, k)[2] +
						coef->GetY(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[2]);
					temp.z_zeta_x = 0.5 * (coef->GetZ(Idx(i + 1, j, k))[2] * node->GetCoord(i + 1, j, k)[0] +
						coef->GetZ(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[0]);
					temp.x_zeta_y = 0.5 * (coef->GetX(Idx(i + 1, j, k))[2] * node->GetCoord(i + 1, j, k)[1] +
						coef->GetX(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[1]);
					temp.y_eta_z = 0.5 * (coef->GetY(Idx(i + 1, j, k))[1] * node->GetCoord(i + 1, j, k)[2] +
						coef->GetY(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[2]);
					temp.z_eta_x = 0.5 * (coef->GetZ(Idx(i + 1, j, k))[1] * node->GetCoord(i + 1, j, k)[0] +
						coef->GetZ(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[0]);
					temp.x_eta_y = 0.5 * (coef->GetX(Idx(i + 1, j, k))[1] * node->GetCoord(i + 1, j, k)[1] +
						coef->GetX(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[1]);
				}
			}
		}
		// i=1/2；(ni-1)-1/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				int idx = Idx(0, j, k);
				temp_xi[idx].y_zeta_z = 0.5 * (3.0 * coef->GetY(Idx(1, j, k))[2] * node->GetCoord(1, j, k)[2] -
					1.0 * coef->GetY(Idx(2, j, k))[2] * node->GetCoord(2, j, k)[2]);
				temp_xi[idx].z_zeta_x = 0.5 * (3.0 * coef->GetZ(Idx(1, j, k))[2] * node->GetCoord(1, j, k)[0] -
					1.0 * coef->GetZ(Idx(2, j, k))[2] * node->GetCoord(2, j, k)[0]);
				temp_xi[idx].x_zeta_y = 0.5 * (3.0 * coef->GetX(Idx(1, j, k))[2] * node->GetCoord(1, j, k)[1] -
					1.0 * coef->GetX(Idx(2, j, k))[2] * node->GetCoord(2, j, k)[1]);
				temp_xi[idx].y_eta_z = 0.5 * (3.0 * coef->GetY(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[2] -
					1.0 * coef->GetY(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[2]);
				temp_xi[idx].z_eta_x = 0.5 * (3.0 * coef->GetZ(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[0] -
					1.0 * coef->GetZ(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[0]);
				temp_xi[idx].x_eta_y = 0.5 * (3.0 * coef->GetX(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[1] -
					1.0 * coef->GetX(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[1]);
				idx = Idx(ni - 2, j, k);
				temp_xi[idx].y_zeta_z = 0.5 * (3.0 * coef->GetY(Idx(ni - 2, j, k))[2] * node->GetCoord(ni - 2, j, k)[2] -
					1.0 * coef->GetY(Idx(ni - 3, j, k))[2] * node->GetCoord(ni - 3, j, k)[2]);
				temp_xi[idx].z_zeta_x = 0.5 * (3.0 * coef->GetZ(Idx(ni - 2, j, k))[2] * node->GetCoord(ni - 2, j, k)[0] -
					1.0 * coef->GetZ(Idx(ni - 3, j, k))[2] * node->GetCoord(ni - 3, j, k)[0]);
				temp_xi[idx].x_zeta_y = 0.5 * (3.0 * coef->GetX(Idx(ni - 2, j, k))[2] * node->GetCoord(ni - 2, j, k)[1] -
					1.0 * coef->GetX(Idx(ni - 3, j, k))[2] * node->GetCoord(ni - 3, j, k)[1]);
				temp_xi[idx].y_eta_z = 0.5 * (3.0 * coef->GetY(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[2] -
					1.0 * coef->GetY(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[2]);
				temp_xi[idx].z_eta_x = 0.5 * (3.0 * coef->GetZ(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[0] -
					1.0 * coef->GetZ(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[0]);
				temp_xi[idx].x_eta_y = 0.5 * (3.0 * coef->GetX(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[1] -
					1.0 * coef->GetX(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[1]);
			}
		}
		// 求出eta方向的临时变量：逆变换度量系数乘以坐标
		struct TempEta
		{
			double y_zeta_z, z_zeta_x, x_zeta_y, y_xi_z, z_xi_x, x_xi_y;
		};
		dynamic_array<TempEta> temp_eta(ni * nj * nk);
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					int idx = Idx(i, j, k);
					auto& temp = temp_eta[idx];
					temp.y_zeta_z = 0.5 * (coef->GetY(Idx(i, j + 1, k))[2] * node->GetCoord(i, j + 1, k)[2] +
						coef->GetY(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[2]);
					temp.z_zeta_x = 0.5 * (coef->GetZ(Idx(i, j + 1, k))[2] * node->GetCoord(i, j + 1, k)[0] +
						coef->GetZ(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[0]);
					temp.x_zeta_y = 0.5 * (coef->GetX(Idx(i, j + 1, k))[2] * node->GetCoord(i, j + 1, k)[1] +
						coef->GetX(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[1]);
					temp.y_xi_z = 0.5 * (coef->GetY(Idx(i, j + 1, k))[0] * node->GetCoord(i, j + 1, k)[2] +
						coef->GetY(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[2]);
					temp.z_xi_x = 0.5 * (coef->GetZ(Idx(i, j + 1, k))[0] * node->GetCoord(i, j + 1, k)[0] +
						coef->GetZ(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[0]);
					temp.x_xi_y = 0.5 * (coef->GetX(Idx(i, j + 1, k))[0] * node->GetCoord(i, j + 1, k)[1] +
						coef->GetX(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[1]);
				}
			}
		}
		// j=1/2；(nj-1)-1/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				int idx = Idx(i, 0, k);
				temp_eta[idx].y_zeta_z = 0.5 * (3.0 * coef->GetY(Idx(i, 1, k))[2] * node->GetCoord(i, 1, k)[2] -
					1.0 * coef->GetY(Idx(i, 2, k))[2] * node->GetCoord(i, 2, k)[2]);
				temp_eta[idx].z_zeta_x = 0.5 * (3.0 * coef->GetZ(Idx(i, 1, k))[2] * node->GetCoord(i, 1, k)[0] -
					1.0 * coef->GetZ(Idx(i, 2, k))[2] * node->GetCoord(i, 2, k)[0]);
				temp_eta[idx].x_zeta_y = 0.5 * (3.0 * coef->GetX(Idx(i, 1, k))[2] * node->GetCoord(i, 1, k)[1] -
					1.0 * coef->GetX(Idx(i, 2, k))[2] * node->GetCoord(i, 2, k)[1]);
				temp_eta[idx].y_xi_z = 0.5 * (3.0 * coef->GetY(Idx(i, 1, k))[0] * node->GetCoord(i, 1, k)[2] -
					1.0 * coef->GetY(Idx(i, 2, k))[0] * node->GetCoord(i, 2, k)[2]);
				temp_eta[idx].z_xi_x = 0.5 * (3.0 * coef->GetZ(Idx(i, 1, k))[0] * node->GetCoord(i, 1, k)[0] -
					1.0 * coef->GetZ(Idx(i, 2, k))[0] * node->GetCoord(i, 2, k)[0]);
				temp_eta[idx].x_xi_y = 0.5 * (3.0 * coef->GetX(Idx(i, 1, k))[0] * node->GetCoord(i, 1, k)[1] -
					1.0 * coef->GetX(Idx(i, 2, k))[0] * node->GetCoord(i, 2, k)[1]);
				idx = Idx(i, nj - 2, k);
				temp_eta[idx].y_zeta_z = 0.5 * (3.0 * coef->GetY(Idx(i, nj - 2, k))[2] * node->GetCoord(i, nj - 2, k)[2] -
					1.0 * coef->GetY(Idx(i, nj - 3, k))[2] * node->GetCoord(i, nj - 3, k)[2]);
				temp_eta[idx].z_zeta_x = 0.5 * (3.0 * coef->GetZ(Idx(i, nj - 2, k))[2] * node->GetCoord(i, nj - 2, k)[0] -
					1.0 * coef->GetZ(Idx(i, nj - 3, k))[2] * node->GetCoord(i, nj - 3, k)[0]);
				temp_eta[idx].x_zeta_y = 0.5 * (3.0 * coef->GetX(Idx(i, nj - 2, k))[2] * node->GetCoord(i, nj - 2, k)[1] -
					1.0 * coef->GetX(Idx(i, nj - 3, k))[2] * node->GetCoord(i, nj - 3, k)[1]);
				temp_eta[idx].y_xi_z = 0.5 * (3.0 * coef->GetY(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[2] -
					1.0 * coef->GetY(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[2]);
				temp_eta[idx].z_xi_x = 0.5 * (3.0 * coef->GetZ(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[0] -
					1.0 * coef->GetZ(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[0]);
				temp_eta[idx].x_xi_y = 0.5 * (3.0 * coef->GetX(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[1] -
					1.0 * coef->GetX(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[1]);
			}
		}

		// 求出zeta方向的临时变量：逆变换度量系数乘以坐标
		struct TempZeta
		{
			double y_eta_z, z_eta_x, x_eta_y, y_xi_z, z_xi_x, x_xi_y;
		};
		dynamic_array<TempZeta> temp_zeta(ni * nj * nk);
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					int idx = Idx(i, j, k);
					auto& temp = temp_zeta[idx];
					temp.y_eta_z = 0.5 * (coef->GetY(Idx(i, j, k + 1))[1] * node->GetCoord(i, j, k + 1)[2] +
						coef->GetY(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[2]);
					temp.z_eta_x = 0.5 * (coef->GetZ(Idx(i, j, k + 1))[1] * node->GetCoord(i, j, k + 1)[0] +
						coef->GetZ(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[0]);
					temp.x_eta_y = 0.5 * (coef->GetX(Idx(i, j, k + 1))[1] * node->GetCoord(i, j, k + 1)[1] +
						coef->GetX(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[1]);
					temp.y_xi_z = 0.5 * (coef->GetY(Idx(i, j, k + 1))[0] * node->GetCoord(i, j, k + 1)[2] +
						coef->GetY(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[2]);
					temp.z_xi_x = 0.5 * (coef->GetZ(Idx(i, j, k + 1))[0] * node->GetCoord(i, j, k + 1)[0] +
						coef->GetZ(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[0]);
					temp.x_xi_y = 0.5 * (coef->GetX(Idx(i, j, k + 1))[0] * node->GetCoord(i, j, k + 1)[1] +
						coef->GetX(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[1]);
				}
			}
		}
		// k=1/2；(nk-1)-1/2
		for (int j = 1; j < nj - 1; ++j)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				int idx = Idx(i, j, 0);
				temp_zeta[idx].y_eta_z = 0.5 * (3.0 * coef->GetY(Idx(i, j, 1))[1] * node->GetCoord(i, j, 1)[2] -
					1.0 * coef->GetY(Idx(i, j, 2))[1] * node->GetCoord(i, j, 2)[2]);
				temp_zeta[idx].z_eta_x = 0.5 * (3.0 * coef->GetZ(Idx(i, j, 1))[1] * node->GetCoord(i, j, 1)[0] -
					1.0 * coef->GetZ(Idx(i, j, 2))[1] * node->GetCoord(i, j, 2)[0]);
				temp_zeta[idx].x_eta_y = 0.5 * (3.0 * coef->GetX(Idx(i, j, 1))[1] * node->GetCoord(i, j, 1)[1] -
					1.0 * coef->GetX(Idx(i, j, 2))[1] * node->GetCoord(i, j, 2)[1]);
				temp_zeta[idx].y_xi_z = 0.5 * (3.0 * coef->GetY(Idx(i, j, 1))[0] * node->GetCoord(i, j, 1)[2] -
					1.0 * coef->GetY(Idx(i, j, 2))[0] * node->GetCoord(i, j, 2)[2]);
				temp_zeta[idx].z_xi_x = 0.5 * (3.0 * coef->GetZ(Idx(i, j, 1))[0] * node->GetCoord(i, j, 1)[0] -
					1.0 * coef->GetZ(Idx(i, j, 2))[0] * node->GetCoord(i, j, 2)[0]);
				temp_zeta[idx].x_xi_y = 0.5 * (3.0 * coef->GetX(Idx(i, j, 1))[0] * node->GetCoord(i, j, 1)[1] -
					1.0 * coef->GetX(Idx(i, j, 2))[0] * node->GetCoord(i, j, 2)[1]);
				idx = Idx(i, j, nk - 2);
				temp_zeta[idx].y_eta_z = 0.5 * (3.0 * coef->GetY(Idx(i, j, nk - 2))[1] * node->GetCoord(i, j, nk - 2)[2] -
					1.0 * coef->GetY(Idx(i, j, nk - 3))[1] * node->GetCoord(i, j, nk - 3)[2]);
				temp_zeta[idx].z_eta_x = 0.5 * (3.0 * coef->GetZ(Idx(i, j, nk - 2))[1] * node->GetCoord(i, j, nk - 2)[0] -
					1.0 * coef->GetZ(Idx(i, j, nk - 3))[1] * node->GetCoord(i, j, nk - 3)[0]);
				temp_zeta[idx].x_eta_y = 0.5 * (3.0 * coef->GetX(Idx(i, j, nk - 2))[1] * node->GetCoord(i, j, nk - 2)[1] -
					1.0 * coef->GetX(Idx(i, j, nk - 3))[1] * node->GetCoord(i, j, nk - 3)[1]);
				temp_zeta[idx].y_xi_z = 0.5 * (3.0 * coef->GetY(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[2] -
					1.0 * coef->GetY(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[2]);
				temp_zeta[idx].z_xi_x = 0.5 * (3.0 * coef->GetZ(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[0] -
					1.0 * coef->GetZ(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[0]);
				temp_zeta[idx].x_xi_y = 0.5 * (3.0 * coef->GetX(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[1] -
					1.0 * coef->GetX(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[1]);
			}
		}

		// 第四步：根据半点坐标和半点逆变换度量系数使用守恒形式计算整点度量系数（CMM1)
		double temp[6][2];
		// 获取下标的lamda函数
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					// i direction
					auto coef_xi = coef->GetXi(Idx(i, j, k));
					auto coef_eta = coef->GetEta(Idx(i, j, k));
					auto coef_zeta = coef->GetZeta(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 2; iTemp++)
					{
						int idx = Idx(i - 1 + iTemp, j, k);
						temp[0][iTemp] = temp_xi[idx].y_zeta_z; // y_zeta*z
						temp[1][iTemp] = temp_xi[idx].z_zeta_x; // z_zeta*x
						temp[2][iTemp] = temp_xi[idx].x_zeta_y; // x_zeta*y
						temp[3][iTemp] = -temp_xi[idx].y_eta_z; // -y_eta*z
						temp[4][iTemp] = -temp_xi[idx].z_eta_x; // -z_eta*x
						temp[5][iTemp] = -temp_xi[idx].x_eta_y; // -x_eta*y
					}
					coef_eta[0] += temp[0][1] - temp[0][0];
					coef_eta[1] += temp[1][1] - temp[1][0];
					coef_eta[2] += temp[2][1] - temp[2][0];
					coef_zeta[0] += temp[3][1] - temp[3][0];
					coef_zeta[1] += temp[4][1] - temp[4][0];
					coef_zeta[2] += temp[5][1] - temp[5][0];
					for (int iTemp = 0; iTemp < 2; iTemp++)
					{
						int idx = Idx(i, j - 1 + iTemp, k);
						temp[0][iTemp] = -temp_eta[idx].y_zeta_z; // -y_zeta*z
						temp[1][iTemp] = -temp_eta[idx].z_zeta_x; // -z_zeta*x
						temp[2][iTemp] = -temp_eta[idx].x_zeta_y; // -x_zeta*y
						temp[3][iTemp] = temp_eta[idx].y_xi_z;    // y_xi*z
						temp[4][iTemp] = temp_eta[idx].z_xi_x;    // z_xi*x
						temp[5][iTemp] = temp_eta[idx].x_xi_y;    // x_xi*y
					}
					coef_xi[0] += temp[0][1] - temp[0][0];
					coef_xi[1] += temp[1][1] - temp[1][0];
					coef_xi[2] += temp[2][1] - temp[2][0];
					coef_zeta[0] += temp[3][1] - temp[3][0];
					coef_zeta[1] += temp[4][1] - temp[4][0];
					coef_zeta[2] += temp[5][1] - temp[5][0];
					for (int iTemp = 0; iTemp < 2; iTemp++)
					{
						int idx = Idx(i, j, k - 1 + iTemp);
						temp[0][iTemp] = temp_zeta[idx].y_eta_z; // y_eta*z
						temp[1][iTemp] = temp_zeta[idx].z_eta_x; // z_eta*x
						temp[2][iTemp] = temp_zeta[idx].x_eta_y; // x_eta*y
						temp[3][iTemp] = -temp_zeta[idx].y_xi_z; // -y_xi*z
						temp[4][iTemp] = -temp_zeta[idx].z_xi_x; // -z_xi*x
						temp[5][iTemp] = -temp_zeta[idx].x_xi_y; // -x_xi*y
					}
					coef_xi[0] += temp[0][1] - temp[0][0];
					coef_xi[1] += temp[1][1] - temp[1][0];
					coef_xi[2] += temp[2][1] - temp[2][0];
					coef_eta[0] += temp[3][1] - temp[3][0];
					coef_eta[1] += temp[4][1] - temp[4][0];
					coef_eta[2] += temp[5][1] - temp[5][0];
					if (grid->GetDim() == TWO_DIM)
					{
						coef_zeta[0] = 0.0;
						coef_zeta[1] = 0.0;
						coef_zeta[2] = 1.0;
					}
					coef_xi[3] = coef_eta[3] = coef_zeta[3] = 0.0;
				}
			}
		}
	}
	void NSSolverStruct::CalcMetricsS1_6th()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto coef = GetNodeMetrics();
		Metric* coef_mid[3] = { GetMidMetricsI(), GetMidMetricsJ(), GetMidMetricsK() };
		auto idx_proxy = GetIdxProxy();
		auto Idx = [&](int i, int j, int k) { return idx_proxy->GetIdx(i, j, k); };
		int ni = grid->GetNi();
		int nj = grid->GetNj();
		int nk = grid->GetNk();
		double inter_temp[6];
		// i+1/2,j+1/2,k+1/2处的坐标
		dynamic_array<dynamic_array<double>> coord_i(ni * nj * nk, dynamic_array<double>(3)),
			coord_j(ni * nj * nk, dynamic_array<double>(3)), coord_k(ni * nj * nk, dynamic_array<double>(3));
		// 计算之前，先把度量系数赋值为0
		for (int i = 0; i < ni; ++i)
		{
			for (int j = 0; j < nj; ++j)
			{
				for (int k = 0; k < nk; ++k)
				{
					for (int iDim = 0; iDim < 4; ++iDim)
					{
						coef->GetX(Idx(i, j, k))[iDim] = 0.0;
						coef->GetY(Idx(i, j, k))[iDim] = 0.0;
						coef->GetZ(Idx(i, j, k))[iDim] = 0.0;
						coef->GetXi(Idx(i, j, k))[iDim] = 0.0;
						coef->GetEta(Idx(i, j, k))[iDim] = 0.0;
						coef->GetZeta(Idx(i, j, k))[iDim] = 0.0;
						for (int jDim = 0; jDim < 3; ++jDim)
						{
							coef_mid[jDim]->GetX(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetY(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetZ(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetXi(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetEta(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetZeta(Idx(i, j, k))[iDim] = 0.0;
						}
					}
				}
			}
		}
		// 第一步：计算半点坐标
		// i+1/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 3; i < ni - 4; ++i)
				{
					for (int iDim = 0; iDim < 3; ++iDim)
					{
						for (int iTemp = 0; iTemp < 6; iTemp++)
						{
							inter_temp[iTemp] = node->GetCoord(i + iTemp - 2, j, k)[iDim];
						}
						coord_i[Idx(i, j, k)][iDim] = MidNodeInter6th(inter_temp);
					}
				}
			}
		}
		// i=1/2,3/2,5/2;ni-3/2,ni-5/2,ni-7/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int iDim = 0; iDim < 3; ++iDim)
				{
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = node->GetCoord(1 + iTemp, j, k)[iDim];
					}
					coord_i[Idx(0, j, k)][iDim] = MidNodeInter4thRight1(inter_temp);
					coord_i[Idx(1, j, k)][iDim] = MidNodeInter4thRight2(inter_temp);
					coord_i[Idx(2, j, k)][iDim] = MidNodeInter4th(inter_temp);
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = node->GetCoord(ni - 5 + iTemp, j, k)[iDim];
					}
					coord_i[Idx(ni - 2, j, k)][iDim] = MidNodeInter4thLeft1(inter_temp);
					coord_i[Idx(ni - 3, j, k)][iDim] = MidNodeInter4thLeft2(inter_temp);
					coord_i[Idx(ni - 4, j, k)][iDim] = MidNodeInter4th(inter_temp);
				}
			}
		}
		// j+1/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 3; j < nj - 4; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					for (int iDim = 0; iDim < 3; ++iDim)
					{
						for (int iTemp = 0; iTemp < 6; iTemp++)
						{
							inter_temp[iTemp] = node->GetCoord(i, j + iTemp - 2, k)[iDim];
						}
						coord_j[Idx(i, j, k)][iDim] = MidNodeInter6th(inter_temp);
					}
				}
			}
		}
		// j=1/2,3/2,5/2;nj-3/2,nj-5/2,nj-7/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				for (int iDim = 0; iDim < 3; ++iDim)
				{
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = node->GetCoord(i, 1 + iTemp, k)[iDim];
					}
					coord_j[Idx(i, 0, k)][iDim] = MidNodeInter4thRight1(inter_temp);
					coord_j[Idx(i, 1, k)][iDim] = MidNodeInter4thRight2(inter_temp);
					coord_j[Idx(i, 2, k)][iDim] = MidNodeInter4th(inter_temp);
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = node->GetCoord(i, nj - 5 + iTemp, k)[iDim];
					}
					coord_j[Idx(i, nj - 2, k)][iDim] = MidNodeInter4thLeft1(inter_temp);
					coord_j[Idx(i, nj - 3, k)][iDim] = MidNodeInter4thLeft2(inter_temp);
					coord_j[Idx(i, nj - 4, k)][iDim] = MidNodeInter4th(inter_temp);
				}
			}
		}
		// k+1/2
		for (int k = 3; k < nk - 4; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					for (int iDim = 0; iDim < 3; ++iDim)
					{
						for (int iTemp = 0; iTemp < 6; iTemp++)
						{
							inter_temp[iTemp] = node->GetCoord(i, j, k + iTemp - 2)[iDim];
						}
						coord_k[Idx(i, j, k)][iDim] = MidNodeInter6th(inter_temp);
					}
				}
			}
		}
		// k=1/2,3/2,5/2;nk-3/2,nk-5/2,nk-7/2
		for (int j = 1; j < nj - 1; ++j)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				for (int iDim = 0; iDim < 3; ++iDim)
				{
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = node->GetCoord(i, j, 1 + iTemp)[iDim];
					}
					coord_k[Idx(i, j, 0)][iDim] = MidNodeInter4thRight1(inter_temp);
					coord_k[Idx(i, j, 1)][iDim] = MidNodeInter4thRight2(inter_temp);
					coord_k[Idx(i, j, 2)][iDim] = MidNodeInter4th(inter_temp);
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = node->GetCoord(i, j, nk - 5 + iTemp)[iDim];
					}
					coord_k[Idx(i, j, nk - 2)][iDim] = MidNodeInter4thLeft1(inter_temp);
					coord_k[Idx(i, j, nk - 3)][iDim] = MidNodeInter4thLeft2(inter_temp);
					coord_k[Idx(i, j, nk - 4)][iDim] = MidNodeInter4th(inter_temp);
				}
			}
		}
		// 第二步：根据半点坐标计算整点逆变换度量系数
		double diff_temp[6];
		// i direction
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 3; i < ni - 3; ++i)
				{
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_i[Idx(i + iTemp - 3, j, k)][0];
					}
					auto coef_x = coef->GetX(Idx(i, j, k));
					coef_x[0] = NodeDifferece6th(diff_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_i[Idx(i + iTemp - 3, j, k)][1];
					}
					auto coef_y = coef->GetY(Idx(i, j, k));
					coef_y[0] = NodeDifferece6th(diff_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_i[Idx(i + iTemp - 3, j, k)][2];
					}
					auto coef_z = coef->GetZ(Idx(i, j, k));
					coef_z[0] = NodeDifferece6th(diff_temp);
				}
			}
		}
		// j direction
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 3; j < nj - 3; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					auto coef_x = coef->GetX(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_j[Idx(i, j + iTemp - 3, k)][0];
					}
					coef_x[1] = NodeDifferece6th(diff_temp);
					auto coef_y = coef->GetY(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_j[Idx(i, j + iTemp - 3, k)][1];
					}
					coef_y[1] = NodeDifferece6th(diff_temp);
					auto coef_z = coef->GetZ(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_j[Idx(i, j + iTemp - 3, k)][2];
					}
					coef_z[1] = NodeDifferece6th(diff_temp);
				}
			}
		}
		// k direction
		for (int k = 3; k < nk - 3; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					auto coef_x = coef->GetX(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_k[Idx(i, j, k + iTemp - 3)][0];
					}
					coef_x[2] = NodeDifferece6th(diff_temp);
					auto coef_y = coef->GetY(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_k[Idx(i, j, k + iTemp - 3)][1];
					}
					coef_y[2] = NodeDifferece6th(diff_temp);
					auto coef_z = coef->GetZ(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_k[Idx(i, j, k + iTemp - 3)][2];
					}
					coef_z[2] = NodeDifferece6th(diff_temp);
					if (grid->GetDim() == TWO_DIM)
					{
						coef_x[2] = 0.0;
						coef_y[2] = 0.0;
						coef_z[2] = 1.0;
					}
				}
			}
		}
		// 计算边界整点的度量系数
		// i=1,2;ni-3,ni-2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_i[Idx(0 + iTemp, j, k)][0];
				}
				coef->GetX(Idx(1, j, k))[0] = NodeDifferece4thRight(diff_temp);
				coef->GetX(Idx(2, j, k))[0] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_i[Idx(0 + iTemp, j, k)][1];
				}
				coef->GetY(Idx(1, j, k))[0] = NodeDifferece4thRight(diff_temp);
				coef->GetY(Idx(2, j, k))[0] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_i[Idx(0 + iTemp, j, k)][2];
				}
				coef->GetZ(Idx(1, j, k))[0] = NodeDifferece4thRight(diff_temp);
				coef->GetZ(Idx(2, j, k))[0] = NodeDifferece4th(diff_temp);

				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_i[Idx(ni - 2 - iTemp, j, k)][0];
				}
				coef->GetX(Idx(ni - 2, j, k))[0] = NodeDifferece4thLeft(diff_temp);
				coef->GetX(Idx(ni - 3, j, k))[0] = -NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_i[Idx(ni - 2 - iTemp, j, k)][1];
				}
				coef->GetY(Idx(ni - 2, j, k))[0] = NodeDifferece4thLeft(diff_temp);
				coef->GetY(Idx(ni - 3, j, k))[0] = -NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_i[Idx(ni - 2 - iTemp, j, k)][2];
				}
				coef->GetZ(Idx(ni - 2, j, k))[0] = NodeDifferece4thLeft(diff_temp);
				coef->GetZ(Idx(ni - 3, j, k))[0] = -NodeDifferece4th(diff_temp);
			}
		}
		// j=1,2;nj-3,nj-2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_j[Idx(i, 0 + iTemp, k)][0];
				}
				coef->GetX(Idx(i, 1, k))[1] = NodeDifferece4thRight(diff_temp);
				coef->GetX(Idx(i, 2, k))[1] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_j[Idx(i, 0 + iTemp, k)][1];
				}
				coef->GetY(Idx(i, 1, k))[1] = NodeDifferece4thRight(diff_temp);
				coef->GetY(Idx(i, 2, k))[1] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_j[Idx(i, 0 + iTemp, k)][2];
				}
				coef->GetZ(Idx(i, 1, k))[1] = NodeDifferece4thRight(diff_temp);
				coef->GetZ(Idx(i, 2, k))[1] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_j[Idx(i, nj - 2 - iTemp, k)][0];
				}
				coef->GetX(Idx(i, nj - 2, k))[1] = NodeDifferece4thLeft(diff_temp);
				coef->GetX(Idx(i, nj - 3, k))[1] = -NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_j[Idx(i, nj - 2 - iTemp, k)][1];
				}
				coef->GetY(Idx(i, nj - 2, k))[1] = NodeDifferece4thLeft(diff_temp);
				coef->GetY(Idx(i, nj - 3, k))[1] = -NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_j[Idx(i, nj - 2 - iTemp, k)][2];
				}
				coef->GetZ(Idx(i, nj - 2, k))[1] = NodeDifferece4thLeft(diff_temp);
				coef->GetZ(Idx(i, nj - 3, k))[1] = -NodeDifferece4th(diff_temp);
			}
		}
		// k=1,2;nk-3,nk-2
		for (int j = 1; j < nj - 1; ++j)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_k[Idx(i, j, 0 + iTemp)][0];
				}
				coef->GetX(Idx(i, j, 1))[2] = NodeDifferece4thRight(diff_temp);
				coef->GetX(Idx(i, j, 2))[2] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_k[Idx(i, j, 0 + iTemp)][1];
				}
				coef->GetY(Idx(i, j, 1))[2] = NodeDifferece4thRight(diff_temp);
				coef->GetY(Idx(i, j, 2))[2] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_k[Idx(i, j, 0 + iTemp)][2];
				}
				coef->GetZ(Idx(i, j, 1))[2] = NodeDifferece4thRight(diff_temp);
				coef->GetZ(Idx(i, j, 2))[2] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_k[Idx(i, j, nk - 2 - iTemp)][0];
				}
				coef->GetX(Idx(i, j, nk - 2))[2] = NodeDifferece4thLeft(diff_temp);
				coef->GetX(Idx(i, j, nk - 3))[2] = -NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_k[Idx(i, j, nk - 2 - iTemp)][1];
				}
				coef->GetY(Idx(i, j, nk - 2))[2] = NodeDifferece4thLeft(diff_temp);
				coef->GetY(Idx(i, j, nk - 3))[2] = -NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_k[Idx(i, j, nk - 2 - iTemp)][2];
				}
				coef->GetZ(Idx(i, j, nk - 2))[2] = NodeDifferece4thLeft(diff_temp);
				coef->GetZ(Idx(i, j, nk - 3))[2] = -NodeDifferece4th(diff_temp);
				if (grid->GetDim() == TWO_DIM)
				{
					coef->GetZ(Idx(i, j, 1))[0] = coef->GetZ(Idx(i, j, 1))[1] = 0.0;
					coef->GetZ(Idx(i, j, 2))[0] = coef->GetZ(Idx(i, j, 2))[1] = 0.0;
					coef->GetZ(Idx(i, j, nk - 2))[0] = coef->GetZ(Idx(i, j, nk - 2))[1] = 0.0;
					coef->GetZ(Idx(i, j, nk - 3))[0] = coef->GetZ(Idx(i, j, nk - 3))[1] = 0.0;
					coef->GetZ(Idx(i, j, 1))[2] = 1.0;
					coef->GetZ(Idx(i, j, 2))[2] = 1.0;
					coef->GetZ(Idx(i, j, nk - 2))[2] = 1.0;
					coef->GetZ(Idx(i, j, nk - 3))[2] = 1.0;
				}
			}
		}

		// 第三步：计算半点处逆变换度量系数和坐标乘积
		//  求出xi方向的临时变量：逆变换度量系数乘以坐标
		struct TempXi
		{
			double y_zeta_z, z_zeta_x, x_zeta_y, y_eta_z, z_eta_x, x_eta_y;
		};
		dynamic_array<TempXi> temp_xi(ni * nj * nk);
		// i+1/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 3; i < ni - 4; ++i)
				{
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetY(Idx(i + iTemp - 2, j, k))[2] * node->GetCoord(i + iTemp - 2, j, k)[2];
					}
					temp_xi[Idx(i, j, k)].y_zeta_z = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetZ(Idx(i + iTemp - 2, j, k))[2] * node->GetCoord(i + iTemp - 2, j, k)[0];
					}
					temp_xi[Idx(i, j, k)].z_zeta_x = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetX(Idx(i + iTemp - 2, j, k))[2] * node->GetCoord(i + iTemp - 2, j, k)[1];
					}
					temp_xi[Idx(i, j, k)].x_zeta_y = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetY(Idx(i + iTemp - 2, j, k))[1] * node->GetCoord(i + iTemp - 2, j, k)[2];
					}
					temp_xi[Idx(i, j, k)].y_eta_z = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetZ(Idx(i + iTemp - 2, j, k))[1] * node->GetCoord(i + iTemp - 2, j, k)[0];
					}
					temp_xi[Idx(i, j, k)].z_eta_x = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetX(Idx(i + iTemp - 2, j, k))[1] * node->GetCoord(i + iTemp - 2, j, k)[1];
					}
					temp_xi[Idx(i, j, k)].x_eta_y = MidNodeInter6th(inter_temp);
				}
			}
		}
		// i=1/2,3/2,5/2;ni-3/2,ni-5/2,ni-7/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				// i=1/2，3/2，5/2
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(1 + iTemp, j, k))[2] * node->GetCoord(1 + iTemp, j, k)[2];
				}
				temp_xi[Idx(0, j, k)].y_zeta_z = MidNodeInter4thRight1(inter_temp);
				temp_xi[Idx(1, j, k)].y_zeta_z = MidNodeInter4thRight2(inter_temp);
				temp_xi[Idx(2, j, k)].y_zeta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(1 + iTemp, j, k))[2] * node->GetCoord(1 + iTemp, j, k)[0];
				}
				temp_xi[Idx(0, j, k)].z_zeta_x = MidNodeInter4thRight1(inter_temp);
				temp_xi[Idx(1, j, k)].z_zeta_x = MidNodeInter4thRight2(inter_temp);
				temp_xi[Idx(2, j, k)].z_zeta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(1 + iTemp, j, k))[2] * node->GetCoord(1 + iTemp, j, k)[1];
				}
				temp_xi[Idx(0, j, k)].x_zeta_y = MidNodeInter4thRight1(inter_temp);
				temp_xi[Idx(1, j, k)].x_zeta_y = MidNodeInter4thRight2(inter_temp);
				temp_xi[Idx(2, j, k)].x_zeta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(1 + iTemp, j, k))[1] * node->GetCoord(1 + iTemp, j, k)[2];
				}
				temp_xi[Idx(0, j, k)].y_eta_z = MidNodeInter4thRight1(inter_temp);
				temp_xi[Idx(1, j, k)].y_eta_z = MidNodeInter4thRight2(inter_temp);
				temp_xi[Idx(2, j, k)].y_eta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(1 + iTemp, j, k))[1] * node->GetCoord(1 + iTemp, j, k)[0];
				}
				temp_xi[Idx(0, j, k)].z_eta_x = MidNodeInter4thRight1(inter_temp);
				temp_xi[Idx(1, j, k)].z_eta_x = MidNodeInter4thRight2(inter_temp);
				temp_xi[Idx(2, j, k)].z_eta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(1 + iTemp, j, k))[1] * node->GetCoord(1 + iTemp, j, k)[1];
				}
				temp_xi[Idx(0, j, k)].x_eta_y = MidNodeInter4thRight1(inter_temp);
				temp_xi[Idx(1, j, k)].x_eta_y = MidNodeInter4thRight2(inter_temp);
				temp_xi[Idx(2, j, k)].x_eta_y = MidNodeInter4th(inter_temp);
				// i=ni-3/2,ni-5/2,ni-7/2
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(ni - 5 + iTemp, j, k))[2] * node->GetCoord(ni - 5 + iTemp, j, k)[2];
				}
				temp_xi[Idx(ni - 2, j, k)].y_zeta_z = MidNodeInter4thLeft1(inter_temp);
				temp_xi[Idx(ni - 3, j, k)].y_zeta_z = MidNodeInter4thLeft2(inter_temp);
				temp_xi[Idx(ni - 4, j, k)].y_zeta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(ni - 5 + iTemp, j, k))[2] * node->GetCoord(ni - 5 + iTemp, j, k)[0];
				}
				temp_xi[Idx(ni - 2, j, k)].z_zeta_x = MidNodeInter4thLeft1(inter_temp);
				temp_xi[Idx(ni - 3, j, k)].z_zeta_x = MidNodeInter4thLeft2(inter_temp);
				temp_xi[Idx(ni - 4, j, k)].z_zeta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(ni - 5 + iTemp, j, k))[2] * node->GetCoord(ni - 5 + iTemp, j, k)[1];
				}
				temp_xi[Idx(ni - 2, j, k)].x_zeta_y = MidNodeInter4thLeft1(inter_temp);
				temp_xi[Idx(ni - 3, j, k)].x_zeta_y = MidNodeInter4thLeft2(inter_temp);
				temp_xi[Idx(ni - 4, j, k)].x_zeta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(ni - 5 + iTemp, j, k))[1] * node->GetCoord(ni - 5 + iTemp, j, k)[2];
				}
				temp_xi[Idx(ni - 2, j, k)].y_eta_z = MidNodeInter4thLeft1(inter_temp);
				temp_xi[Idx(ni - 3, j, k)].y_eta_z = MidNodeInter4thLeft2(inter_temp);
				temp_xi[Idx(ni - 4, j, k)].y_eta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(ni - 5 + iTemp, j, k))[1] * node->GetCoord(ni - 5 + iTemp, j, k)[0];
				}
				temp_xi[Idx(ni - 2, j, k)].z_eta_x = MidNodeInter4thLeft1(inter_temp);
				temp_xi[Idx(ni - 3, j, k)].z_eta_x = MidNodeInter4thLeft2(inter_temp);
				temp_xi[Idx(ni - 4, j, k)].z_eta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(ni - 5 + iTemp, j, k))[1] * node->GetCoord(ni - 5 + iTemp, j, k)[1];
				}
				temp_xi[Idx(ni - 2, j, k)].x_eta_y = MidNodeInter4thLeft1(inter_temp);
				temp_xi[Idx(ni - 3, j, k)].x_eta_y = MidNodeInter4thLeft2(inter_temp);
				temp_xi[Idx(ni - 4, j, k)].x_eta_y = MidNodeInter4th(inter_temp);
			}
		}
		// 求出eta方向的临时变量：逆变换度量系数乘以坐标
		struct TempEta
		{
			double y_zeta_z, z_zeta_x, x_zeta_y, y_xi_z, z_xi_x, x_xi_y;
		};
		dynamic_array<TempEta> temp_eta(ni * nj * nk);
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 3; j < nj - 4; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetY(Idx(i, j + iTemp - 2, k))[2] * node->GetCoord(i, j + iTemp - 2, k)[2];
					}
					temp_eta[Idx(i, j, k)].y_zeta_z = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetZ(Idx(i, j + iTemp - 2, k))[2] * node->GetCoord(i, j + iTemp - 2, k)[0];
					}
					temp_eta[Idx(i, j, k)].z_zeta_x = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetX(Idx(i, j + iTemp - 2, k))[2] * node->GetCoord(i, j + iTemp - 2, k)[1];
					}
					temp_eta[Idx(i, j, k)].x_zeta_y = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetY(Idx(i, j + iTemp - 2, k))[0] * node->GetCoord(i, j + iTemp - 2, k)[2];
					}
					temp_eta[Idx(i, j, k)].y_xi_z = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetZ(Idx(i, j + iTemp - 2, k))[0] * node->GetCoord(i, j + iTemp - 2, k)[0];
					}
					temp_eta[Idx(i, j, k)].z_xi_x = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetX(Idx(i, j + iTemp - 2, k))[0] * node->GetCoord(i, j + iTemp - 2, k)[1];
					}
					temp_eta[Idx(i, j, k)].x_xi_y = MidNodeInter6th(inter_temp);
				}
			}
		}
		// j=1/2,3/2,5/2;nj-3/2,nj-5/2,nj-7/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				// j=1/2，3/2，5/2
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, 1 + iTemp, k))[2] * node->GetCoord(i, 1 + iTemp, k)[2];
				}
				temp_eta[Idx(i, 0, k)].y_zeta_z = MidNodeInter4thRight1(inter_temp);
				temp_eta[Idx(i, 1, k)].y_zeta_z = MidNodeInter4thRight2(inter_temp);
				temp_eta[Idx(i, 2, k)].y_zeta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, 1 + iTemp, k))[2] * node->GetCoord(i, 1 + iTemp, k)[0];
				}
				temp_eta[Idx(i, 0, k)].z_zeta_x = MidNodeInter4thRight1(inter_temp);
				temp_eta[Idx(i, 1, k)].z_zeta_x = MidNodeInter4thRight2(inter_temp);
				temp_eta[Idx(i, 2, k)].z_zeta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, 1 + iTemp, k))[2] * node->GetCoord(i, 1 + iTemp, k)[1];
				}
				temp_eta[Idx(i, 0, k)].x_zeta_y = MidNodeInter4thRight1(inter_temp);
				temp_eta[Idx(i, 1, k)].x_zeta_y = MidNodeInter4thRight2(inter_temp);
				temp_eta[Idx(i, 2, k)].x_zeta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, 1 + iTemp, k))[0] * node->GetCoord(i, 1 + iTemp, k)[2];
				}
				temp_eta[Idx(i, 0, k)].y_xi_z = MidNodeInter4thRight1(inter_temp);
				temp_eta[Idx(i, 1, k)].y_xi_z = MidNodeInter4thRight2(inter_temp);
				temp_eta[Idx(i, 2, k)].y_xi_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, 1 + iTemp, k))[0] * node->GetCoord(i, 1 + iTemp, k)[0];
				}
				temp_eta[Idx(i, 0, k)].z_xi_x = MidNodeInter4thRight1(inter_temp);
				temp_eta[Idx(i, 1, k)].z_xi_x = MidNodeInter4thRight2(inter_temp);
				temp_eta[Idx(i, 2, k)].z_xi_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, 1 + iTemp, k))[0] * node->GetCoord(i, 1 + iTemp, k)[1];
				}
				temp_eta[Idx(i, 0, k)].x_xi_y = MidNodeInter4thRight1(inter_temp);
				temp_eta[Idx(i, 1, k)].x_xi_y = MidNodeInter4thRight2(inter_temp);
				temp_eta[Idx(i, 2, k)].x_xi_y = MidNodeInter4th(inter_temp);
				// j=nj-3/2,nj-5/2,nj-7/2
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, nj - 5 + iTemp, k))[2] * node->GetCoord(i, nj - 5 + iTemp, k)[2];
				}
				temp_eta[Idx(i, nj - 2, k)].y_zeta_z = MidNodeInter4thLeft1(inter_temp);
				temp_eta[Idx(i, nj - 3, k)].y_zeta_z = MidNodeInter4thLeft2(inter_temp);
				temp_eta[Idx(i, nj - 4, k)].y_zeta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, nj - 5 + iTemp, k))[2] * node->GetCoord(i, nj - 5 + iTemp, k)[0];
				}
				temp_eta[Idx(i, nj - 2, k)].z_zeta_x = MidNodeInter4thLeft1(inter_temp);
				temp_eta[Idx(i, nj - 3, k)].z_zeta_x = MidNodeInter4thLeft2(inter_temp);
				temp_eta[Idx(i, nj - 4, k)].z_zeta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, nj - 5 + iTemp, k))[2] * node->GetCoord(i, nj - 5 + iTemp, k)[1];
				}
				temp_eta[Idx(i, nj - 2, k)].x_zeta_y = MidNodeInter4thLeft1(inter_temp);
				temp_eta[Idx(i, nj - 3, k)].x_zeta_y = MidNodeInter4thLeft2(inter_temp);
				temp_eta[Idx(i, nj - 4, k)].x_zeta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, nj - 5 + iTemp, k))[0] * node->GetCoord(i, nj - 5 + iTemp, k)[2];
				}
				temp_eta[Idx(i, nj - 2, k)].y_xi_z = MidNodeInter4thLeft1(inter_temp);
				temp_eta[Idx(i, nj - 3, k)].y_xi_z = MidNodeInter4thLeft2(inter_temp);
				temp_eta[Idx(i, nj - 4, k)].y_xi_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, nj - 5 + iTemp, k))[0] * node->GetCoord(i, nj - 5 + iTemp, k)[0];
				}
				temp_eta[Idx(i, nj - 2, k)].z_xi_x = MidNodeInter4thLeft1(inter_temp);
				temp_eta[Idx(i, nj - 3, k)].z_xi_x = MidNodeInter4thLeft2(inter_temp);
				temp_eta[Idx(i, nj - 4, k)].z_xi_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, nj - 5 + iTemp, k))[0] * node->GetCoord(i, nj - 5 + iTemp, k)[1];
				}
				temp_eta[Idx(i, nj - 2, k)].x_xi_y = MidNodeInter4thLeft1(inter_temp);
				temp_eta[Idx(i, nj - 3, k)].x_xi_y = MidNodeInter4thLeft2(inter_temp);
				temp_eta[Idx(i, nj - 4, k)].x_xi_y = MidNodeInter4th(inter_temp);
			}
		}

		// 求出zeta方向的临时变量：逆变换度量系数乘以坐标
		struct TempZeta
		{
			double y_eta_z, z_eta_x, x_eta_y, y_xi_z, z_xi_x, x_xi_y;
		};
		dynamic_array<TempZeta> temp_zeta(ni * nj * nk);
		for (int k = 3; k < nk - 4; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetY(Idx(i, j, k + iTemp - 2))[1] * node->GetCoord(i, j, k + iTemp - 2)[2];
					}
					temp_zeta[Idx(i, j, k)].y_eta_z = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetZ(Idx(i, j, k + iTemp - 2))[1] * node->GetCoord(i, j, k + iTemp - 2)[0];
					}
					temp_zeta[Idx(i, j, k)].z_eta_x = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetX(Idx(i, j, k + iTemp - 2))[1] * node->GetCoord(i, j, k + iTemp - 2)[1];
					}
					temp_zeta[Idx(i, j, k)].x_eta_y = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetY(Idx(i, j, k + iTemp - 2))[0] * node->GetCoord(i, j, k + iTemp - 2)[2];
					}
					temp_zeta[Idx(i, j, k)].y_xi_z = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetZ(Idx(i, j, k + iTemp - 2))[0] * node->GetCoord(i, j, k + iTemp - 2)[0];
					}
					temp_zeta[Idx(i, j, k)].z_xi_x = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetX(Idx(i, j, k + iTemp - 2))[0] * node->GetCoord(i, j, k + iTemp - 2)[1];
					}
					temp_zeta[Idx(i, j, k)].x_xi_y = MidNodeInter6th(inter_temp);
				}
			}
		}
		// k=1/2,3/2,5/2;nk-3/2,nk-5/2,nk-7/2
		for (int j = 1; j < nj - 1; ++j)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				// k=1/2，3/2，5/2
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, j, 1 + iTemp))[1] * node->GetCoord(i, j, 1 + iTemp)[2];
				}
				temp_zeta[Idx(i, j, 0)].y_eta_z = MidNodeInter4thRight1(inter_temp);
				temp_zeta[Idx(i, j, 1)].y_eta_z = MidNodeInter4thRight2(inter_temp);
				temp_zeta[Idx(i, j, 2)].y_eta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, j, 1 + iTemp))[1] * node->GetCoord(i, j, 1 + iTemp)[0];
				}
				temp_zeta[Idx(i, j, 0)].z_eta_x = MidNodeInter4thRight1(inter_temp);
				temp_zeta[Idx(i, j, 1)].z_eta_x = MidNodeInter4thRight2(inter_temp);
				temp_zeta[Idx(i, j, 2)].z_eta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, j, 1 + iTemp))[1] * node->GetCoord(i, j, 1 + iTemp)[1];
				}
				temp_zeta[Idx(i, j, 0)].x_eta_y = MidNodeInter4thRight1(inter_temp);
				temp_zeta[Idx(i, j, 1)].x_eta_y = MidNodeInter4thRight2(inter_temp);
				temp_zeta[Idx(i, j, 2)].x_eta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, j, 1 + iTemp))[0] * node->GetCoord(i, j, 1 + iTemp)[2];
				}
				temp_zeta[Idx(i, j, 0)].y_xi_z = MidNodeInter4thRight1(inter_temp);
				temp_zeta[Idx(i, j, 1)].y_xi_z = MidNodeInter4thRight2(inter_temp);
				temp_zeta[Idx(i, j, 2)].y_xi_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, j, 1 + iTemp))[0] * node->GetCoord(i, j, 1 + iTemp)[0];
				}
				temp_zeta[Idx(i, j, 0)].z_xi_x = MidNodeInter4thRight1(inter_temp);
				temp_zeta[Idx(i, j, 1)].z_xi_x = MidNodeInter4thRight2(inter_temp);
				temp_zeta[Idx(i, j, 2)].z_xi_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, j, 1 + iTemp))[0] * node->GetCoord(i, j, 1 + iTemp)[1];
				}
				temp_zeta[Idx(i, j, 0)].x_xi_y = MidNodeInter4thRight1(inter_temp);
				temp_zeta[Idx(i, j, 1)].x_xi_y = MidNodeInter4thRight2(inter_temp);
				temp_zeta[Idx(i, j, 2)].x_xi_y = MidNodeInter4th(inter_temp);
				// k=nk-3/2,nk-5/2,nk-7/2
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, j, nk - 5 + iTemp))[1] * node->GetCoord(i, j, nk - 5 + iTemp)[2];
				}
				temp_zeta[Idx(i, j, nk - 2)].y_eta_z = MidNodeInter4thLeft1(inter_temp);
				temp_zeta[Idx(i, j, nk - 3)].y_eta_z = MidNodeInter4thLeft2(inter_temp);
				temp_zeta[Idx(i, j, nk - 4)].y_eta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, j, nk - 5 + iTemp))[1] * node->GetCoord(i, j, nk - 5 + iTemp)[0];
				}
				temp_zeta[Idx(i, j, nk - 2)].z_eta_x = MidNodeInter4thLeft1(inter_temp);
				temp_zeta[Idx(i, j, nk - 3)].z_eta_x = MidNodeInter4thLeft2(inter_temp);
				temp_zeta[Idx(i, j, nk - 4)].z_eta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, j, nk - 5 + iTemp))[1] * node->GetCoord(i, j, nk - 5 + iTemp)[1];
				}
				temp_zeta[Idx(i, j, nk - 2)].x_eta_y = MidNodeInter4thLeft1(inter_temp);
				temp_zeta[Idx(i, j, nk - 3)].x_eta_y = MidNodeInter4thLeft2(inter_temp);
				temp_zeta[Idx(i, j, nk - 4)].x_eta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, j, nk - 5 + iTemp))[0] * node->GetCoord(i, j, nk - 5 + iTemp)[2];
				}
				temp_zeta[Idx(i, j, nk - 2)].y_xi_z = MidNodeInter4thLeft1(inter_temp);
				temp_zeta[Idx(i, j, nk - 3)].y_xi_z = MidNodeInter4thLeft2(inter_temp);
				temp_zeta[Idx(i, j, nk - 4)].y_xi_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, j, nk - 5 + iTemp))[0] * node->GetCoord(i, j, nk - 5 + iTemp)[0];
				}
				temp_zeta[Idx(i, j, nk - 2)].z_xi_x = MidNodeInter4thLeft1(inter_temp);
				temp_zeta[Idx(i, j, nk - 3)].z_xi_x = MidNodeInter4thLeft2(inter_temp);
				temp_zeta[Idx(i, j, nk - 4)].z_xi_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, j, nk - 5 + iTemp))[0] * node->GetCoord(i, j, nk - 5 + iTemp)[1];
				}
				temp_zeta[Idx(i, j, nk - 2)].x_xi_y = MidNodeInter4thLeft1(inter_temp);
				temp_zeta[Idx(i, j, nk - 3)].x_xi_y = MidNodeInter4thLeft2(inter_temp);
				temp_zeta[Idx(i, j, nk - 4)].x_xi_y = MidNodeInter4th(inter_temp);
			}
		}
		// 第五步：根据半点坐标和半点逆变换度量系数使用守恒形式计算整点度量系数（CMM1)
		// i direction
		double temp[6][6];
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 3; i < ni - 3; ++i)
				{
					auto coef_xi = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, k));
					auto coef_eta = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, k));
					auto coef_zeta = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, k));
					// i direction
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						int idx = Idx(i - 3 + iTemp, j, k);
						temp[0][iTemp] = temp_xi[idx].y_zeta_z; // y_zeta*z
						temp[1][iTemp] = temp_xi[idx].z_zeta_x; // z_zeta*x
						temp[2][iTemp] = temp_xi[idx].x_zeta_y; // x_zeta*y
						temp[3][iTemp] = -temp_xi[idx].y_eta_z; // -y_eta*z
						temp[4][iTemp] = -temp_xi[idx].z_eta_x; // -z_eta*x
						temp[5][iTemp] = -temp_xi[idx].x_eta_y; // -x_eta*y
					}
					coef_eta[0] += NodeDifferece6th(temp[0]);
					coef_eta[1] += NodeDifferece6th(temp[1]);
					coef_eta[2] += NodeDifferece6th(temp[2]);
					coef_zeta[0] += NodeDifferece6th(temp[3]);
					coef_zeta[1] += NodeDifferece6th(temp[4]);
					coef_zeta[2] += NodeDifferece6th(temp[5]);
				}
			}
		}
		// j direction
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 3; j < nj - 3; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					auto coef_xi = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, k));
					auto coef_eta = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, k));
					auto coef_zeta = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						int idx = Idx(i, j - 3 + iTemp, k);
						temp[0][iTemp] = -temp_eta[idx].y_zeta_z; // -y_zeta*z
						temp[1][iTemp] = -temp_eta[idx].z_zeta_x; // -z_zeta*x
						temp[2][iTemp] = -temp_eta[idx].x_zeta_y; // -x_zeta*y
						temp[3][iTemp] = temp_eta[idx].y_xi_z;    // y_xi*z
						temp[4][iTemp] = temp_eta[idx].z_xi_x;    // z_xi*x
						temp[5][iTemp] = temp_eta[idx].x_xi_y;    // x_xi*y
					}
					coef_xi[0] += NodeDifferece6th(temp[0]);
					coef_xi[1] += NodeDifferece6th(temp[1]);
					coef_xi[2] += NodeDifferece6th(temp[2]);
					coef_zeta[0] += NodeDifferece6th(temp[3]);
					coef_zeta[1] += NodeDifferece6th(temp[4]);
					coef_zeta[2] += NodeDifferece6th(temp[5]);
				}
			}
		}
		// k direction
		for (int k = 3; k < nk - 3; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					auto coef_xi = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, k));
					auto coef_eta = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, k));
					auto coef_zeta = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						int idx = Idx(i, j, k - 3 + iTemp);
						temp[0][iTemp] = temp_zeta[idx].y_eta_z; // y_eta*z
						temp[1][iTemp] = temp_zeta[idx].z_eta_x; // z_eta*x
						temp[2][iTemp] = temp_zeta[idx].x_eta_y; // x_eta*y
						temp[3][iTemp] = -temp_zeta[idx].y_xi_z; // -y_xi*z
						temp[4][iTemp] = -temp_zeta[idx].z_xi_x; // -z_xi*x
						temp[5][iTemp] = -temp_zeta[idx].x_xi_y; // -x_xi*y
					}
					coef_xi[0] += NodeDifferece6th(temp[0]);
					coef_xi[1] += NodeDifferece6th(temp[1]);
					coef_xi[2] += NodeDifferece6th(temp[2]);
					coef_eta[0] += NodeDifferece6th(temp[3]);
					coef_eta[1] += NodeDifferece6th(temp[4]);
					coef_eta[2] += NodeDifferece6th(temp[5]);
				}
			}
		}
		// 第六步：计算边界点度量系数
		// i=1,2;ni-3,ni-2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				auto coef_xi1 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(1, j, k));
				auto coef_eta1 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(1, j, k));
				auto coef_zeta1 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(1, j, k));
				auto coef_xi2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(2, j, k));
				auto coef_eta2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(2, j, k));
				auto coef_zeta2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(2, j, k));
				for (int iTemp = 0; iTemp < 5; iTemp++)
				{
					int idx = Idx(0 + iTemp, j, k);
					temp[0][iTemp] = temp_xi[idx].y_zeta_z; // y_zeta*z
					temp[1][iTemp] = temp_xi[idx].z_zeta_x; // z_zeta*x
					temp[2][iTemp] = temp_xi[idx].x_zeta_y; // x_zeta*y
					temp[3][iTemp] = -temp_xi[idx].y_eta_z; // -y_eta*z
					temp[4][iTemp] = -temp_xi[idx].z_eta_x; // -z_eta*x
					temp[5][iTemp] = -temp_xi[idx].x_eta_y; // -x_eta*y
				}
				coef_eta1[0] += NodeDifferece4thRight(temp[0]);
				coef_eta1[1] += NodeDifferece4thRight(temp[1]);
				coef_eta1[2] += NodeDifferece4thRight(temp[2]);
				coef_zeta1[0] += NodeDifferece4thRight(temp[3]);
				coef_zeta1[1] += NodeDifferece4thRight(temp[4]);
				coef_zeta1[2] += NodeDifferece4thRight(temp[5]);
				coef_eta2[0] += NodeDifferece4th(temp[0]);
				coef_eta2[1] += NodeDifferece4th(temp[1]);
				coef_eta2[2] += NodeDifferece4th(temp[2]);
				coef_zeta2[0] += NodeDifferece4th(temp[3]);
				coef_zeta2[1] += NodeDifferece4th(temp[4]);
				coef_zeta2[2] = NodeDifferece4th(temp[5]);
				auto coef_xi_n2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(ni - 2, j, k));
				auto coef_eta_n2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(ni - 2, j, k));
				auto coef_zeta_n2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(ni - 2, j, k));
				auto coef_xi_n3 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(ni - 3, j, k));
				auto coef_eta_n3 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(ni - 3, j, k));
				auto coef_zeta_n3 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(ni - 3, j, k));
				for (int iTemp = 0; iTemp < 5; iTemp++)
				{
					int idx = Idx(ni - 2 - iTemp, j, k);
					temp[0][iTemp] = temp_xi[idx].y_zeta_z; // y_zeta*z
					temp[1][iTemp] = temp_xi[idx].z_zeta_x; // z_zeta*x
					temp[2][iTemp] = temp_xi[idx].x_zeta_y; // x_zeta*y
					temp[3][iTemp] = -temp_xi[idx].y_eta_z; // -y_eta*z
					temp[4][iTemp] = -temp_xi[idx].z_eta_x; // -z_eta*x
					temp[5][iTemp] = -temp_xi[idx].x_eta_y; // -x_eta*y
				}
				coef_eta_n2[0] += NodeDifferece4thLeft(temp[0]);
				coef_eta_n2[1] += NodeDifferece4thLeft(temp[1]);
				coef_eta_n2[2] += NodeDifferece4thLeft(temp[2]);
				coef_zeta_n2[0] += NodeDifferece4thLeft(temp[3]);
				coef_zeta_n2[1] += NodeDifferece4thLeft(temp[4]);
				coef_zeta_n2[2] += NodeDifferece4thLeft(temp[5]);
				coef_eta_n3[0] += -NodeDifferece4th(temp[0]);
				coef_eta_n3[1] += -NodeDifferece4th(temp[1]);
				coef_eta_n3[2] += -NodeDifferece4th(temp[2]);
				coef_zeta_n3[0] += -NodeDifferece4th(temp[3]);
				coef_zeta_n3[1] += -NodeDifferece4th(temp[4]);
				coef_zeta_n3[2] += -NodeDifferece4th(temp[5]);
			}
		}
		// j=1,2;nj-3,nj-2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				auto coef_xi1 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, 1, k));
				auto coef_eta1 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, 1, k));
				auto coef_zeta1 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, 1, k));
				auto coef_xi2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, 2, k));
				auto coef_eta2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, 2, k));
				auto coef_zeta2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, 2, k));
				for (int iTemp = 0; iTemp < 5; iTemp++)
				{
					int idx = Idx(i, 0 + iTemp, k);
					temp[0][iTemp] = -temp_eta[idx].y_zeta_z; // -y_zeta*z
					temp[1][iTemp] = -temp_eta[idx].z_zeta_x; // -z_zeta*x
					temp[2][iTemp] = -temp_eta[idx].x_zeta_y; // -x_zeta*y
					temp[3][iTemp] = temp_eta[idx].y_xi_z;    // y_xi*z
					temp[4][iTemp] = temp_eta[idx].z_xi_x;    // z_xi*x
					temp[5][iTemp] = temp_eta[idx].x_xi_y;    // x_xi*y
				}
				coef_xi1[0] += NodeDifferece4thRight(temp[0]);
				coef_xi1[1] += NodeDifferece4thRight(temp[1]);
				coef_xi1[2] += NodeDifferece4thRight(temp[2]);
				coef_zeta1[0] += NodeDifferece4thRight(temp[3]);
				coef_zeta1[1] += NodeDifferece4thRight(temp[4]);
				coef_zeta1[2] += NodeDifferece4thRight(temp[5]);
				coef_xi2[0] += NodeDifferece4th(temp[0]);
				coef_xi2[1] += NodeDifferece4th(temp[1]);
				coef_xi2[2] += NodeDifferece4th(temp[2]);
				coef_zeta2[0] += NodeDifferece4th(temp[3]);
				coef_zeta2[1] += NodeDifferece4th(temp[4]);
				coef_zeta2[2] += NodeDifferece4th(temp[5]);
				auto coef_xi_n2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, nj - 2, k));
				auto coef_eta_n2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, nj - 2, k));
				auto coef_zeta_n2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, nj - 2, k));
				auto coef_xi_n3 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, nj - 3, k));
				auto coef_eta_n3 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, nj - 3, k));
				auto coef_zeta_n3 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, nj - 3, k));
				for (int iTemp = 0; iTemp < 5; iTemp++)
				{
					int idx = Idx(i, nj - 2 - iTemp, k);
					temp[0][iTemp] = -temp_eta[idx].y_zeta_z; // -y_zeta*z
					temp[1][iTemp] = -temp_eta[idx].z_zeta_x; // -z_zeta*x
					temp[2][iTemp] = -temp_eta[idx].x_zeta_y; // -x_zeta*y
					temp[3][iTemp] = temp_eta[idx].y_xi_z;    // y_xi*z
					temp[4][iTemp] = temp_eta[idx].z_xi_x;    // z_xi*x
					temp[5][iTemp] = temp_eta[idx].x_xi_y;    // x_xi*y
				}
				coef_xi_n2[0] += NodeDifferece4thLeft(temp[0]);
				coef_xi_n2[1] += NodeDifferece4thLeft(temp[1]);
				coef_xi_n2[2] += NodeDifferece4thLeft(temp[2]);
				coef_zeta_n2[0] += NodeDifferece4thLeft(temp[3]);
				coef_zeta_n2[1] += NodeDifferece4thLeft(temp[4]);
				coef_zeta_n2[2] += NodeDifferece4thLeft(temp[5]);
				coef_xi_n3[0] += -NodeDifferece4th(temp[0]);
				coef_xi_n3[1] += -NodeDifferece4th(temp[1]);
				coef_xi_n3[2] += -NodeDifferece4th(temp[2]);
				coef_zeta_n3[0] += -NodeDifferece4th(temp[3]);
				coef_zeta_n3[1] += -NodeDifferece4th(temp[4]);
				coef_zeta_n3[2] += -NodeDifferece4th(temp[5]);
			}
		}
		// k=1,2;nk-3,nk-2
		for (int j = 1; j < nj - 1; ++j)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				auto coef_xi1 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, 1));
				auto coef_eta1 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, 1));
				auto coef_zeta1 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, 1));
				auto coef_xi2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, 2));
				auto coef_eta2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, 2));
				auto coef_zeta2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, 2));
				for (int iTemp = 0; iTemp < 5; iTemp++)
				{
					int idx = Idx(i, j, 0 + iTemp);
					temp[0][iTemp] = temp_zeta[idx].y_eta_z; // y_eta*z
					temp[1][iTemp] = temp_zeta[idx].z_eta_x; // z_eta*x
					temp[2][iTemp] = temp_zeta[idx].x_eta_y; // x_eta*y
					temp[3][iTemp] = -temp_zeta[idx].y_xi_z; // -y_xi*z
					temp[4][iTemp] = -temp_zeta[idx].z_xi_x; // -z_xi*x
					temp[5][iTemp] = -temp_zeta[idx].x_xi_y; // -x_xi*y
				}
				coef_xi1[0] += NodeDifferece4thRight(temp[0]);
				coef_xi1[1] += NodeDifferece4thRight(temp[1]);
				coef_xi1[2] += NodeDifferece4thRight(temp[2]);
				coef_eta1[0] += NodeDifferece4thRight(temp[3]);
				coef_eta1[1] += NodeDifferece4thRight(temp[4]);
				coef_eta1[2] += NodeDifferece4thRight(temp[5]);
				coef_xi2[0] += NodeDifferece4th(temp[0]);
				coef_xi2[1] += NodeDifferece4th(temp[1]);
				coef_xi2[2] += NodeDifferece4th(temp[2]);
				coef_eta2[0] += NodeDifferece4th(temp[3]);
				coef_eta2[1] += NodeDifferece4th(temp[4]);
				coef_eta2[2] += NodeDifferece4th(temp[5]);
				auto coef_xi_n2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, nk - 2));
				auto coef_eta_n2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, nk - 2));
				auto coef_zeta_n2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, nk - 2));
				auto coef_xi_n3 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, nk - 3));
				auto coef_eta_n3 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, nk - 3));
				auto coef_zeta_n3 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, nk - 3));
				for (int iTemp = 0; iTemp < 5; iTemp++)
				{
					int idx = Idx(i, j, nk - 2 - iTemp);
					temp[0][iTemp] = temp_zeta[idx].y_eta_z; // y_eta*z
					temp[1][iTemp] = temp_zeta[idx].z_eta_x; // z_eta*x
					temp[2][iTemp] = temp_zeta[idx].x_eta_y; // x_eta*y
					temp[3][iTemp] = -temp_zeta[idx].y_xi_z; // -y_xi*z
					temp[4][iTemp] = -temp_zeta[idx].z_xi_x; // -z_xi*x
					temp[5][iTemp] = -temp_zeta[idx].x_xi_y; // -x_xi*y
				}
				coef_xi_n2[0] += NodeDifferece4thLeft(temp[0]);
				coef_xi_n2[1] += NodeDifferece4thLeft(temp[1]);
				coef_xi_n2[2] += NodeDifferece4thLeft(temp[2]);
				coef_eta_n2[0] += NodeDifferece4thLeft(temp[3]);
				coef_eta_n2[1] += NodeDifferece4thLeft(temp[4]);
				coef_eta_n2[2] += NodeDifferece4thLeft(temp[5]);
				coef_xi_n3[0] += -NodeDifferece4th(temp[0]);
				coef_xi_n3[1] += -NodeDifferece4th(temp[1]);
				coef_xi_n3[2] += -NodeDifferece4th(temp[2]);
				coef_eta_n3[0] += -NodeDifferece4th(temp[3]);
				coef_eta_n3[1] += -NodeDifferece4th(temp[4]);
				coef_eta_n3[2] += -NodeDifferece4th(temp[5]);
			}
		}
	}
	void NSSolverStruct::CalcMetricsS2()
	{
		auto para = GetPara();
		auto metric_diff_scheme = para->GetMetricDifferenceScheme();
		if (metric_diff_scheme == MetricDifferenceScheme::SecondOrder)
		{
			CalcMetricsS2_2nd();
		}
		else if (metric_diff_scheme == MetricDifferenceScheme::SixthOrder)
		{
			CalcMetricsS2_6th();
		}
		else
		{
			Log::warn("FluxDifferenceScheme is not defined, use SecondOrder as default");
			CalcMetricsS2_2nd();
		}
	}
	void NSSolverStruct::CalcMetricsS2_2nd()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto coef = GetNodeMetrics();
		Metric* coef_mid[3] = { GetMidMetricsI(), GetMidMetricsJ(), GetMidMetricsK() };
		auto idx_proxy = GetIdxProxy();
		auto Idx = [&](int i, int j, int k) { return idx_proxy->GetIdx(i, j, k); };
		int ni = grid->GetNi();
		int nj = grid->GetNj();
		int nk = grid->GetNk();
		int idx_temp[3];
		// i+1/2,j+1/2,k+1/2处的坐标
		dynamic_array<dynamic_array<double>> coord_i(ni * nj * nk, dynamic_array<double>(3)),
			coord_j(ni * nj * nk, dynamic_array<double>(3)), coord_k(ni * nj * nk, dynamic_array<double>(3));
		// 计算之前，先把度量系数赋值为0
		for (int i = 0; i < ni; ++i)
		{
			for (int j = 0; j < nj; ++j)
			{
				for (int k = 0; k < nk; ++k)
				{
					for (int iDim = 0; iDim < 4; ++iDim)
					{
						coef->GetX(Idx(i, j, k))[iDim] = 0.0;
						coef->GetY(Idx(i, j, k))[iDim] = 0.0;
						coef->GetZ(Idx(i, j, k))[iDim] = 0.0;
						coef->GetXi(Idx(i, j, k))[iDim] = 0.0;
						coef->GetEta(Idx(i, j, k))[iDim] = 0.0;
						coef->GetZeta(Idx(i, j, k))[iDim] = 0.0;
						for (int jDim = 0; jDim < 3; ++jDim)
						{
							coef_mid[jDim]->GetX(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetY(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetZ(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetXi(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetEta(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetZeta(Idx(i, j, k))[iDim] = 0.0;
						}
					}
				}
			}
		}
		// 第一步：计算半点坐标
		for (int k = 0; k < nk - 1; ++k)
		{
			for (int j = 0; j < nj - 1; ++j)
			{
				for (int i = 0; i < ni - 1; ++i)
				{
					for (int iDim = 0; iDim < 3; ++iDim)
					{
						coord_i[Idx(i, j, k)][iDim] =
							0.5 * (node->GetCoord(i, j, k)[iDim] + node->GetCoord(i + 1, j, k)[iDim]);
						coord_j[Idx(i, j, k)][iDim] =
							0.5 * (node->GetCoord(i, j, k)[iDim] + node->GetCoord(i, j + 1, k)[iDim]);
						coord_k[Idx(i, j, k)][iDim] =
							0.5 * (node->GetCoord(i, j, k)[iDim] + node->GetCoord(i, j, k + 1)[iDim]);
					}
				}
			}
		}
		// 第二步：根据半点坐标计算整点逆变换度量系数
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					int idx = idx_proxy->GetIdx(i, j, k);
					auto coef_x = coef->GetX(idx);
					auto coef_y = coef->GetY(idx);
					auto coef_z = coef->GetZ(idx);
					for (int iDim = 0; iDim < 3; ++iDim)
					{
						// i direction
						coef_x[iDim] = coord_i[Idx(i, j, k)][iDim] - coord_i[Idx(i - 1, j, k)][iDim];
						// j direction
						coef_y[iDim] = coord_j[Idx(i, j, k)][iDim] - coord_j[Idx(i, j - 1, k)][iDim];
						// k direction
						coef_z[iDim] = coord_k[Idx(i, j, k)][iDim] - coord_k[Idx(i, j, k - 1)][iDim];
					}
					coef_x[3] = 0.0;
					coef_y[3] = 0.0;
					coef_z[3] = 0.0;
					if (grid->GetDim() == TWO_DIM)
					{
						coef_z[0] = 0.0;
						coef_z[1] = 0.0;
						coef_z[2] = 1.0;
					}
				}
			}
		}
		struct TempXi
		{
			double z_zeta_y, x_zeta_z, y_zeta_x, z_eta_y, x_eta_z, y_eta_x;
		};
		dynamic_array<TempXi> temp_xi(ni * nj * nk);
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					int idx = Idx(i, j, k);
					auto& temp = temp_xi[idx];
					temp.z_zeta_y = 0.5 * (coef->GetZ(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[1] +
						coef->GetZ(Idx(i + 1, j, k))[2] * node->GetCoord(i + 1, j, k)[1]);
					temp.x_zeta_z = 0.5 * (coef->GetX(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[2] +
						coef->GetX(Idx(i + 1, j, k))[2] * node->GetCoord(i + 1, j, k)[2]);
					temp.y_zeta_x = 0.5 * (coef->GetY(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[0] +
						coef->GetY(Idx(i + 1, j, k))[2] * node->GetCoord(i + 1, j, k)[0]);
					temp.z_eta_y = 0.5 * (coef->GetZ(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[1] +
						coef->GetZ(Idx(i + 1, j, k))[1] * node->GetCoord(i + 1, j, k)[1]);
					temp.x_eta_z = 0.5 * (coef->GetX(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[2] +
						coef->GetX(Idx(i + 1, j, k))[1] * node->GetCoord(i + 1, j, k)[2]);
					temp.y_eta_x = 0.5 * (coef->GetY(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[0] +
						coef->GetY(Idx(i + 1, j, k))[1] * node->GetCoord(i + 1, j, k)[0]);
				}
			}
		}
		// i=1/2;ni-3/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				temp_xi[Idx(0, j, k)].z_zeta_y = 0.5 * (3.0 * coef->GetZ(Idx(1, j, k))[2] * node->GetCoord(1, j, k)[1] -
					1.0 * coef->GetZ(Idx(2, j, k))[2] * node->GetCoord(2, j, k)[1]);
				temp_xi[Idx(0, j, k)].x_zeta_z = 0.5 * (3.0 * coef->GetX(Idx(1, j, k))[2] * node->GetCoord(1, j, k)[2] -
					1.0 * coef->GetX(Idx(2, j, k))[2] * node->GetCoord(2, j, k)[2]);
				temp_xi[Idx(0, j, k)].y_zeta_x = 0.5 * (3.0 * coef->GetY(Idx(1, j, k))[2] * node->GetCoord(1, j, k)[0] -
					1.0 * coef->GetY(Idx(2, j, k))[2] * node->GetCoord(2, j, k)[0]);
				temp_xi[Idx(0, j, k)].z_eta_y = 0.5 * (3.0 * coef->GetZ(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[1] -
					1.0 * coef->GetZ(Idx(2, j, k))[1] * node->GetCoord(2, j, k)[1]);
				temp_xi[Idx(0, j, k)].x_eta_z = 0.5 * (3.0 * coef->GetX(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[2] -
					1.0 * coef->GetX(Idx(2, j, k))[1] * node->GetCoord(2, j, k)[2]);
				temp_xi[Idx(0, j, k)].y_eta_x = 0.5 * (3.0 * coef->GetY(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[0] -
					1.0 * coef->GetY(Idx(2, j, k))[1] * node->GetCoord(2, j, k)[0]);
				temp_xi[Idx(ni - 2, j, k)].z_zeta_y =
					0.5 * (3.0 * coef->GetZ(Idx(ni - 2, j, k))[2] * node->GetCoord(ni - 2, j, k)[1] -
						1.0 * coef->GetZ(Idx(ni - 3, j, k))[2] * node->GetCoord(ni - 3, j, k)[1]);
				temp_xi[Idx(ni - 2, j, k)].x_zeta_z =
					0.5 * (3.0 * coef->GetX(Idx(ni - 2, j, k))[2] * node->GetCoord(ni - 2, j, k)[2] -
						1.0 * coef->GetX(Idx(ni - 3, j, k))[2] * node->GetCoord(ni - 3, j, k)[2]);
				temp_xi[Idx(ni - 2, j, k)].y_zeta_x =
					0.5 * (3.0 * coef->GetY(Idx(ni - 2, j, k))[2] * node->GetCoord(ni - 2, j, k)[0] -
						1.0 * coef->GetY(Idx(ni - 3, j, k))[2] * node->GetCoord(ni - 3, j, k)[0]);
				temp_xi[Idx(ni - 2, j, k)].z_eta_y =
					0.5 * (3.0 * coef->GetZ(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[1] -
						1.0 * coef->GetZ(Idx(ni - 3, j, k))[1] * node->GetCoord(ni - 3, j, k)[1]);
				temp_xi[Idx(ni - 2, j, k)].x_eta_z =
					0.5 * (3.0 * coef->GetX(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[2] -
						1.0 * coef->GetX(Idx(ni - 3, j, k))[1] * node->GetCoord(ni - 3, j, k)[2]);
				temp_xi[Idx(ni - 2, j, k)].y_eta_x =
					0.5 * (3.0 * coef->GetY(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[0] -
						1.0 * coef->GetY(Idx(ni - 3, j, k))[1] * node->GetCoord(ni - 3, j, k)[0]);
			}
		}
		// 求出eta方向的临时变量：逆变换度量系数乘以坐标
		struct TempEta
		{
			double z_zeta_y, x_zeta_z, y_zeta_x, z_xi_y, x_xi_z, y_xi_x;
		};
		dynamic_array<TempEta> temp_eta(ni * nj * nk);
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					int idx = Idx(i, j, k);
					auto& temp = temp_eta[idx];
					temp.z_zeta_y = 0.5 * (coef->GetZ(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[1] +
						coef->GetZ(Idx(i, j + 1, k))[2] * node->GetCoord(i, j + 1, k)[1]);
					temp.x_zeta_z = 0.5 * (coef->GetX(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[2] +
						coef->GetX(Idx(i, j + 1, k))[2] * node->GetCoord(i, j + 1, k)[2]);
					temp.y_zeta_x = 0.5 * (coef->GetY(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[0] +
						coef->GetY(Idx(i, j + 1, k))[2] * node->GetCoord(i, j + 1, k)[0]);
					temp.z_xi_y = 0.5 * (coef->GetZ(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[1] +
						coef->GetZ(Idx(i, j + 1, k))[0] * node->GetCoord(i, j + 1, k)[1]);
					temp.x_xi_z = 0.5 * (coef->GetX(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[2] +
						coef->GetX(Idx(i, j + 1, k))[0] * node->GetCoord(i, j + 1, k)[2]);
					temp.y_xi_x = 0.5 * (coef->GetY(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[0] +
						coef->GetY(Idx(i, j + 1, k))[0] * node->GetCoord(i, j + 1, k)[0]);
				}
			}
		}
		// j=1/2;nj-3/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				temp_eta[Idx(i, 0, k)].z_zeta_y = 0.5 * (3.0 * coef->GetZ(Idx(i, 1, k))[2] * node->GetCoord(i, 1, k)[1] -
					1.0 * coef->GetZ(Idx(i, 2, k))[2] * node->GetCoord(i, 2, k)[1]);
				temp_eta[Idx(i, 0, k)].x_zeta_z = 0.5 * (3.0 * coef->GetX(Idx(i, 1, k))[2] * node->GetCoord(i, 1, k)[2] -
					1.0 * coef->GetX(Idx(i, 2, k))[2] * node->GetCoord(i, 2, k)[2]);
				temp_eta[Idx(i, 0, k)].y_zeta_x = 0.5 * (3.0 * coef->GetY(Idx(i, 1, k))[2] * node->GetCoord(i, 1, k)[0] -
					1.0 * coef->GetY(Idx(i, 2, k))[2] * node->GetCoord(i, 2, k)[0]);
				temp_eta[Idx(i, 0, k)].z_xi_y = 0.5 * (3.0 * coef->GetZ(Idx(i, 1, k))[0] * node->GetCoord(i, 1, k)[1] -
					1.0 * coef->GetZ(Idx(i, 2, k))[0] * node->GetCoord(i, 2, k)[1]);
				temp_eta[Idx(i, 0, k)].x_xi_z = 0.5 * (3.0 * coef->GetX(Idx(i, 1, k))[0] * node->GetCoord(i, 1, k)[2] -
					1.0 * coef->GetX(Idx(i, 2, k))[0] * node->GetCoord(i, 2, k)[2]);
				temp_eta[Idx(i, 0, k)].y_xi_x = 0.5 * (3.0 * coef->GetY(Idx(i, 1, k))[0] * node->GetCoord(i, 1, k)[0] -
					1.0 * coef->GetY(Idx(i, 2, k))[0] * node->GetCoord(i, 2, k)[0]);
				temp_eta[Idx(i, nj - 2, k)].z_zeta_y =
					0.5 * (3.0 * coef->GetZ(Idx(i, nj - 2, k))[2] * node->GetCoord(i, nj - 2, k)[1] -
						1.0 * coef->GetZ(Idx(i, nj - 3, k))[2] * node->GetCoord(i, nj - 3, k)[1]);
				temp_eta[Idx(i, nj - 2, k)].x_zeta_z =
					0.5 * (3.0 * coef->GetX(Idx(i, nj - 2, k))[2] * node->GetCoord(i, nj - 2, k)[2] -
						1.0 * coef->GetX(Idx(i, nj - 3, k))[2] * node->GetCoord(i, nj - 3, k)[2]);
				temp_eta[Idx(i, nj - 2, k)].y_zeta_x =
					0.5 * (3.0 * coef->GetY(Idx(i, nj - 2, k))[2] * node->GetCoord(i, nj - 2, k)[0] -
						1.0 * coef->GetY(Idx(i, nj - 3, k))[2] * node->GetCoord(i, nj - 3, k)[0]);
				temp_eta[Idx(i, nj - 2, k)].z_xi_y =
					0.5 * (3.0 * coef->GetZ(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[1] -
						1.0 * coef->GetZ(Idx(i, nj - 3, k))[0] * node->GetCoord(i, nj - 3, k)[1]);
				temp_eta[Idx(i, nj - 2, k)].x_xi_z =
					0.5 * (3.0 * coef->GetX(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[2] -
						1.0 * coef->GetX(Idx(i, nj - 3, k))[0] * node->GetCoord(i, nj - 3, k)[2]);
				temp_eta[Idx(i, nj - 2, k)].y_xi_x =
					0.5 * (3.0 * coef->GetY(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[0] -
						1.0 * coef->GetY(Idx(i, nj - 3, k))[0] * node->GetCoord(i, nj - 3, k)[0]);
			}
		}
		// 求出zeta方向的临时变量：逆变换度量系数乘以坐标
		struct TempZeta
		{
			double z_eta_y, x_eta_z, y_eta_x, z_xi_y, x_xi_z, y_xi_x;
		};
		dynamic_array<TempZeta> temp_zeta(ni * nj * nk);
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					int idx = Idx(i, j, k);
					auto& temp = temp_zeta[idx];
					temp.z_eta_y = 0.5 * (coef->GetZ(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[1] +
						coef->GetZ(Idx(i, j, k + 1))[1] * node->GetCoord(i, j, k + 1)[1]);
					temp.x_eta_z = 0.5 * (coef->GetX(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[2] +
						coef->GetX(Idx(i, j, k + 1))[1] * node->GetCoord(i, j, k + 1)[2]);
					temp.y_eta_x = 0.5 * (coef->GetY(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[0] +
						coef->GetY(Idx(i, j, k + 1))[1] * node->GetCoord(i, j, k + 1)[0]);
					temp.z_xi_y = 0.5 * (coef->GetZ(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[1] +
						coef->GetZ(Idx(i, j, k + 1))[0] * node->GetCoord(i, j, k + 1)[1]);
					temp.x_xi_z = 0.5 * (coef->GetX(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[2] +
						coef->GetX(Idx(i, j, k + 1))[0] * node->GetCoord(i, j, k + 1)[2]);
					temp.y_xi_x = 0.5 * (coef->GetY(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[0] +
						coef->GetY(Idx(i, j, k + 1))[0] * node->GetCoord(i, j, k + 1)[0]);
				}
			}
		}
		// k=1/2;nk-3/2
		for (int j = 1; j < nj - 1; ++j)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				temp_zeta[Idx(i, j, 0)].z_eta_y = 0.5 * (3.0 * coef->GetZ(Idx(i, j, 1))[1] * node->GetCoord(i, j, 1)[1] -
					1.0 * coef->GetZ(Idx(i, j, 2))[1] * node->GetCoord(i, j, 2)[1]);
				temp_zeta[Idx(i, j, 0)].x_eta_z = 0.5 * (3.0 * coef->GetX(Idx(i, j, 1))[1] * node->GetCoord(i, j, 1)[2] -
					1.0 * coef->GetX(Idx(i, j, 2))[1] * node->GetCoord(i, j, 2)[2]);
				temp_zeta[Idx(i, j, 0)].y_eta_x = 0.5 * (3.0 * coef->GetY(Idx(i, j, 1))[1] * node->GetCoord(i, j, 1)[0] -
					1.0 * coef->GetY(Idx(i, j, 2))[1] * node->GetCoord(i, j, 2)[0]);
				temp_zeta[Idx(i, j, 0)].z_xi_y = 0.5 * (3.0 * coef->GetZ(Idx(i, j, 1))[0] * node->GetCoord(i, j, 1)[1] -
					1.0 * coef->GetZ(Idx(i, j, 2))[0] * node->GetCoord(i, j, 2)[1]);
				temp_zeta[Idx(i, j, 0)].x_xi_z = 0.5 * (3.0 * coef->GetX(Idx(i, j, 1))[0] * node->GetCoord(i, j, 1)[2] -
					1.0 * coef->GetX(Idx(i, j, 2))[0] * node->GetCoord(i, j, 2)[2]);
				temp_zeta[Idx(i, j, 0)].y_xi_x = 0.5 * (3.0 * coef->GetY(Idx(i, j, 1))[0] * node->GetCoord(i, j, 1)[0] -
					1.0 * coef->GetY(Idx(i, j, 2))[0] * node->GetCoord(i, j, 2)[0]);
				temp_zeta[Idx(i, j, nk - 2)].z_eta_y =
					0.5 * (3.0 * coef->GetZ(Idx(i, j, nk - 2))[1] * node->GetCoord(i, j, nk - 2)[1] -
						1.0 * coef->GetZ(Idx(i, j, nk - 3))[1] * node->GetCoord(i, j, nk - 3)[1]);
				temp_zeta[Idx(i, j, nk - 2)].x_eta_z =
					0.5 * (3.0 * coef->GetX(Idx(i, j, nk - 2))[1] * node->GetCoord(i, j, nk - 2)[2] -
						1.0 * coef->GetX(Idx(i, j, nk - 3))[1] * node->GetCoord(i, j, nk - 3)[2]);
				temp_zeta[Idx(i, j, nk - 2)].y_eta_x =
					0.5 * (3.0 * coef->GetY(Idx(i, j, nk - 2))[1] * node->GetCoord(i, j, nk - 2)[0] -
						1.0 * coef->GetY(Idx(i, j, nk - 3))[1] * node->GetCoord(i, j, nk - 3)[0]);
				temp_zeta[Idx(i, j, nk - 2)].z_xi_y =
					0.5 * (3.0 * coef->GetZ(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[1] -
						1.0 * coef->GetZ(Idx(i, j, nk - 3))[0] * node->GetCoord(i, j, nk - 3)[1]);
				temp_zeta[Idx(i, j, nk - 2)].x_xi_z =
					0.5 * (3.0 * coef->GetX(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[2] -
						1.0 * coef->GetX(Idx(i, j, nk - 3))[0] * node->GetCoord(i, j, nk - 3)[2]);
				temp_zeta[Idx(i, j, nk - 2)].y_xi_x =
					0.5 * (3.0 * coef->GetY(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[0] -
						1.0 * coef->GetY(Idx(i, j, nk - 3))[0] * node->GetCoord(i, j, nk - 3)[0]);
			}
		}
		// 第五步：根据半点坐标和半点逆变换度量系数使用守恒形式计算整点度量系数（CMM2)
		double temp[6][2];
		// 获取下标的lamda函数
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					// i direction
					auto coef_xi = coef->GetXi(Idx(i, j, k));
					auto coef_eta = coef->GetEta(Idx(i, j, k));
					auto coef_zeta = coef->GetZeta(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 2; iTemp++)
					{
						int idx = Idx(i - 1 + iTemp, j, k);
						temp[0][iTemp] = -temp_xi[idx].z_zeta_y; // -z_zeta*y
						temp[1][iTemp] = -temp_xi[idx].x_zeta_z; // -x_zeta*z
						temp[2][iTemp] = -temp_xi[idx].y_zeta_x; // -y_zeta*x
						temp[3][iTemp] = temp_xi[idx].z_eta_y;   // z_eta*y
						temp[4][iTemp] = temp_xi[idx].x_eta_z;   // x_eta*z
						temp[5][iTemp] = temp_xi[idx].y_eta_x;   // y_eta*x
					}
					coef_eta[0] += temp[0][1] - temp[0][0];
					coef_eta[1] += temp[1][1] - temp[1][0];
					coef_eta[2] += temp[2][1] - temp[2][0];
					coef_zeta[0] += temp[3][1] - temp[3][0];
					coef_zeta[1] += temp[4][1] - temp[4][0];
					coef_zeta[2] += temp[5][1] - temp[5][0];
					for (int iTemp = 0; iTemp < 2; iTemp++)
					{
						int idx = Idx(i, j - 1 + iTemp, k);
						temp[0][iTemp] = temp_eta[idx].z_zeta_y; // z_zeta*y
						temp[1][iTemp] = temp_eta[idx].x_zeta_z; // x_zeta*z
						temp[2][iTemp] = temp_eta[idx].y_zeta_x; // y_zeta*x
						temp[3][iTemp] = -temp_eta[idx].z_xi_y;  // -z_xi*y
						temp[4][iTemp] = -temp_eta[idx].x_xi_z;  // -x_xi*z
						temp[5][iTemp] = -temp_eta[idx].y_xi_x;  // -y_xi*x
					}
					coef_xi[0] += temp[0][1] - temp[0][0];
					coef_xi[1] += temp[1][1] - temp[1][0];
					coef_xi[2] += temp[2][1] - temp[2][0];
					coef_zeta[0] += temp[3][1] - temp[3][0];
					coef_zeta[1] += temp[4][1] - temp[4][0];
					coef_zeta[2] += temp[5][1] - temp[5][0];
					for (int iTemp = 0; iTemp < 2; iTemp++)
					{
						int idx = Idx(i, j, k - 1 + iTemp);
						temp[0][iTemp] = -temp_zeta[idx].z_eta_y; // -z_eta*y
						temp[1][iTemp] = -temp_zeta[idx].x_eta_z; // -x_eta*z
						temp[2][iTemp] = -temp_zeta[idx].y_eta_x; // -y_eta*x
						temp[3][iTemp] = temp_zeta[idx].z_xi_y;   // z_xi*y
						temp[4][iTemp] = temp_zeta[idx].x_xi_z;   // x_xi*z
						temp[5][iTemp] = temp_zeta[idx].y_xi_x;   // y_xi*x
					}
					coef_xi[0] += temp[0][1] - temp[0][0];
					coef_xi[1] += temp[1][1] - temp[1][0];
					coef_xi[2] += temp[2][1] - temp[2][0];
					coef_eta[0] += temp[3][1] - temp[3][0];
					coef_eta[1] += temp[4][1] - temp[4][0];
					coef_eta[2] += temp[5][1] - temp[5][0];
					if (grid->GetDim() == TWO_DIM)
					{
						coef_zeta[0] = 0.0;
						coef_zeta[1] = 0.0;
						coef_zeta[2] = 1.0;
					}
					coef_xi[3] = coef_eta[3] = coef_zeta[3] = 0.0;
				}
			}
		}
	}
	void NSSolverStruct::CalcMetricsS2_6th()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto coef = GetNodeMetrics();
		Metric* coef_mid[3] = { GetMidMetricsI(), GetMidMetricsJ(), GetMidMetricsK() };
		auto idx_proxy = GetIdxProxy();
		auto Idx = [&](int i, int j, int k) { return idx_proxy->GetIdx(i, j, k); };
		int ni = grid->GetNi();
		int nj = grid->GetNj();
		int nk = grid->GetNk();
		double inter_temp[6];
		// i+1/2,j+1/2,k+1/2处的坐标
		dynamic_array<dynamic_array<double>> coord_i(ni * nj * nk, dynamic_array<double>(3)),
			coord_j(ni * nj * nk, dynamic_array<double>(3)), coord_k(ni * nj * nk, dynamic_array<double>(3));
		// 计算之前，先把度量系数赋值为0
		for (int i = 0; i < ni; ++i)
		{
			for (int j = 0; j < nj; ++j)
			{
				for (int k = 0; k < nk; ++k)
				{
					for (int iDim = 0; iDim < 4; ++iDim)
					{
						coef->GetX(Idx(i, j, k))[iDim] = 0.0;
						coef->GetY(Idx(i, j, k))[iDim] = 0.0;
						coef->GetZ(Idx(i, j, k))[iDim] = 0.0;
						coef->GetXi(Idx(i, j, k))[iDim] = 0.0;
						coef->GetEta(Idx(i, j, k))[iDim] = 0.0;
						coef->GetZeta(Idx(i, j, k))[iDim] = 0.0;
						for (int jDim = 0; jDim < 3; ++jDim)
						{
							coef_mid[jDim]->GetX(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetY(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetZ(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetXi(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetEta(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetZeta(Idx(i, j, k))[iDim] = 0.0;
						}
					}
				}
			}
		}
		// 第一步：计算半点坐标
		// i+1/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 3; i < ni - 4; ++i)
				{
					for (int iDim = 0; iDim < 3; ++iDim)
					{
						for (int iTemp = 0; iTemp < 6; iTemp++)
						{
							inter_temp[iTemp] = node->GetCoord(i + iTemp - 2, j, k)[iDim];
						}
						coord_i[Idx(i, j, k)][iDim] = MidNodeInter6th(inter_temp);
					}
				}
			}
		}
		// i=1/2,3/2,5/2;ni-3/2,ni-5/2,ni-7/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int iDim = 0; iDim < 3; ++iDim)
				{
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = node->GetCoord(1 + iTemp, j, k)[iDim];
					}
					coord_i[Idx(0, j, k)][iDim] = MidNodeInter4thRight1(inter_temp);
					coord_i[Idx(1, j, k)][iDim] = MidNodeInter4thRight2(inter_temp);
					coord_i[Idx(2, j, k)][iDim] = MidNodeInter4th(inter_temp);
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = node->GetCoord(ni - 5 + iTemp, j, k)[iDim];
					}
					coord_i[Idx(ni - 2, j, k)][iDim] = MidNodeInter4thLeft1(inter_temp);
					coord_i[Idx(ni - 3, j, k)][iDim] = MidNodeInter4thLeft2(inter_temp);
					coord_i[Idx(ni - 4, j, k)][iDim] = MidNodeInter4th(inter_temp);
				}
			}
		}
		// j+1/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 3; j < nj - 4; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					for (int iDim = 0; iDim < 3; ++iDim)
					{
						for (int iTemp = 0; iTemp < 6; iTemp++)
						{
							inter_temp[iTemp] = node->GetCoord(i, j + iTemp - 2, k)[iDim];
						}
						coord_j[Idx(i, j, k)][iDim] = MidNodeInter6th(inter_temp);
					}
				}
			}
		}
		// j=1/2,3/2,5/2;nj-3/2,nj-5/2,nj-7/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				for (int iDim = 0; iDim < 3; ++iDim)
				{
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = node->GetCoord(i, 1 + iTemp, k)[iDim];
					}
					coord_j[Idx(i, 0, k)][iDim] = MidNodeInter4thRight1(inter_temp);
					coord_j[Idx(i, 1, k)][iDim] = MidNodeInter4thRight2(inter_temp);
					coord_j[Idx(i, 2, k)][iDim] = MidNodeInter4th(inter_temp);
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = node->GetCoord(i, nj - 5 + iTemp, k)[iDim];
					}
					coord_j[Idx(i, nj - 2, k)][iDim] = MidNodeInter4thLeft1(inter_temp);
					coord_j[Idx(i, nj - 3, k)][iDim] = MidNodeInter4thLeft2(inter_temp);
					coord_j[Idx(i, nj - 4, k)][iDim] = MidNodeInter4th(inter_temp);
				}
			}
		}
		// k+1/2
		for (int k = 3; k < nk - 4; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					for (int iDim = 0; iDim < 3; ++iDim)
					{
						for (int iTemp = 0; iTemp < 6; iTemp++)
						{
							inter_temp[iTemp] = node->GetCoord(i, j, k + iTemp - 2)[iDim];
						}
						coord_k[Idx(i, j, k)][iDim] = MidNodeInter6th(inter_temp);
					}
				}
			}
		}
		// k=1/2,3/2,5/2;nk-3/2,nk-5/2,nk-7/2
		for (int j = 1; j < nj - 1; ++j)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				for (int iDim = 0; iDim < 3; ++iDim)
				{
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = node->GetCoord(i, j, 1 + iTemp)[iDim];
					}
					coord_k[Idx(i, j, 0)][iDim] = MidNodeInter4thRight1(inter_temp);
					coord_k[Idx(i, j, 1)][iDim] = MidNodeInter4thRight2(inter_temp);
					coord_k[Idx(i, j, 2)][iDim] = MidNodeInter4th(inter_temp);
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = node->GetCoord(i, j, nk - 5 + iTemp)[iDim];
					}
					coord_k[Idx(i, j, nk - 2)][iDim] = MidNodeInter4thLeft1(inter_temp);
					coord_k[Idx(i, j, nk - 3)][iDim] = MidNodeInter4thLeft2(inter_temp);
					coord_k[Idx(i, j, nk - 4)][iDim] = MidNodeInter4th(inter_temp);
				}
			}
		}
		// 第二步：根据半点坐标计算整点逆变换度量系数
		double diff_temp[6];
		// i direction
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 3; i < ni - 3; ++i)
				{
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_i[Idx(i + iTemp - 3, j, k)][0];
					}
					auto coef_x = coef->GetX(Idx(i, j, k));
					coef_x[0] = NodeDifferece6th(diff_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_i[Idx(i + iTemp - 3, j, k)][1];
					}
					auto coef_y = coef->GetY(Idx(i, j, k));
					coef_y[0] = NodeDifferece6th(diff_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_i[Idx(i + iTemp - 3, j, k)][2];
					}
					auto coef_z = coef->GetZ(Idx(i, j, k));
					coef_z[0] = NodeDifferece6th(diff_temp);
				}
			}
		}
		// j direction
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 3; j < nj - 3; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					auto coef_x = coef->GetX(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_j[Idx(i, j + iTemp - 3, k)][0];
					}
					coef_x[1] = NodeDifferece6th(diff_temp);
					auto coef_y = coef->GetY(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_j[Idx(i, j + iTemp - 3, k)][1];
					}
					coef_y[1] = NodeDifferece6th(diff_temp);
					auto coef_z = coef->GetZ(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_j[Idx(i, j + iTemp - 3, k)][2];
					}
					coef_z[1] = NodeDifferece6th(diff_temp);
				}
			}
		}
		// k direction
		for (int k = 3; k < nk - 3; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					auto coef_x = coef->GetX(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_k[Idx(i, j, k + iTemp - 3)][0];
					}
					coef_x[2] = NodeDifferece6th(diff_temp);
					auto coef_y = coef->GetY(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_k[Idx(i, j, k + iTemp - 3)][1];
					}
					coef_y[2] = NodeDifferece6th(diff_temp);
					auto coef_z = coef->GetZ(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_k[Idx(i, j, k + iTemp - 3)][2];
					}
					coef_z[2] = NodeDifferece6th(diff_temp);
					if (grid->GetDim() == TWO_DIM)
					{
						coef_x[2] = 0.0;
						coef_y[2] = 0.0;
						coef_z[2] = 1.0;
					}
				}
			}
		}
		// 计算边界整点的度量系数
		// i=1,2;ni-3,ni-2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_i[Idx(0 + iTemp, j, k)][0];
				}
				coef->GetX(Idx(1, j, k))[0] = NodeDifferece4thRight(diff_temp);
				coef->GetX(Idx(2, j, k))[0] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_i[Idx(0 + iTemp, j, k)][1];
				}
				coef->GetY(Idx(1, j, k))[0] = NodeDifferece4thRight(diff_temp);
				coef->GetY(Idx(2, j, k))[0] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_i[Idx(0 + iTemp, j, k)][2];
				}
				coef->GetZ(Idx(1, j, k))[0] = NodeDifferece4thRight(diff_temp);
				coef->GetZ(Idx(2, j, k))[0] = NodeDifferece4th(diff_temp);

				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_i[Idx(ni - 2 - iTemp, j, k)][0];
				}
				coef->GetX(Idx(ni - 2, j, k))[0] = NodeDifferece4thLeft(diff_temp);
				coef->GetX(Idx(ni - 3, j, k))[0] = -NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_i[Idx(ni - 2 - iTemp, j, k)][1];
				}
				coef->GetY(Idx(ni - 2, j, k))[0] = NodeDifferece4thLeft(diff_temp);
				coef->GetY(Idx(ni - 3, j, k))[0] = -NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_i[Idx(ni - 2 - iTemp, j, k)][2];
				}
				coef->GetZ(Idx(ni - 2, j, k))[0] = NodeDifferece4thLeft(diff_temp);
				coef->GetZ(Idx(ni - 3, j, k))[0] = -NodeDifferece4th(diff_temp);
			}
		}
		// j=1,2;nj-3,nj-2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_j[Idx(i, 0 + iTemp, k)][0];
				}
				coef->GetX(Idx(i, 1, k))[1] = NodeDifferece4thRight(diff_temp);
				coef->GetX(Idx(i, 2, k))[1] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_j[Idx(i, 0 + iTemp, k)][1];
				}
				coef->GetY(Idx(i, 1, k))[1] = NodeDifferece4thRight(diff_temp);
				coef->GetY(Idx(i, 2, k))[1] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_j[Idx(i, 0 + iTemp, k)][2];
				}
				coef->GetZ(Idx(i, 1, k))[1] = NodeDifferece4thRight(diff_temp);
				coef->GetZ(Idx(i, 2, k))[1] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_j[Idx(i, nj - 2 - iTemp, k)][0];
				}
				coef->GetX(Idx(i, nj - 2, k))[1] = NodeDifferece4thLeft(diff_temp);
				coef->GetX(Idx(i, nj - 3, k))[1] = -NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_j[Idx(i, nj - 2 - iTemp, k)][1];
				}
				coef->GetY(Idx(i, nj - 2, k))[1] = NodeDifferece4thLeft(diff_temp);
				coef->GetY(Idx(i, nj - 3, k))[1] = -NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_j[Idx(i, nj - 2 - iTemp, k)][2];
				}
				coef->GetZ(Idx(i, nj - 2, k))[1] = NodeDifferece4thLeft(diff_temp);
				coef->GetZ(Idx(i, nj - 3, k))[1] = -NodeDifferece4th(diff_temp);
			}
		}
		// k=1,2;nk-3,nk-2
		for (int j = 1; j < nj - 1; ++j)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_k[Idx(i, j, 0 + iTemp)][0];
				}
				coef->GetX(Idx(i, j, 1))[2] = NodeDifferece4thRight(diff_temp);
				coef->GetX(Idx(i, j, 2))[2] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_k[Idx(i, j, 0 + iTemp)][1];
				}
				coef->GetY(Idx(i, j, 1))[2] = NodeDifferece4thRight(diff_temp);
				coef->GetY(Idx(i, j, 2))[2] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_k[Idx(i, j, 0 + iTemp)][2];
				}
				coef->GetZ(Idx(i, j, 1))[2] = NodeDifferece4thRight(diff_temp);
				coef->GetZ(Idx(i, j, 2))[2] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_k[Idx(i, j, nk - 2 - iTemp)][0];
				}
				coef->GetX(Idx(i, j, nk - 2))[2] = NodeDifferece4thLeft(diff_temp);
				coef->GetX(Idx(i, j, nk - 3))[2] = -NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_k[Idx(i, j, nk - 2 - iTemp)][1];
				}
				coef->GetY(Idx(i, j, nk - 2))[2] = NodeDifferece4thLeft(diff_temp);
				coef->GetY(Idx(i, j, nk - 3))[2] = -NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_k[Idx(i, j, nk - 2 - iTemp)][2];
				}
				coef->GetZ(Idx(i, j, nk - 2))[2] = NodeDifferece4thLeft(diff_temp);
				coef->GetZ(Idx(i, j, nk - 3))[2] = -NodeDifferece4th(diff_temp);
				if (grid->GetDim() == TWO_DIM)
				{
					coef->GetZ(Idx(i, j, 1))[0] = coef->GetZ(Idx(i, j, 1))[1] = 0.0;
					coef->GetZ(Idx(i, j, 2))[0] = coef->GetZ(Idx(i, j, 2))[1] = 0.0;
					coef->GetZ(Idx(i, j, nk - 2))[0] = coef->GetZ(Idx(i, j, nk - 2))[1] = 0.0;
					coef->GetZ(Idx(i, j, nk - 3))[0] = coef->GetZ(Idx(i, j, nk - 3))[1] = 0.0;
					coef->GetZ(Idx(i, j, 1))[2] = 1.0;
					coef->GetZ(Idx(i, j, 2))[2] = 1.0;
					coef->GetZ(Idx(i, j, nk - 2))[2] = 1.0;
					coef->GetZ(Idx(i, j, nk - 3))[2] = 1.0;
				}
			}
		}
		// 第三步：计算半点处逆变换度量系数和坐标乘积
		//  求出xi方向的临时变量：逆变换度量系数乘以坐标
		struct TempXi
		{
			double z_zeta_y, x_zeta_z, y_zeta_x, z_eta_y, x_eta_z, y_eta_x;
		};
		dynamic_array<TempXi> temp_xi(ni * nj * nk);
		// i+1/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 3; i < ni - 4; ++i)
				{
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetZ(Idx(i - 2 + iTemp, j, k))[2] * node->GetCoord(i - 2 + iTemp, j, k)[1];
					}
					temp_xi[Idx(i, j, k)].z_zeta_y = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetX(Idx(i - 2 + iTemp, j, k))[2] * node->GetCoord(i - 2 + iTemp, j, k)[2];
					}
					temp_xi[Idx(i, j, k)].x_zeta_z = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetY(Idx(i - 2 + iTemp, j, k))[2] * node->GetCoord(i - 2 + iTemp, j, k)[0];
					}
					temp_xi[Idx(i, j, k)].y_zeta_x = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetZ(Idx(i - 2 + iTemp, j, k))[1] * node->GetCoord(i - 2 + iTemp, j, k)[1];
					}
					temp_xi[Idx(i, j, k)].z_eta_y = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetX(Idx(i - 2 + iTemp, j, k))[1] * node->GetCoord(i - 2 + iTemp, j, k)[2];
					}
					temp_xi[Idx(i, j, k)].x_eta_z = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetY(Idx(i - 2 + iTemp, j, k))[1] * node->GetCoord(i - 2 + iTemp, j, k)[0];
					}
					temp_xi[Idx(i, j, k)].y_eta_x = MidNodeInter6th(inter_temp);
				}
			}
		}
		// i=1/2,3/2,5/2;ni-3/2,ni-5/2,ni-7/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				// i=1/2，3/2，5/2
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(1 + iTemp, j, k))[2] * node->GetCoord(1 + iTemp, j, k)[1];
				}
				temp_xi[Idx(0, j, k)].z_zeta_y = MidNodeInter4thRight1(inter_temp);
				temp_xi[Idx(1, j, k)].z_zeta_y = MidNodeInter4thRight2(inter_temp);
				temp_xi[Idx(2, j, k)].z_zeta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(1 + iTemp, j, k))[2] * node->GetCoord(1 + iTemp, j, k)[2];
				}
				temp_xi[Idx(0, j, k)].x_zeta_z = MidNodeInter4thRight1(inter_temp);
				temp_xi[Idx(1, j, k)].x_zeta_z = MidNodeInter4thRight2(inter_temp);
				temp_xi[Idx(2, j, k)].x_zeta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(1 + iTemp, j, k))[2] * node->GetCoord(1 + iTemp, j, k)[0];
				}
				temp_xi[Idx(0, j, k)].y_zeta_x = MidNodeInter4thRight1(inter_temp);
				temp_xi[Idx(1, j, k)].y_zeta_x = MidNodeInter4thRight2(inter_temp);
				temp_xi[Idx(2, j, k)].y_zeta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(1 + iTemp, j, k))[1] * node->GetCoord(1 + iTemp, j, k)[1];
				}
				temp_xi[Idx(0, j, k)].z_eta_y = MidNodeInter4thRight1(inter_temp);
				temp_xi[Idx(1, j, k)].z_eta_y = MidNodeInter4thRight2(inter_temp);
				temp_xi[Idx(2, j, k)].z_eta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(1 + iTemp, j, k))[1] * node->GetCoord(1 + iTemp, j, k)[2];
				}
				temp_xi[Idx(0, j, k)].x_eta_z = MidNodeInter4thRight1(inter_temp);
				temp_xi[Idx(1, j, k)].x_eta_z = MidNodeInter4thRight2(inter_temp);
				temp_xi[Idx(2, j, k)].x_eta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(1 + iTemp, j, k))[1] * node->GetCoord(1 + iTemp, j, k)[0];
				}
				temp_xi[Idx(0, j, k)].y_eta_x = MidNodeInter4thRight1(inter_temp);
				temp_xi[Idx(1, j, k)].y_eta_x = MidNodeInter4thRight2(inter_temp);
				temp_xi[Idx(2, j, k)].y_eta_x = MidNodeInter4th(inter_temp);
				// i=ni-3/2,ni-5/2,ni-7/2
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(ni - 5 + iTemp, j, k))[2] * node->GetCoord(ni - 5 + iTemp, j, k)[1];
				}
				temp_xi[Idx(ni - 2, j, k)].z_zeta_y = MidNodeInter4thLeft1(inter_temp);
				temp_xi[Idx(ni - 3, j, k)].z_zeta_y = MidNodeInter4thLeft2(inter_temp);
				temp_xi[Idx(ni - 4, j, k)].z_zeta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(ni - 5 + iTemp, j, k))[2] * node->GetCoord(ni - 5 + iTemp, j, k)[2];
				}
				temp_xi[Idx(ni - 2, j, k)].x_zeta_z = MidNodeInter4thLeft1(inter_temp);
				temp_xi[Idx(ni - 3, j, k)].x_zeta_z = MidNodeInter4thLeft2(inter_temp);
				temp_xi[Idx(ni - 4, j, k)].x_zeta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(ni - 5 + iTemp, j, k))[2] * node->GetCoord(ni - 5 + iTemp, j, k)[0];
				}
				temp_xi[Idx(ni - 2, j, k)].y_zeta_x = MidNodeInter4thLeft1(inter_temp);
				temp_xi[Idx(ni - 3, j, k)].y_zeta_x = MidNodeInter4thLeft2(inter_temp);
				temp_xi[Idx(ni - 4, j, k)].y_zeta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(ni - 5 + iTemp, j, k))[1] * node->GetCoord(ni - 5 + iTemp, j, k)[1];
				}
				temp_xi[Idx(ni - 2, j, k)].z_eta_y = MidNodeInter4thLeft1(inter_temp);
				temp_xi[Idx(ni - 3, j, k)].z_eta_y = MidNodeInter4thLeft2(inter_temp);
				temp_xi[Idx(ni - 4, j, k)].z_eta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(ni - 5 + iTemp, j, k))[1] * node->GetCoord(ni - 5 + iTemp, j, k)[2];
				}
				temp_xi[Idx(ni - 2, j, k)].x_eta_z = MidNodeInter4thLeft1(inter_temp);
				temp_xi[Idx(ni - 3, j, k)].x_eta_z = MidNodeInter4thLeft2(inter_temp);
				temp_xi[Idx(ni - 4, j, k)].x_eta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(ni - 5 + iTemp, j, k))[1] * node->GetCoord(ni - 5 + iTemp, j, k)[0];
				}
				temp_xi[Idx(ni - 2, j, k)].y_eta_x = MidNodeInter4thLeft1(inter_temp);
				temp_xi[Idx(ni - 3, j, k)].y_eta_x = MidNodeInter4thLeft2(inter_temp);
				temp_xi[Idx(ni - 4, j, k)].y_eta_x = MidNodeInter4th(inter_temp);
			}
		}
		// 求出eta方向的临时变量：逆变换度量系数乘以坐标
		struct TempEta
		{
			double z_zeta_y, x_zeta_z, y_zeta_x, z_xi_y, x_xi_z, y_xi_x;
		};
		dynamic_array<TempEta> temp_eta(ni * nj * nk);
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 3; j < nj - 4; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetZ(Idx(i, j - 2 + iTemp, k))[2] * node->GetCoord(i, j - 2 + iTemp, k)[1];
					}
					temp_eta[Idx(i, j, k)].z_zeta_y = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetX(Idx(i, j - 2 + iTemp, k))[2] * node->GetCoord(i, j - 2 + iTemp, k)[2];
					}
					temp_eta[Idx(i, j, k)].x_zeta_z = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetY(Idx(i, j - 2 + iTemp, k))[2] * node->GetCoord(i, j - 2 + iTemp, k)[0];
					}
					temp_eta[Idx(i, j, k)].y_zeta_x = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetZ(Idx(i, j - 2 + iTemp, k))[0] * node->GetCoord(i, j - 2 + iTemp, k)[1];
					}
					temp_eta[Idx(i, j, k)].z_xi_y = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetX(Idx(i, j - 2 + iTemp, k))[0] * node->GetCoord(i, j - 2 + iTemp, k)[2];
					}
					temp_eta[Idx(i, j, k)].x_xi_z = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetY(Idx(i, j - 2 + iTemp, k))[0] * node->GetCoord(i, j - 2 + iTemp, k)[0];
					}
					temp_eta[Idx(i, j, k)].y_xi_x = MidNodeInter6th(inter_temp);
				}
			}
		}
		// j=1/2,3/2,5/2;nj-3/2,nj-5/2,nj-7/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				// j=1/2,3/2,5/2
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, 1 + iTemp, k))[2] * node->GetCoord(i, 1 + iTemp, k)[1];
				}
				temp_eta[Idx(i, 0, k)].z_zeta_y = MidNodeInter4thRight1(inter_temp);
				temp_eta[Idx(i, 1, k)].z_zeta_y = MidNodeInter4thRight2(inter_temp);
				temp_eta[Idx(i, 2, k)].z_zeta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, 1 + iTemp, k))[2] * node->GetCoord(i, 1 + iTemp, k)[2];
				}
				temp_eta[Idx(i, 0, k)].x_zeta_z = MidNodeInter4thRight1(inter_temp);
				temp_eta[Idx(i, 1, k)].x_zeta_z = MidNodeInter4thRight2(inter_temp);
				temp_eta[Idx(i, 2, k)].x_zeta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, 1 + iTemp, k))[2] * node->GetCoord(i, 1 + iTemp, k)[0];
				}
				temp_eta[Idx(i, 0, k)].y_zeta_x = MidNodeInter4thRight1(inter_temp);
				temp_eta[Idx(i, 1, k)].y_zeta_x = MidNodeInter4thRight2(inter_temp);
				temp_eta[Idx(i, 2, k)].y_zeta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, 1 + iTemp, k))[0] * node->GetCoord(i, 1 + iTemp, k)[1];
				}
				temp_eta[Idx(i, 0, k)].z_xi_y = MidNodeInter4thRight1(inter_temp);
				temp_eta[Idx(i, 1, k)].z_xi_y = MidNodeInter4thRight2(inter_temp);
				temp_eta[Idx(i, 2, k)].z_xi_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, 1 + iTemp, k))[0] * node->GetCoord(i, 1 + iTemp, k)[2];
				}
				temp_eta[Idx(i, 0, k)].x_xi_z = MidNodeInter4thRight1(inter_temp);
				temp_eta[Idx(i, 1, k)].x_xi_z = MidNodeInter4thRight2(inter_temp);
				temp_eta[Idx(i, 2, k)].x_xi_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, 1 + iTemp, k))[0] * node->GetCoord(i, 1 + iTemp, k)[0];
				}
				temp_eta[Idx(i, 0, k)].y_xi_x = MidNodeInter4thRight1(inter_temp);
				temp_eta[Idx(i, 1, k)].y_xi_x = MidNodeInter4thRight2(inter_temp);
				temp_eta[Idx(i, 2, k)].y_xi_x = MidNodeInter4th(inter_temp);
				// j=nj-3/2,nj-5/2,nj-7/2
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, nj - 5 + iTemp, k))[2] * node->GetCoord(i, nj - 5 + iTemp, k)[1];
				}
				temp_eta[Idx(i, nj - 2, k)].z_zeta_y = MidNodeInter4thLeft1(inter_temp);
				temp_eta[Idx(i, nj - 3, k)].z_zeta_y = MidNodeInter4thLeft2(inter_temp);
				temp_eta[Idx(i, nj - 4, k)].z_zeta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, nj - 5 + iTemp, k))[2] * node->GetCoord(i, nj - 5 + iTemp, k)[2];
				}
				temp_eta[Idx(i, nj - 2, k)].x_zeta_z = MidNodeInter4thLeft1(inter_temp);
				temp_eta[Idx(i, nj - 3, k)].x_zeta_z = MidNodeInter4thLeft2(inter_temp);
				temp_eta[Idx(i, nj - 4, k)].x_zeta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, nj - 5 + iTemp, k))[2] * node->GetCoord(i, nj - 5 + iTemp, k)[0];
				}
				temp_eta[Idx(i, nj - 2, k)].y_zeta_x = MidNodeInter4thLeft1(inter_temp);
				temp_eta[Idx(i, nj - 3, k)].y_zeta_x = MidNodeInter4thLeft2(inter_temp);
				temp_eta[Idx(i, nj - 4, k)].y_zeta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, nj - 5 + iTemp, k))[0] * node->GetCoord(i, nj - 5 + iTemp, k)[1];
				}
				temp_eta[Idx(i, nj - 2, k)].z_xi_y = MidNodeInter4thLeft1(inter_temp);
				temp_eta[Idx(i, nj - 3, k)].z_xi_y = MidNodeInter4thLeft2(inter_temp);
				temp_eta[Idx(i, nj - 4, k)].z_xi_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, nj - 5 + iTemp, k))[0] * node->GetCoord(i, nj - 5 + iTemp, k)[2];
				}
				temp_eta[Idx(i, nj - 2, k)].x_xi_z = MidNodeInter4thLeft1(inter_temp);
				temp_eta[Idx(i, nj - 3, k)].x_xi_z = MidNodeInter4thLeft2(inter_temp);
				temp_eta[Idx(i, nj - 4, k)].x_xi_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, nj - 5 + iTemp, k))[0] * node->GetCoord(i, nj - 5 + iTemp, k)[0];
				}
				temp_eta[Idx(i, nj - 2, k)].y_xi_x = MidNodeInter4thLeft1(inter_temp);
				temp_eta[Idx(i, nj - 3, k)].y_xi_x = MidNodeInter4thLeft2(inter_temp);
				temp_eta[Idx(i, nj - 4, k)].y_xi_x = MidNodeInter4th(inter_temp);
			}
		}
		// 求出zeta方向的临时变量：逆变换度量系数乘以坐标
		struct TempZeta
		{
			double z_eta_y, x_eta_z, y_eta_x, z_xi_y, x_xi_z, y_xi_x;
		};
		dynamic_array<TempZeta> temp_zeta(ni * nj * nk);
		for (int k = 3; k < nk - 4; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetZ(Idx(i, j, k - 2 + iTemp))[1] * node->GetCoord(i, j, k - 2 + iTemp)[1];
					}
					temp_zeta[Idx(i, j, k)].z_eta_y = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetX(Idx(i, j, k - 2 + iTemp))[1] * node->GetCoord(i, j, k - 2 + iTemp)[2];
					}
					temp_zeta[Idx(i, j, k)].x_eta_z = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetY(Idx(i, j, k - 2 + iTemp))[1] * node->GetCoord(i, j, k - 2 + iTemp)[0];
					}
					temp_zeta[Idx(i, j, k)].y_eta_x = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetZ(Idx(i, j, k - 2 + iTemp))[0] * node->GetCoord(i, j, k - 2 + iTemp)[1];
					}
					temp_zeta[Idx(i, j, k)].z_xi_y = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetX(Idx(i, j, k - 2 + iTemp))[0] * node->GetCoord(i, j, k - 2 + iTemp)[2];
					}
					temp_zeta[Idx(i, j, k)].x_xi_z = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetY(Idx(i, j, k - 2 + iTemp))[0] * node->GetCoord(i, j, k - 2 + iTemp)[0];
					}
					temp_zeta[Idx(i, j, k)].y_xi_x = MidNodeInter6th(inter_temp);
				}
			}
		}
		// k=1/2,3/2,5/2;nk-3/2,nk-5/2,nk-7/2
		for (int j = 1; j < nj - 1; ++j)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				// k=1/2,3/2,5/2
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, j, 1 + iTemp))[1] * node->GetCoord(i, j, 1 + iTemp)[1];
				}
				temp_zeta[Idx(i, j, 0)].z_eta_y = MidNodeInter4thRight1(inter_temp);
				temp_zeta[Idx(i, j, 1)].z_eta_y = MidNodeInter4thRight2(inter_temp);
				temp_zeta[Idx(i, j, 2)].z_eta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, j, 1 + iTemp))[1] * node->GetCoord(i, j, 1 + iTemp)[2];
				}
				temp_zeta[Idx(i, j, 0)].x_eta_z = MidNodeInter4thRight1(inter_temp);
				temp_zeta[Idx(i, j, 1)].x_eta_z = MidNodeInter4thRight2(inter_temp);
				temp_zeta[Idx(i, j, 2)].x_eta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, j, 1 + iTemp))[1] * node->GetCoord(i, j, 1 + iTemp)[0];
				}
				temp_zeta[Idx(i, j, 0)].y_eta_x = MidNodeInter4thRight1(inter_temp);
				temp_zeta[Idx(i, j, 1)].y_eta_x = MidNodeInter4thRight2(inter_temp);
				temp_zeta[Idx(i, j, 2)].y_eta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, j, 1 + iTemp))[0] * node->GetCoord(i, j, 1 + iTemp)[1];
				}
				temp_zeta[Idx(i, j, 0)].z_xi_y = MidNodeInter4thRight1(inter_temp);
				temp_zeta[Idx(i, j, 1)].z_xi_y = MidNodeInter4thRight2(inter_temp);
				temp_zeta[Idx(i, j, 2)].z_xi_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, j, 1 + iTemp))[0] * node->GetCoord(i, j, 1 + iTemp)[2];
				}
				temp_zeta[Idx(i, j, 0)].x_xi_z = MidNodeInter4thRight1(inter_temp);
				temp_zeta[Idx(i, j, 1)].x_xi_z = MidNodeInter4thRight2(inter_temp);
				temp_zeta[Idx(i, j, 2)].x_xi_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, j, 1 + iTemp))[0] * node->GetCoord(i, j, 1 + iTemp)[0];
				}
				temp_zeta[Idx(i, j, 0)].y_xi_x = MidNodeInter4thRight1(inter_temp);
				temp_zeta[Idx(i, j, 1)].y_xi_x = MidNodeInter4thRight2(inter_temp);
				temp_zeta[Idx(i, j, 2)].y_xi_x = MidNodeInter4th(inter_temp);
				// k=nk-3/2,nk-5/2,nk-7/2
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, j, nk - 5 + iTemp))[1] * node->GetCoord(i, j, nk - 5 + iTemp)[1];
				}
				temp_zeta[Idx(i, j, nk - 2)].z_eta_y = MidNodeInter4thLeft1(inter_temp);
				temp_zeta[Idx(i, j, nk - 3)].z_eta_y = MidNodeInter4thLeft2(inter_temp);
				temp_zeta[Idx(i, j, nk - 4)].z_eta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, j, nk - 5 + iTemp))[1] * node->GetCoord(i, j, nk - 5 + iTemp)[2];
				}
				temp_zeta[Idx(i, j, nk - 2)].x_eta_z = MidNodeInter4thLeft1(inter_temp);
				temp_zeta[Idx(i, j, nk - 3)].x_eta_z = MidNodeInter4thLeft2(inter_temp);
				temp_zeta[Idx(i, j, nk - 4)].x_eta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, j, nk - 5 + iTemp))[1] * node->GetCoord(i, j, nk - 5 + iTemp)[0];
				}
				temp_zeta[Idx(i, j, nk - 2)].y_eta_x = MidNodeInter4thLeft1(inter_temp);
				temp_zeta[Idx(i, j, nk - 3)].y_eta_x = MidNodeInter4thLeft2(inter_temp);
				temp_zeta[Idx(i, j, nk - 4)].y_eta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, j, nk - 5 + iTemp))[0] * node->GetCoord(i, j, nk - 5 + iTemp)[1];
				}
				temp_zeta[Idx(i, j, nk - 2)].z_xi_y = MidNodeInter4thLeft1(inter_temp);
				temp_zeta[Idx(i, j, nk - 3)].z_xi_y = MidNodeInter4thLeft2(inter_temp);
				temp_zeta[Idx(i, j, nk - 4)].z_xi_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, j, nk - 5 + iTemp))[0] * node->GetCoord(i, j, nk - 5 + iTemp)[2];
				}
				temp_zeta[Idx(i, j, nk - 2)].x_xi_z = MidNodeInter4thLeft1(inter_temp);
				temp_zeta[Idx(i, j, nk - 3)].x_xi_z = MidNodeInter4thLeft2(inter_temp);
				temp_zeta[Idx(i, j, nk - 4)].x_xi_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, j, nk - 5 + iTemp))[0] * node->GetCoord(i, j, nk - 5 + iTemp)[0];
				}
				temp_zeta[Idx(i, j, nk - 2)].y_xi_x = MidNodeInter4thLeft1(inter_temp);
				temp_zeta[Idx(i, j, nk - 3)].y_xi_x = MidNodeInter4thLeft2(inter_temp);
				temp_zeta[Idx(i, j, nk - 4)].y_xi_x = MidNodeInter4th(inter_temp);
			}
		}

		// 第五步：根据半点坐标和半点逆变换度量系数使用守恒形式计算整点度量系数（CMM1)
		// i direction
		double temp[6][6];
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 3; i < ni - 3; ++i)
				{
					auto coef_xi = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, k));
					auto coef_eta = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, k));
					auto coef_zeta = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, k));
					// i direction
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						int idx = Idx(i - 3 + iTemp, j, k);
						temp[0][iTemp] = -temp_xi[idx].z_zeta_y; // -z_zeta*y
						temp[1][iTemp] = -temp_xi[idx].x_zeta_z; // -x_zeta*z
						temp[2][iTemp] = -temp_xi[idx].y_zeta_x; // -y_zeta*x
						temp[3][iTemp] = temp_xi[idx].z_eta_y;   // z_eta*y
						temp[4][iTemp] = temp_xi[idx].x_eta_z;   // x_eta*z
						temp[5][iTemp] = temp_xi[idx].y_eta_x;   // y_eta*x
					}
					coef_eta[0] += NodeDifferece6th(temp[0]);
					coef_eta[1] += NodeDifferece6th(temp[1]);
					coef_eta[2] += NodeDifferece6th(temp[2]);
					coef_zeta[0] += NodeDifferece6th(temp[3]);
					coef_zeta[1] += NodeDifferece6th(temp[4]);
					coef_zeta[2] += NodeDifferece6th(temp[5]);
				}
			}
		}
		// j direction
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 3; j < nj - 3; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					auto coef_xi = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, k));
					auto coef_eta = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, k));
					auto coef_zeta = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						int idx = Idx(i, j - 3 + iTemp, k);
						temp[0][iTemp] = temp_eta[idx].z_zeta_y; // z_zeta*y
						temp[1][iTemp] = temp_eta[idx].x_zeta_z; // x_zeta*z
						temp[2][iTemp] = temp_eta[idx].y_zeta_x; // y_zeta*x
						temp[3][iTemp] = -temp_eta[idx].z_xi_y;  // -z_xi*y
						temp[4][iTemp] = -temp_eta[idx].x_xi_z;  // -x_xi*z
						temp[5][iTemp] = -temp_eta[idx].y_xi_x;  // -y_xi*x
					}
					coef_xi[0] += NodeDifferece6th(temp[0]);
					coef_xi[1] += NodeDifferece6th(temp[1]);
					coef_xi[2] += NodeDifferece6th(temp[2]);
					coef_zeta[0] += NodeDifferece6th(temp[3]);
					coef_zeta[1] += NodeDifferece6th(temp[4]);
					coef_zeta[2] += NodeDifferece6th(temp[5]);
				}
			}
		}
		// k direction
		for (int k = 3; k < nk - 3; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					auto coef_xi = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, k));
					auto coef_eta = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, k));
					auto coef_zeta = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						int idx = Idx(i, j, k - 3 + iTemp);
						temp[0][iTemp] = -temp_zeta[idx].z_eta_y; // -z_eta*y
						temp[1][iTemp] = -temp_zeta[idx].x_eta_z; // -x_eta*z
						temp[2][iTemp] = -temp_zeta[idx].y_eta_x; // -y_eta*x
						temp[3][iTemp] = temp_zeta[idx].z_xi_y;   // z_xi*y
						temp[4][iTemp] = temp_zeta[idx].x_xi_z;   // x_xi*z
						temp[5][iTemp] = temp_zeta[idx].y_xi_x;   // y_xi*x
					}
					coef_xi[0] += NodeDifferece6th(temp[0]);
					coef_xi[1] += NodeDifferece6th(temp[1]);
					coef_xi[2] += NodeDifferece6th(temp[2]);
					coef_eta[0] += NodeDifferece6th(temp[3]);
					coef_eta[1] += NodeDifferece6th(temp[4]);
					coef_eta[2] += NodeDifferece6th(temp[5]);
				}
			}
		}
		// 第六步：计算边界点度量系数
		// i=1,2;ni-3,ni-2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				auto coef_xi1 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(1, j, k));
				auto coef_eta1 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(1, j, k));
				auto coef_zeta1 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(1, j, k));
				auto coef_xi2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(2, j, k));
				auto coef_eta2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(2, j, k));
				auto coef_zeta2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(2, j, k));
				for (int iTemp = 0; iTemp < 5; iTemp++)
				{
					int idx = Idx(0 + iTemp, j, k);
					temp[0][iTemp] = -temp_xi[idx].z_zeta_y; // -z_zeta*y
					temp[1][iTemp] = -temp_xi[idx].x_zeta_z; // -x_zeta*z
					temp[2][iTemp] = -temp_xi[idx].y_zeta_x; // -y_zeta*x
					temp[3][iTemp] = temp_xi[idx].z_eta_y;   // z_eta*y
					temp[4][iTemp] = temp_xi[idx].x_eta_z;   // x_eta*z
					temp[5][iTemp] = temp_xi[idx].y_eta_x;   // y_eta*x
				}
				coef_eta1[0] += NodeDifferece4thRight(temp[0]);
				coef_eta1[1] += NodeDifferece4thRight(temp[1]);
				coef_eta1[2] += NodeDifferece4thRight(temp[2]);
				coef_zeta1[0] += NodeDifferece4thRight(temp[3]);
				coef_zeta1[1] += NodeDifferece4thRight(temp[4]);
				coef_zeta1[2] += NodeDifferece4thRight(temp[5]);
				coef_eta2[0] += NodeDifferece4th(temp[0]);
				coef_eta2[1] += NodeDifferece4th(temp[1]);
				coef_eta2[2] += NodeDifferece4th(temp[2]);
				coef_zeta2[0] += NodeDifferece4th(temp[3]);
				coef_zeta2[1] += NodeDifferece4th(temp[4]);
				coef_zeta2[2] = NodeDifferece4th(temp[5]);
				auto coef_xi_n2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(ni - 2, j, k));
				auto coef_eta_n2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(ni - 2, j, k));
				auto coef_zeta_n2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(ni - 2, j, k));
				auto coef_xi_n3 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(ni - 3, j, k));
				auto coef_eta_n3 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(ni - 3, j, k));
				auto coef_zeta_n3 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(ni - 3, j, k));
				for (int iTemp = 0; iTemp < 5; iTemp++)
				{
					int idx = Idx(ni - 2 - iTemp, j, k);
					temp[0][iTemp] = -temp_xi[idx].z_zeta_y; // -z_zeta*y
					temp[1][iTemp] = -temp_xi[idx].x_zeta_z; // -x_zeta*z
					temp[2][iTemp] = -temp_xi[idx].y_zeta_x; // -y_zeta*x
					temp[3][iTemp] = temp_xi[idx].z_eta_y;   // z_eta*y
					temp[4][iTemp] = temp_xi[idx].x_eta_z;   // x_eta*z
					temp[5][iTemp] = temp_xi[idx].y_eta_x;   // y_eta*x
				}
				coef_eta_n2[0] += NodeDifferece4thLeft(temp[0]);
				coef_eta_n2[1] += NodeDifferece4thLeft(temp[1]);
				coef_eta_n2[2] += NodeDifferece4thLeft(temp[2]);
				coef_zeta_n2[0] += NodeDifferece4thLeft(temp[3]);
				coef_zeta_n2[1] += NodeDifferece4thLeft(temp[4]);
				coef_zeta_n2[2] += NodeDifferece4thLeft(temp[5]);
				coef_eta_n3[0] += -NodeDifferece4th(temp[0]);
				coef_eta_n3[1] += -NodeDifferece4th(temp[1]);
				coef_eta_n3[2] += -NodeDifferece4th(temp[2]);
				coef_zeta_n3[0] += -NodeDifferece4th(temp[3]);
				coef_zeta_n3[1] += -NodeDifferece4th(temp[4]);
				coef_zeta_n3[2] += -NodeDifferece4th(temp[5]);
			}
		}
		// j=1,2;nj-3,nj-2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				auto coef_xi1 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, 1, k));
				auto coef_eta1 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, 1, k));
				auto coef_zeta1 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, 1, k));
				auto coef_xi2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, 2, k));
				auto coef_eta2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, 2, k));
				auto coef_zeta2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, 2, k));
				for (int iTemp = 0; iTemp < 5; iTemp++)
				{
					int idx = Idx(i, 0 + iTemp, k);
					temp[0][iTemp] = temp_eta[idx].z_zeta_y; // z_zeta*y
					temp[1][iTemp] = temp_eta[idx].x_zeta_z; // x_zeta*z
					temp[2][iTemp] = temp_eta[idx].y_zeta_x; // y_zeta*x
					temp[3][iTemp] = -temp_eta[idx].z_xi_y;  // -z_xi*y
					temp[4][iTemp] = -temp_eta[idx].x_xi_z;  // -x_xi*z
					temp[5][iTemp] = -temp_eta[idx].y_xi_x;  // -y_xi*x
				}
				coef_xi1[0] += NodeDifferece4thRight(temp[0]);
				coef_xi1[1] += NodeDifferece4thRight(temp[1]);
				coef_xi1[2] += NodeDifferece4thRight(temp[2]);
				coef_zeta1[0] += NodeDifferece4thRight(temp[3]);
				coef_zeta1[1] += NodeDifferece4thRight(temp[4]);
				coef_zeta1[2] += NodeDifferece4thRight(temp[5]);
				coef_xi2[0] += NodeDifferece4th(temp[0]);
				coef_xi2[1] += NodeDifferece4th(temp[1]);
				coef_xi2[2] += NodeDifferece4th(temp[2]);
				coef_zeta2[0] += NodeDifferece4th(temp[3]);
				coef_zeta2[1] += NodeDifferece4th(temp[4]);
				coef_zeta2[2] += NodeDifferece4th(temp[5]);
				auto coef_xi_n2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, nj - 2, k));
				auto coef_eta_n2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, nj - 2, k));
				auto coef_zeta_n2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, nj - 2, k));
				auto coef_xi_n3 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, nj - 3, k));
				auto coef_eta_n3 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, nj - 3, k));
				auto coef_zeta_n3 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, nj - 3, k));
				for (int iTemp = 0; iTemp < 5; iTemp++)
				{
					int idx = Idx(i, nj - 2 - iTemp, k);
					temp[0][iTemp] = temp_eta[idx].z_zeta_y; // z_zeta*y
					temp[1][iTemp] = temp_eta[idx].x_zeta_z; // x_zeta*z
					temp[2][iTemp] = temp_eta[idx].y_zeta_x; // y_zeta*x
					temp[3][iTemp] = -temp_eta[idx].z_xi_y;  // -z_xi*y
					temp[4][iTemp] = -temp_eta[idx].x_xi_z;  // -x_xi*z
					temp[5][iTemp] = -temp_eta[idx].y_xi_x;  // -y_xi*x
				}
				coef_xi_n2[0] += NodeDifferece4thLeft(temp[0]);
				coef_xi_n2[1] += NodeDifferece4thLeft(temp[1]);
				coef_xi_n2[2] += NodeDifferece4thLeft(temp[2]);
				coef_zeta_n2[0] += NodeDifferece4thLeft(temp[3]);
				coef_zeta_n2[1] += NodeDifferece4thLeft(temp[4]);
				coef_zeta_n2[2] += NodeDifferece4thLeft(temp[5]);
				coef_xi_n3[0] += -NodeDifferece4th(temp[0]);
				coef_xi_n3[1] += -NodeDifferece4th(temp[1]);
				coef_xi_n3[2] += -NodeDifferece4th(temp[2]);
				coef_zeta_n3[0] += -NodeDifferece4th(temp[3]);
				coef_zeta_n3[1] += -NodeDifferece4th(temp[4]);
				coef_zeta_n3[2] += -NodeDifferece4th(temp[5]);
			}
		}
		// k=1,2;nk-3,nk-2
		for (int j = 1; j < nj - 1; ++j)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				auto coef_xi1 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, 1));
				auto coef_eta1 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, 1));
				auto coef_zeta1 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, 1));
				auto coef_xi2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, 2));
				auto coef_eta2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, 2));
				auto coef_zeta2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, 2));
				for (int iTemp = 0; iTemp < 5; iTemp++)
				{
					int idx = Idx(i, j, 0 + iTemp);
					temp[0][iTemp] = -temp_zeta[idx].z_eta_y; // -z_eta*y
					temp[1][iTemp] = -temp_zeta[idx].x_eta_z; // -x_eta*z
					temp[2][iTemp] = -temp_zeta[idx].y_eta_x; // -y_eta*x
					temp[3][iTemp] = temp_zeta[idx].z_xi_y;   // z_xi*y
					temp[4][iTemp] = temp_zeta[idx].x_xi_z;   // x_xi*z
					temp[5][iTemp] = temp_zeta[idx].y_xi_x;   // y_xi*x
				}
				coef_xi1[0] += NodeDifferece4thRight(temp[0]);
				coef_xi1[1] += NodeDifferece4thRight(temp[1]);
				coef_xi1[2] += NodeDifferece4thRight(temp[2]);
				coef_eta1[0] += NodeDifferece4thRight(temp[3]);
				coef_eta1[1] += NodeDifferece4thRight(temp[4]);
				coef_eta1[2] += NodeDifferece4thRight(temp[5]);
				coef_xi2[0] += NodeDifferece4th(temp[0]);
				coef_xi2[1] += NodeDifferece4th(temp[1]);
				coef_xi2[2] += NodeDifferece4th(temp[2]);
				coef_eta2[0] += NodeDifferece4th(temp[3]);
				coef_eta2[1] += NodeDifferece4th(temp[4]);
				coef_eta2[2] += NodeDifferece4th(temp[5]);
				auto coef_xi_n2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, nk - 2));
				auto coef_eta_n2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, nk - 2));
				auto coef_zeta_n2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, nk - 2));
				auto coef_xi_n3 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, nk - 3));
				auto coef_eta_n3 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, nk - 3));
				auto coef_zeta_n3 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, nk - 3));
				for (int iTemp = 0; iTemp < 5; iTemp++)
				{
					int idx = Idx(i, j, nk - 2 - iTemp);
					temp[0][iTemp] = -temp_zeta[idx].z_eta_y; // -z_eta*y
					temp[1][iTemp] = -temp_zeta[idx].x_eta_z; // -x_eta*z
					temp[2][iTemp] = -temp_zeta[idx].y_eta_x; // -y_eta*x
					temp[3][iTemp] = temp_zeta[idx].z_xi_y;   // z_xi*y
					temp[4][iTemp] = temp_zeta[idx].x_xi_z;   // x_xi*z
					temp[5][iTemp] = temp_zeta[idx].y_xi_x;   // y_xi*x
				}
				coef_xi_n2[0] += NodeDifferece4thLeft(temp[0]);
				coef_xi_n2[1] += NodeDifferece4thLeft(temp[1]);
				coef_xi_n2[2] += NodeDifferece4thLeft(temp[2]);
				coef_eta_n2[0] += NodeDifferece4thLeft(temp[3]);
				coef_eta_n2[1] += NodeDifferece4thLeft(temp[4]);
				coef_eta_n2[2] += NodeDifferece4thLeft(temp[5]);
				coef_xi_n3[0] += -NodeDifferece4th(temp[0]);
				coef_xi_n3[1] += -NodeDifferece4th(temp[1]);
				coef_xi_n3[2] += -NodeDifferece4th(temp[2]);
				coef_eta_n3[0] += -NodeDifferece4th(temp[3]);
				coef_eta_n3[1] += -NodeDifferece4th(temp[4]);
				coef_eta_n3[2] += -NodeDifferece4th(temp[5]);
			}
		}
	}
	void NSSolverStruct::CalcMetricsS3()
	{
		auto para = GetPara();
		auto metric_diff_scheme = para->GetMetricDifferenceScheme();
		if (metric_diff_scheme == MetricDifferenceScheme::SecondOrder)
		{
			CalcMetricsS3_2nd();
		}
		else if (metric_diff_scheme == MetricDifferenceScheme::SixthOrder)
		{
			CalcMetricsS3_6th();
		}
		else
		{
			Log::warn("FluxDifferenceScheme is not defined, use SecondOrder as default");
			CalcMetricsS3_2nd();
		}
	}
	void NSSolverStruct::CalcMetricsS3_2nd()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto coef = GetNodeMetrics();
		Metric* coef_mid[3] = { GetMidMetricsI(), GetMidMetricsJ(), GetMidMetricsK() };
		auto idx_proxy = GetIdxProxy();
		auto Idx = [&](int i, int j, int k) { return idx_proxy->GetIdx(i, j, k); };
		int ni = grid->GetNi();
		int nj = grid->GetNj();
		int nk = grid->GetNk();
		int idx_temp[3];
		// i+1/2,j+1/2,k+1/2处的坐标
		dynamic_array<dynamic_array<double>> coord_i(ni * nj * nk, dynamic_array<double>(3)),
			coord_j(ni * nj * nk, dynamic_array<double>(3)), coord_k(ni * nj * nk, dynamic_array<double>(3));
		// 计算之前，先把度量系数赋值为0
		for (int i = 0; i < ni; ++i)
		{
			for (int j = 0; j < nj; ++j)
			{
				for (int k = 0; k < nk; ++k)
				{
					for (int iDim = 0; iDim < 4; ++iDim)
					{
						coef->GetX(Idx(i, j, k))[iDim] = 0.0;
						coef->GetY(Idx(i, j, k))[iDim] = 0.0;
						coef->GetZ(Idx(i, j, k))[iDim] = 0.0;
						coef->GetXi(Idx(i, j, k))[iDim] = 0.0;
						coef->GetEta(Idx(i, j, k))[iDim] = 0.0;
						coef->GetZeta(Idx(i, j, k))[iDim] = 0.0;
						for (int jDim = 0; jDim < 3; ++jDim)
						{
							coef_mid[jDim]->GetX(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetY(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetZ(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetXi(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetEta(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetZeta(Idx(i, j, k))[iDim] = 0.0;
						}
					}
				}
			}
		}
		// 第一步：计算半点坐标
		for (int k = 0; k < nk - 1; ++k)
		{
			for (int j = 0; j < nj - 1; ++j)
			{
				for (int i = 0; i < ni - 1; ++i)
				{
					for (int iDim = 0; iDim < 3; ++iDim)
					{
						coord_i[Idx(i, j, k)][iDim] =
							0.5 * (node->GetCoord(i, j, k)[iDim] + node->GetCoord(i + 1, j, k)[iDim]);
						coord_j[Idx(i, j, k)][iDim] =
							0.5 * (node->GetCoord(i, j, k)[iDim] + node->GetCoord(i, j + 1, k)[iDim]);
						coord_k[Idx(i, j, k)][iDim] =
							0.5 * (node->GetCoord(i, j, k)[iDim] + node->GetCoord(i, j, k + 1)[iDim]);
					}
				}
			}
		}
		// 第二步：根据半点坐标计算整点逆变换度量系数
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					int idx = idx_proxy->GetIdx(i, j, k);
					auto coef_x = coef->GetX(idx);
					auto coef_y = coef->GetY(idx);
					auto coef_z = coef->GetZ(idx);
					coef_x[0] = coord_i[Idx(i, j, k)][0] - coord_i[Idx(i - 1, j, k)][0];
					coef_x[1] = coord_j[Idx(i, j, k)][0] - coord_j[Idx(i, j - 1, k)][0];
					coef_x[2] = coord_k[Idx(i, j, k)][0] - coord_k[Idx(i, j, k - 1)][0];
					coef_y[0] = coord_i[Idx(i, j, k)][1] - coord_i[Idx(i - 1, j, k)][1];
					coef_y[1] = coord_j[Idx(i, j, k)][1] - coord_j[Idx(i, j - 1, k)][1];
					coef_y[2] = coord_k[Idx(i, j, k)][1] - coord_k[Idx(i, j, k - 1)][1];
					coef_z[0] = coord_i[Idx(i, j, k)][2] - coord_i[Idx(i - 1, j, k)][2];
					coef_z[1] = coord_j[Idx(i, j, k)][2] - coord_j[Idx(i, j - 1, k)][2];
					coef_z[2] = coord_k[Idx(i, j, k)][2] - coord_k[Idx(i, j, k - 1)][2];
					coef_x[3] = 0.0;
					coef_y[3] = 0.0;
					coef_z[3] = 0.0;
					if (grid->GetDim() == TWO_DIM)
					{
						coef_z[0] = 0.0;
						coef_z[1] = 0.0;
						coef_z[2] = 1.0;
					}
				}
			}
		}
		// 第三步：计算半点处逆变换度量系数和坐标乘积
		//  求出xi方向的临时变量：逆变换度量系数乘以坐标
		struct TempXi
		{
			double y_zeta_z, z_zeta_x, x_zeta_y, y_eta_z, z_eta_x, x_eta_y, z_zeta_y, x_zeta_z, y_zeta_x, z_eta_y, x_eta_z,
				y_eta_x;
		};
		dynamic_array<TempXi> temp_xi(ni * nj * nk);
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					int idx = Idx(i, j, k);
					auto& temp = temp_xi[idx];
					temp.y_zeta_z = 0.5 * (coef->GetY(Idx(i + 1, j, k))[2] * node->GetCoord(i + 1, j, k)[2] +
						coef->GetY(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[2]);
					temp.z_zeta_x = 0.5 * (coef->GetZ(Idx(i + 1, j, k))[2] * node->GetCoord(i + 1, j, k)[0] +
						coef->GetZ(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[0]);
					temp.x_zeta_y = 0.5 * (coef->GetX(Idx(i + 1, j, k))[2] * node->GetCoord(i + 1, j, k)[1] +
						coef->GetX(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[1]);
					temp.y_eta_z = 0.5 * (coef->GetY(Idx(i + 1, j, k))[1] * node->GetCoord(i + 1, j, k)[2] +
						coef->GetY(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[2]);
					temp.z_eta_x = 0.5 * (coef->GetZ(Idx(i + 1, j, k))[1] * node->GetCoord(i + 1, j, k)[0] +
						coef->GetZ(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[0]);
					temp.x_eta_y = 0.5 * (coef->GetX(Idx(i + 1, j, k))[1] * node->GetCoord(i + 1, j, k)[1] +
						coef->GetX(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[1]);
					temp.z_zeta_y = 0.5 * (coef->GetZ(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[1] +
						coef->GetZ(Idx(i + 1, j, k))[2] * node->GetCoord(i + 1, j, k)[1]);
					temp.x_zeta_z = 0.5 * (coef->GetX(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[2] +
						coef->GetX(Idx(i + 1, j, k))[2] * node->GetCoord(i + 1, j, k)[2]);
					temp.y_zeta_x = 0.5 * (coef->GetY(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[0] +
						coef->GetY(Idx(i + 1, j, k))[2] * node->GetCoord(i + 1, j, k)[0]);
					temp.z_eta_y = 0.5 * (coef->GetZ(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[1] +
						coef->GetZ(Idx(i + 1, j, k))[1] * node->GetCoord(i + 1, j, k)[1]);
					temp.x_eta_z = 0.5 * (coef->GetX(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[2] +
						coef->GetX(Idx(i + 1, j, k))[1] * node->GetCoord(i + 1, j, k)[2]);
					temp.y_eta_x = 0.5 * (coef->GetY(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[0] +
						coef->GetY(Idx(i + 1, j, k))[1] * node->GetCoord(i + 1, j, k)[0]);
				}
			}
		}
		// i=1/2；(ni-1)-1/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				int idx = Idx(0, j, k);
				temp_xi[idx].y_zeta_z = 0.5 * (3.0 * coef->GetY(Idx(1, j, k))[2] * node->GetCoord(1, j, k)[2] -
					1.0 * coef->GetY(Idx(2, j, k))[2] * node->GetCoord(2, j, k)[2]);
				temp_xi[idx].z_zeta_x = 0.5 * (3.0 * coef->GetZ(Idx(1, j, k))[2] * node->GetCoord(1, j, k)[0] -
					1.0 * coef->GetZ(Idx(2, j, k))[2] * node->GetCoord(2, j, k)[0]);
				temp_xi[idx].x_zeta_y = 0.5 * (3.0 * coef->GetX(Idx(1, j, k))[2] * node->GetCoord(1, j, k)[1] -
					1.0 * coef->GetX(Idx(2, j, k))[2] * node->GetCoord(2, j, k)[1]);
				temp_xi[idx].y_eta_z = 0.5 * (3.0 * coef->GetY(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[2] -
					1.0 * coef->GetY(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[2]);
				temp_xi[idx].z_eta_x = 0.5 * (3.0 * coef->GetZ(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[0] -
					1.0 * coef->GetZ(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[0]);
				temp_xi[idx].x_eta_y = 0.5 * (3.0 * coef->GetX(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[1] -
					1.0 * coef->GetX(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[1]);
				temp_xi[Idx(0, j, k)].z_zeta_y = 0.5 * (3.0 * coef->GetZ(Idx(1, j, k))[2] * node->GetCoord(1, j, k)[1] -
					1.0 * coef->GetZ(Idx(2, j, k))[2] * node->GetCoord(2, j, k)[1]);
				temp_xi[Idx(0, j, k)].x_zeta_z = 0.5 * (3.0 * coef->GetX(Idx(1, j, k))[2] * node->GetCoord(1, j, k)[2] -
					1.0 * coef->GetX(Idx(2, j, k))[2] * node->GetCoord(2, j, k)[2]);
				temp_xi[Idx(0, j, k)].y_zeta_x = 0.5 * (3.0 * coef->GetY(Idx(1, j, k))[2] * node->GetCoord(1, j, k)[0] -
					1.0 * coef->GetY(Idx(2, j, k))[2] * node->GetCoord(2, j, k)[0]);
				temp_xi[Idx(0, j, k)].z_eta_y = 0.5 * (3.0 * coef->GetZ(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[1] -
					1.0 * coef->GetZ(Idx(2, j, k))[1] * node->GetCoord(2, j, k)[1]);
				temp_xi[Idx(0, j, k)].x_eta_z = 0.5 * (3.0 * coef->GetX(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[2] -
					1.0 * coef->GetX(Idx(2, j, k))[1] * node->GetCoord(2, j, k)[2]);
				temp_xi[Idx(0, j, k)].y_eta_x = 0.5 * (3.0 * coef->GetY(Idx(1, j, k))[1] * node->GetCoord(1, j, k)[0] -
					1.0 * coef->GetY(Idx(2, j, k))[1] * node->GetCoord(2, j, k)[0]);

				idx = Idx(ni - 2, j, k);
				temp_xi[idx].y_zeta_z = 0.5 * (3.0 * coef->GetY(Idx(ni - 2, j, k))[2] * node->GetCoord(ni - 2, j, k)[2] -
					1.0 * coef->GetY(Idx(ni - 3, j, k))[2] * node->GetCoord(ni - 3, j, k)[2]);
				temp_xi[idx].z_zeta_x = 0.5 * (3.0 * coef->GetZ(Idx(ni - 2, j, k))[2] * node->GetCoord(ni - 2, j, k)[0] -
					1.0 * coef->GetZ(Idx(ni - 3, j, k))[2] * node->GetCoord(ni - 3, j, k)[0]);
				temp_xi[idx].x_zeta_y = 0.5 * (3.0 * coef->GetX(Idx(ni - 2, j, k))[2] * node->GetCoord(ni - 2, j, k)[1] -
					1.0 * coef->GetX(Idx(ni - 3, j, k))[2] * node->GetCoord(ni - 3, j, k)[1]);
				temp_xi[idx].y_eta_z = 0.5 * (3.0 * coef->GetY(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[2] -
					1.0 * coef->GetY(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[2]);
				temp_xi[idx].z_eta_x = 0.5 * (3.0 * coef->GetZ(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[0] -
					1.0 * coef->GetZ(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[0]);
				temp_xi[idx].x_eta_y = 0.5 * (3.0 * coef->GetX(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[1] -
					1.0 * coef->GetX(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[1]);
				temp_xi[Idx(ni - 2, j, k)].z_zeta_y =
					0.5 * (3.0 * coef->GetZ(Idx(ni - 2, j, k))[2] * node->GetCoord(ni - 2, j, k)[1] -
						1.0 * coef->GetZ(Idx(ni - 3, j, k))[2] * node->GetCoord(ni - 3, j, k)[1]);
				temp_xi[Idx(ni - 2, j, k)].x_zeta_z =
					0.5 * (3.0 * coef->GetX(Idx(ni - 2, j, k))[2] * node->GetCoord(ni - 2, j, k)[2] -
						1.0 * coef->GetX(Idx(ni - 3, j, k))[2] * node->GetCoord(ni - 3, j, k)[2]);
				temp_xi[Idx(ni - 2, j, k)].y_zeta_x =
					0.5 * (3.0 * coef->GetY(Idx(ni - 2, j, k))[2] * node->GetCoord(ni - 2, j, k)[0] -
						1.0 * coef->GetY(Idx(ni - 3, j, k))[2] * node->GetCoord(ni - 3, j, k)[0]);
				temp_xi[Idx(ni - 2, j, k)].z_eta_y =
					0.5 * (3.0 * coef->GetZ(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[1] -
						1.0 * coef->GetZ(Idx(ni - 3, j, k))[1] * node->GetCoord(ni - 3, j, k)[1]);
				temp_xi[Idx(ni - 2, j, k)].x_eta_z =
					0.5 * (3.0 * coef->GetX(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[2] -
						1.0 * coef->GetX(Idx(ni - 3, j, k))[1] * node->GetCoord(ni - 3, j, k)[2]);
				temp_xi[Idx(ni - 2, j, k)].y_eta_x =
					0.5 * (3.0 * coef->GetY(Idx(ni - 2, j, k))[1] * node->GetCoord(ni - 2, j, k)[0] -
						1.0 * coef->GetY(Idx(ni - 3, j, k))[1] * node->GetCoord(ni - 3, j, k)[0]);
			}
		}
		// 求出eta方向的临时变量：逆变换度量系数乘以坐标
		struct TempEta
		{
			double y_zeta_z, z_zeta_x, x_zeta_y, y_xi_z, z_xi_x, x_xi_y, z_zeta_y, x_zeta_z, y_zeta_x, z_xi_y, x_xi_z,
				y_xi_x;
		};
		dynamic_array<TempEta> temp_eta(ni * nj * nk);
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					int idx = Idx(i, j, k);
					auto& temp = temp_eta[idx];
					temp.y_zeta_z = 0.5 * (coef->GetY(Idx(i, j + 1, k))[2] * node->GetCoord(i, j + 1, k)[2] +
						coef->GetY(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[2]);
					temp.z_zeta_x = 0.5 * (coef->GetZ(Idx(i, j + 1, k))[2] * node->GetCoord(i, j + 1, k)[0] +
						coef->GetZ(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[0]);
					temp.x_zeta_y = 0.5 * (coef->GetX(Idx(i, j + 1, k))[2] * node->GetCoord(i, j + 1, k)[1] +
						coef->GetX(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[1]);
					temp.y_xi_z = 0.5 * (coef->GetY(Idx(i, j + 1, k))[0] * node->GetCoord(i, j + 1, k)[2] +
						coef->GetY(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[2]);
					temp.z_xi_x = 0.5 * (coef->GetZ(Idx(i, j + 1, k))[0] * node->GetCoord(i, j + 1, k)[0] +
						coef->GetZ(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[0]);
					temp.x_xi_y = 0.5 * (coef->GetX(Idx(i, j + 1, k))[0] * node->GetCoord(i, j + 1, k)[1] +
						coef->GetX(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[1]);
					temp.z_zeta_y = 0.5 * (coef->GetZ(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[1] +
						coef->GetZ(Idx(i, j + 1, k))[2] * node->GetCoord(i, j + 1, k)[1]);
					temp.x_zeta_z = 0.5 * (coef->GetX(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[2] +
						coef->GetX(Idx(i, j + 1, k))[2] * node->GetCoord(i, j + 1, k)[2]);
					temp.y_zeta_x = 0.5 * (coef->GetY(Idx(i, j, k))[2] * node->GetCoord(i, j, k)[0] +
						coef->GetY(Idx(i, j + 1, k))[2] * node->GetCoord(i, j + 1, k)[0]);
					temp.z_xi_y = 0.5 * (coef->GetZ(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[1] +
						coef->GetZ(Idx(i, j + 1, k))[0] * node->GetCoord(i, j + 1, k)[1]);
					temp.x_xi_z = 0.5 * (coef->GetX(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[2] +
						coef->GetX(Idx(i, j + 1, k))[0] * node->GetCoord(i, j + 1, k)[2]);
					temp.y_xi_x = 0.5 * (coef->GetY(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[0] +
						coef->GetY(Idx(i, j + 1, k))[0] * node->GetCoord(i, j + 1, k)[0]);
				}
			}
		}
		// j=1/2；(nj-1)-1/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				int idx = Idx(i, 0, k);
				temp_eta[idx].y_zeta_z = 0.5 * (3.0 * coef->GetY(Idx(i, 1, k))[2] * node->GetCoord(i, 1, k)[2] -
					1.0 * coef->GetY(Idx(i, 2, k))[2] * node->GetCoord(i, 2, k)[2]);
				temp_eta[idx].z_zeta_x = 0.5 * (3.0 * coef->GetZ(Idx(i, 1, k))[2] * node->GetCoord(i, 1, k)[0] -
					1.0 * coef->GetZ(Idx(i, 2, k))[2] * node->GetCoord(i, 2, k)[0]);
				temp_eta[idx].x_zeta_y = 0.5 * (3.0 * coef->GetX(Idx(i, 1, k))[2] * node->GetCoord(i, 1, k)[1] -
					1.0 * coef->GetX(Idx(i, 2, k))[2] * node->GetCoord(i, 2, k)[1]);
				temp_eta[idx].y_xi_z = 0.5 * (3.0 * coef->GetY(Idx(i, 1, k))[0] * node->GetCoord(i, 1, k)[2] -
					1.0 * coef->GetY(Idx(i, 2, k))[0] * node->GetCoord(i, 2, k)[2]);
				temp_eta[idx].z_xi_x = 0.5 * (3.0 * coef->GetZ(Idx(i, 1, k))[0] * node->GetCoord(i, 1, k)[0] -
					1.0 * coef->GetZ(Idx(i, 2, k))[0] * node->GetCoord(i, 2, k)[0]);
				temp_eta[idx].x_xi_y = 0.5 * (3.0 * coef->GetX(Idx(i, 1, k))[0] * node->GetCoord(i, 1, k)[1] -
					1.0 * coef->GetX(Idx(i, 2, k))[0] * node->GetCoord(i, 2, k)[1]);
				temp_eta[Idx(i, 0, k)].z_zeta_y = 0.5 * (3.0 * coef->GetZ(Idx(i, 1, k))[2] * node->GetCoord(i, 1, k)[1] -
					1.0 * coef->GetZ(Idx(i, 2, k))[2] * node->GetCoord(i, 2, k)[1]);
				temp_eta[Idx(i, 0, k)].x_zeta_z = 0.5 * (3.0 * coef->GetX(Idx(i, 1, k))[2] * node->GetCoord(i, 1, k)[2] -
					1.0 * coef->GetX(Idx(i, 2, k))[2] * node->GetCoord(i, 2, k)[2]);
				temp_eta[Idx(i, 0, k)].y_zeta_x = 0.5 * (3.0 * coef->GetY(Idx(i, 1, k))[2] * node->GetCoord(i, 1, k)[0] -
					1.0 * coef->GetY(Idx(i, 2, k))[2] * node->GetCoord(i, 2, k)[0]);
				temp_eta[Idx(i, 0, k)].z_xi_y = 0.5 * (3.0 * coef->GetZ(Idx(i, 1, k))[0] * node->GetCoord(i, 1, k)[1] -
					1.0 * coef->GetZ(Idx(i, 2, k))[0] * node->GetCoord(i, 2, k)[1]);
				temp_eta[Idx(i, 0, k)].x_xi_z = 0.5 * (3.0 * coef->GetX(Idx(i, 1, k))[0] * node->GetCoord(i, 1, k)[2] -
					1.0 * coef->GetX(Idx(i, 2, k))[0] * node->GetCoord(i, 2, k)[2]);
				temp_eta[Idx(i, 0, k)].y_xi_x = 0.5 * (3.0 * coef->GetY(Idx(i, 1, k))[0] * node->GetCoord(i, 1, k)[0] -
					1.0 * coef->GetY(Idx(i, 2, k))[0] * node->GetCoord(i, 2, k)[0]);

				idx = Idx(i, nj - 2, k);
				temp_eta[idx].y_zeta_z = 0.5 * (3.0 * coef->GetY(Idx(i, nj - 2, k))[2] * node->GetCoord(i, nj - 2, k)[2] -
					1.0 * coef->GetY(Idx(i, nj - 3, k))[2] * node->GetCoord(i, nj - 3, k)[2]);
				temp_eta[idx].z_zeta_x = 0.5 * (3.0 * coef->GetZ(Idx(i, nj - 2, k))[2] * node->GetCoord(i, nj - 2, k)[0] -
					1.0 * coef->GetZ(Idx(i, nj - 3, k))[2] * node->GetCoord(i, nj - 3, k)[0]);
				temp_eta[idx].x_zeta_y = 0.5 * (3.0 * coef->GetX(Idx(i, nj - 2, k))[2] * node->GetCoord(i, nj - 2, k)[1] -
					1.0 * coef->GetX(Idx(i, nj - 3, k))[2] * node->GetCoord(i, nj - 3, k)[1]);
				temp_eta[idx].y_xi_z = 0.5 * (3.0 * coef->GetY(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[2] -
					1.0 * coef->GetY(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[2]);
				temp_eta[idx].z_xi_x = 0.5 * (3.0 * coef->GetZ(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[0] -
					1.0 * coef->GetZ(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[0]);
				temp_eta[idx].x_xi_y = 0.5 * (3.0 * coef->GetX(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[1] -
					1.0 * coef->GetX(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[1]);
				temp_eta[Idx(i, nj - 2, k)].z_zeta_y =
					0.5 * (3.0 * coef->GetZ(Idx(i, nj - 2, k))[2] * node->GetCoord(i, nj - 2, k)[1] -
						1.0 * coef->GetZ(Idx(i, nj - 3, k))[2] * node->GetCoord(i, nj - 3, k)[1]);
				temp_eta[Idx(i, nj - 2, k)].x_zeta_z =
					0.5 * (3.0 * coef->GetX(Idx(i, nj - 2, k))[2] * node->GetCoord(i, nj - 2, k)[2] -
						1.0 * coef->GetX(Idx(i, nj - 3, k))[2] * node->GetCoord(i, nj - 3, k)[2]);
				temp_eta[Idx(i, nj - 2, k)].y_zeta_x =
					0.5 * (3.0 * coef->GetY(Idx(i, nj - 2, k))[2] * node->GetCoord(i, nj - 2, k)[0] -
						1.0 * coef->GetY(Idx(i, nj - 3, k))[2] * node->GetCoord(i, nj - 3, k)[0]);
				temp_eta[Idx(i, nj - 2, k)].z_xi_y =
					0.5 * (3.0 * coef->GetZ(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[1] -
						1.0 * coef->GetZ(Idx(i, nj - 3, k))[0] * node->GetCoord(i, nj - 3, k)[1]);
				temp_eta[Idx(i, nj - 2, k)].x_xi_z =
					0.5 * (3.0 * coef->GetX(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[2] -
						1.0 * coef->GetX(Idx(i, nj - 3, k))[0] * node->GetCoord(i, nj - 3, k)[2]);
				temp_eta[Idx(i, nj - 2, k)].y_xi_x =
					0.5 * (3.0 * coef->GetY(Idx(i, nj - 2, k))[0] * node->GetCoord(i, nj - 2, k)[0] -
						1.0 * coef->GetY(Idx(i, nj - 3, k))[0] * node->GetCoord(i, nj - 3, k)[0]);
			}
		}

		// 求出zeta方向的临时变量：逆变换度量系数乘以坐标
		struct TempZeta
		{
			double y_eta_z, z_eta_x, x_eta_y, y_xi_z, z_xi_x, x_xi_y, z_eta_y, x_eta_z, y_eta_x, z_xi_y, x_xi_z, y_xi_x;
		};
		dynamic_array<TempZeta> temp_zeta(ni * nj * nk);
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					int idx = Idx(i, j, k);
					auto& temp = temp_zeta[idx];
					temp.y_eta_z = 0.5 * (coef->GetY(Idx(i, j, k + 1))[1] * node->GetCoord(i, j, k + 1)[2] +
						coef->GetY(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[2]);
					temp.z_eta_x = 0.5 * (coef->GetZ(Idx(i, j, k + 1))[1] * node->GetCoord(i, j, k + 1)[0] +
						coef->GetZ(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[0]);
					temp.x_eta_y = 0.5 * (coef->GetX(Idx(i, j, k + 1))[1] * node->GetCoord(i, j, k + 1)[1] +
						coef->GetX(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[1]);
					temp.y_xi_z = 0.5 * (coef->GetY(Idx(i, j, k + 1))[0] * node->GetCoord(i, j, k + 1)[2] +
						coef->GetY(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[2]);
					temp.z_xi_x = 0.5 * (coef->GetZ(Idx(i, j, k + 1))[0] * node->GetCoord(i, j, k + 1)[0] +
						coef->GetZ(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[0]);
					temp.x_xi_y = 0.5 * (coef->GetX(Idx(i, j, k + 1))[0] * node->GetCoord(i, j, k + 1)[1] +
						coef->GetX(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[1]);
					temp.z_eta_y = 0.5 * (coef->GetZ(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[1] +
						coef->GetZ(Idx(i, j, k + 1))[1] * node->GetCoord(i, j, k + 1)[1]);
					temp.x_eta_z = 0.5 * (coef->GetX(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[2] +
						coef->GetX(Idx(i, j, k + 1))[1] * node->GetCoord(i, j, k + 1)[2]);
					temp.y_eta_x = 0.5 * (coef->GetY(Idx(i, j, k))[1] * node->GetCoord(i, j, k)[0] +
						coef->GetY(Idx(i, j, k + 1))[1] * node->GetCoord(i, j, k + 1)[0]);
					temp.z_xi_y = 0.5 * (coef->GetZ(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[1] +
						coef->GetZ(Idx(i, j, k + 1))[0] * node->GetCoord(i, j, k + 1)[1]);
					temp.x_xi_z = 0.5 * (coef->GetX(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[2] +
						coef->GetX(Idx(i, j, k + 1))[0] * node->GetCoord(i, j, k + 1)[2]);
					temp.y_xi_x = 0.5 * (coef->GetY(Idx(i, j, k))[0] * node->GetCoord(i, j, k)[0] +
						coef->GetY(Idx(i, j, k + 1))[0] * node->GetCoord(i, j, k + 1)[0]);
				}
			}
		}
		// k=1/2；(nk-1)-1/2
		for (int j = 1; j < nj - 1; ++j)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				int idx = Idx(i, j, 0);
				temp_zeta[idx].y_eta_z = 0.5 * (3.0 * coef->GetY(Idx(i, j, 1))[1] * node->GetCoord(i, j, 1)[2] -
					1.0 * coef->GetY(Idx(i, j, 2))[1] * node->GetCoord(i, j, 2)[2]);
				temp_zeta[idx].z_eta_x = 0.5 * (3.0 * coef->GetZ(Idx(i, j, 1))[1] * node->GetCoord(i, j, 1)[0] -
					1.0 * coef->GetZ(Idx(i, j, 2))[1] * node->GetCoord(i, j, 2)[0]);
				temp_zeta[idx].x_eta_y = 0.5 * (3.0 * coef->GetX(Idx(i, j, 1))[1] * node->GetCoord(i, j, 1)[1] -
					1.0 * coef->GetX(Idx(i, j, 2))[1] * node->GetCoord(i, j, 2)[1]);
				temp_zeta[idx].y_xi_z = 0.5 * (3.0 * coef->GetY(Idx(i, j, 1))[0] * node->GetCoord(i, j, 1)[2] -
					1.0 * coef->GetY(Idx(i, j, 2))[0] * node->GetCoord(i, j, 2)[2]);
				temp_zeta[idx].z_xi_x = 0.5 * (3.0 * coef->GetZ(Idx(i, j, 1))[0] * node->GetCoord(i, j, 1)[0] -
					1.0 * coef->GetZ(Idx(i, j, 2))[0] * node->GetCoord(i, j, 2)[0]);
				temp_zeta[idx].x_xi_y = 0.5 * (3.0 * coef->GetX(Idx(i, j, 1))[0] * node->GetCoord(i, j, 1)[1] -
					1.0 * coef->GetX(Idx(i, j, 2))[0] * node->GetCoord(i, j, 2)[1]);
				temp_zeta[Idx(i, j, 0)].z_eta_y = 0.5 * (3.0 * coef->GetZ(Idx(i, j, 1))[1] * node->GetCoord(i, j, 1)[1] -
					1.0 * coef->GetZ(Idx(i, j, 2))[1] * node->GetCoord(i, j, 2)[1]);
				temp_zeta[Idx(i, j, 0)].x_eta_z = 0.5 * (3.0 * coef->GetX(Idx(i, j, 1))[1] * node->GetCoord(i, j, 1)[2] -
					1.0 * coef->GetX(Idx(i, j, 2))[1] * node->GetCoord(i, j, 2)[2]);
				temp_zeta[Idx(i, j, 0)].y_eta_x = 0.5 * (3.0 * coef->GetY(Idx(i, j, 1))[1] * node->GetCoord(i, j, 1)[0] -
					1.0 * coef->GetY(Idx(i, j, 2))[1] * node->GetCoord(i, j, 2)[0]);
				temp_zeta[Idx(i, j, 0)].z_xi_y = 0.5 * (3.0 * coef->GetZ(Idx(i, j, 1))[0] * node->GetCoord(i, j, 1)[1] -
					1.0 * coef->GetZ(Idx(i, j, 2))[0] * node->GetCoord(i, j, 2)[1]);
				temp_zeta[Idx(i, j, 0)].x_xi_z = 0.5 * (3.0 * coef->GetX(Idx(i, j, 1))[0] * node->GetCoord(i, j, 1)[2] -
					1.0 * coef->GetX(Idx(i, j, 2))[0] * node->GetCoord(i, j, 2)[2]);
				temp_zeta[Idx(i, j, 0)].y_xi_x = 0.5 * (3.0 * coef->GetY(Idx(i, j, 1))[0] * node->GetCoord(i, j, 1)[0] -
					1.0 * coef->GetY(Idx(i, j, 2))[0] * node->GetCoord(i, j, 2)[0]);
				idx = Idx(i, j, nk - 2);
				temp_zeta[idx].y_eta_z = 0.5 * (3.0 * coef->GetY(Idx(i, j, nk - 2))[1] * node->GetCoord(i, j, nk - 2)[2] -
					1.0 * coef->GetY(Idx(i, j, nk - 3))[1] * node->GetCoord(i, j, nk - 3)[2]);
				temp_zeta[idx].z_eta_x = 0.5 * (3.0 * coef->GetZ(Idx(i, j, nk - 2))[1] * node->GetCoord(i, j, nk - 2)[0] -
					1.0 * coef->GetZ(Idx(i, j, nk - 3))[1] * node->GetCoord(i, j, nk - 3)[0]);
				temp_zeta[idx].x_eta_y = 0.5 * (3.0 * coef->GetX(Idx(i, j, nk - 2))[1] * node->GetCoord(i, j, nk - 2)[1] -
					1.0 * coef->GetX(Idx(i, j, nk - 3))[1] * node->GetCoord(i, j, nk - 3)[1]);
				temp_zeta[idx].y_xi_z = 0.5 * (3.0 * coef->GetY(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[2] -
					1.0 * coef->GetY(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[2]);
				temp_zeta[idx].z_xi_x = 0.5 * (3.0 * coef->GetZ(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[0] -
					1.0 * coef->GetZ(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[0]);
				temp_zeta[idx].x_xi_y = 0.5 * (3.0 * coef->GetX(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[1] -
					1.0 * coef->GetX(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[1]);
				temp_zeta[Idx(i, j, nk - 2)].z_eta_y =
					0.5 * (3.0 * coef->GetZ(Idx(i, j, nk - 2))[1] * node->GetCoord(i, j, nk - 2)[1] -
						1.0 * coef->GetZ(Idx(i, j, nk - 3))[1] * node->GetCoord(i, j, nk - 3)[1]);
				temp_zeta[Idx(i, j, nk - 2)].x_eta_z =
					0.5 * (3.0 * coef->GetX(Idx(i, j, nk - 2))[1] * node->GetCoord(i, j, nk - 2)[2] -
						1.0 * coef->GetX(Idx(i, j, nk - 3))[1] * node->GetCoord(i, j, nk - 3)[2]);
				temp_zeta[Idx(i, j, nk - 2)].y_eta_x =
					0.5 * (3.0 * coef->GetY(Idx(i, j, nk - 2))[1] * node->GetCoord(i, j, nk - 2)[0] -
						1.0 * coef->GetY(Idx(i, j, nk - 3))[1] * node->GetCoord(i, j, nk - 3)[0]);
				temp_zeta[Idx(i, j, nk - 2)].z_xi_y =
					0.5 * (3.0 * coef->GetZ(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[1] -
						1.0 * coef->GetZ(Idx(i, j, nk - 3))[0] * node->GetCoord(i, j, nk - 3)[1]);
				temp_zeta[Idx(i, j, nk - 2)].x_xi_z =
					0.5 * (3.0 * coef->GetX(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[2] -
						1.0 * coef->GetX(Idx(i, j, nk - 3))[0] * node->GetCoord(i, j, nk - 3)[2]);
				temp_zeta[Idx(i, j, nk - 2)].y_xi_x =
					0.5 * (3.0 * coef->GetY(Idx(i, j, nk - 2))[0] * node->GetCoord(i, j, nk - 2)[0] -
						1.0 * coef->GetY(Idx(i, j, nk - 3))[0] * node->GetCoord(i, j, nk - 3)[0]);
			}
		}

		// 第四步：根据半点坐标和半点逆变换度量系数使用守恒形式计算整点度量系数（CMM1)
		double temp[6][2];
		// 获取下标的lamda函数
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					// i direction
					auto coef_xi = coef->GetXi(Idx(i, j, k));
					auto coef_eta = coef->GetEta(Idx(i, j, k));
					auto coef_zeta = coef->GetZeta(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 2; iTemp++)
					{
						int idx = Idx(i - 1 + iTemp, j, k);
						temp[0][iTemp] = 0.5 * (temp_xi[idx].y_zeta_z - temp_xi[idx].z_zeta_y); // y_zeta*z-z_zeta*y
						temp[1][iTemp] = 0.5 * (temp_xi[idx].z_zeta_x - temp_xi[idx].x_zeta_z); // z_zeta*x-x_zeta*z
						temp[2][iTemp] = 0.5 * (temp_xi[idx].x_zeta_y - temp_xi[idx].y_zeta_x); // x_zeta*y-y_zeta*x
						temp[3][iTemp] = 0.5 * (temp_xi[idx].z_eta_y - temp_xi[idx].y_eta_z);   // z_eta*y-y_eta*z
						temp[4][iTemp] = 0.5 * (temp_xi[idx].x_eta_z - temp_xi[idx].z_eta_x);   // x_eta*z-z_eta*x
						temp[5][iTemp] = 0.5 * (temp_xi[idx].y_eta_x - temp_xi[idx].x_eta_y);   // y_eta*x-x_eta*y
					}
					coef_eta[0] += temp[0][1] - temp[0][0];
					coef_eta[1] += temp[1][1] - temp[1][0];
					coef_eta[2] += temp[2][1] - temp[2][0];
					coef_zeta[0] += temp[3][1] - temp[3][0];
					coef_zeta[1] += temp[4][1] - temp[4][0];
					coef_zeta[2] += temp[5][1] - temp[5][0];
					for (int iTemp = 0; iTemp < 2; iTemp++)
					{
						int idx = Idx(i, j - 1 + iTemp, k);
						temp[0][iTemp] = 0.5 * (temp_eta[idx].z_zeta_y - temp_eta[idx].y_zeta_z); // z_zeta*y-y_zeta*z
						temp[1][iTemp] = 0.5 * (temp_eta[idx].x_zeta_z - temp_eta[idx].z_zeta_x); // x_zeta*z-z_zeta*x
						temp[2][iTemp] = 0.5 * (temp_eta[idx].y_zeta_x - temp_eta[idx].x_zeta_y); //  y_zeta*x-x_zeta*y
						temp[3][iTemp] = 0.5 * (temp_eta[idx].y_xi_z - temp_eta[idx].z_xi_y);     // y_xi*z-z_xi*y
						temp[4][iTemp] = 0.5 * (temp_eta[idx].z_xi_x - temp_eta[idx].x_xi_z);     // z_xi*x-x_xi*z
						temp[5][iTemp] = 0.5 * (temp_eta[idx].x_xi_y - temp_eta[idx].y_xi_x);     // x_xi*y-y_xi*x
					}
					coef_xi[0] += temp[0][1] - temp[0][0];
					coef_xi[1] += temp[1][1] - temp[1][0];
					coef_xi[2] += temp[2][1] - temp[2][0];
					coef_zeta[0] += temp[3][1] - temp[3][0];
					coef_zeta[1] += temp[4][1] - temp[4][0];
					coef_zeta[2] += temp[5][1] - temp[5][0];
					for (int iTemp = 0; iTemp < 2; iTemp++)
					{
						int idx = Idx(i, j, k - 1 + iTemp);
						temp[0][iTemp] = 0.5 * (temp_zeta[idx].y_eta_z - temp_zeta[idx].z_eta_y); // y_eta*z-z_eta*y
						temp[1][iTemp] = 0.5 * (temp_zeta[idx].z_eta_x - temp_zeta[idx].x_eta_z); // z_eta*x-x_eta*z
						temp[2][iTemp] = 0.5 * (temp_zeta[idx].x_eta_y - temp_zeta[idx].y_eta_x); // x_eta*y-y_eta*x
						temp[3][iTemp] = 0.5 * (temp_zeta[idx].z_xi_y - temp_zeta[idx].y_xi_z);   // z_xi*y-y_xi*z
						temp[4][iTemp] = 0.5 * (temp_zeta[idx].x_xi_z - temp_zeta[idx].z_xi_x);   // x_xi*z-z_xi*x
						temp[5][iTemp] = 0.5 * (temp_zeta[idx].y_xi_x - temp_zeta[idx].x_xi_y);   //  y_xi*x-x_xi*y
					}
					coef_xi[0] += temp[0][1] - temp[0][0];
					coef_xi[1] += temp[1][1] - temp[1][0];
					coef_xi[2] += temp[2][1] - temp[2][0];
					coef_eta[0] += temp[3][1] - temp[3][0];
					coef_eta[1] += temp[4][1] - temp[4][0];
					coef_eta[2] += temp[5][1] - temp[5][0];
					if (grid->GetDim() == TWO_DIM)
					{
						coef_zeta[0] = 0.0;
						coef_zeta[1] = 0.0;
						coef_zeta[2] = 1.0;
					}
					coef_xi[3] = coef_eta[3] = coef_zeta[3] = 0.0;
				}
			}
		}
	}
	void NSSolverStruct::CalcMetricsS3_6th()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto coef = GetNodeMetrics();
		Metric* coef_mid[3] = { GetMidMetricsI(), GetMidMetricsJ(), GetMidMetricsK() };
		auto idx_proxy = GetIdxProxy();
		auto Idx = [&](int i, int j, int k) { return idx_proxy->GetIdx(i, j, k); };
		int ni = grid->GetNi();
		int nj = grid->GetNj();
		int nk = grid->GetNk();
		double inter_temp[6];
		// i+1/2,j+1/2,k+1/2处的坐标
		dynamic_array<dynamic_array<double>> coord_i(ni * nj * nk, dynamic_array<double>(3)),
			coord_j(ni * nj * nk, dynamic_array<double>(3)), coord_k(ni * nj * nk, dynamic_array<double>(3));
		// 计算之前，先把度量系数赋值为0
		for (int i = 0; i < ni; ++i)
		{
			for (int j = 0; j < nj; ++j)
			{
				for (int k = 0; k < nk; ++k)
				{
					for (int iDim = 0; iDim < 4; ++iDim)
					{
						coef->GetX(Idx(i, j, k))[iDim] = 0.0;
						coef->GetY(Idx(i, j, k))[iDim] = 0.0;
						coef->GetZ(Idx(i, j, k))[iDim] = 0.0;
						coef->GetXi(Idx(i, j, k))[iDim] = 0.0;
						coef->GetEta(Idx(i, j, k))[iDim] = 0.0;
						coef->GetZeta(Idx(i, j, k))[iDim] = 0.0;
						for (int jDim = 0; jDim < 3; ++jDim)
						{
							coef_mid[jDim]->GetX(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetY(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetZ(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetXi(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetEta(Idx(i, j, k))[iDim] = 0.0;
							coef_mid[jDim]->GetZeta(Idx(i, j, k))[iDim] = 0.0;
						}
					}
				}
			}
		}
		// 第一步：计算半点坐标
		// i+1/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 3; i < ni - 4; ++i)
				{
					for (int iDim = 0; iDim < 3; ++iDim)
					{
						for (int iTemp = 0; iTemp < 6; iTemp++)
						{
							inter_temp[iTemp] = node->GetCoord(i + iTemp - 2, j, k)[iDim];
						}
						coord_i[Idx(i, j, k)][iDim] = MidNodeInter6th(inter_temp);
					}
				}
			}
		}
		// i=1/2,3/2,5/2;ni-3/2,ni-5/2,ni-7/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int iDim = 0; iDim < 3; ++iDim)
				{
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = node->GetCoord(1 + iTemp, j, k)[iDim];
					}
					coord_i[Idx(0, j, k)][iDim] = MidNodeInter4thRight1(inter_temp);
					coord_i[Idx(1, j, k)][iDim] = MidNodeInter4thRight2(inter_temp);
					coord_i[Idx(2, j, k)][iDim] = MidNodeInter4th(inter_temp);
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = node->GetCoord(ni - 5 + iTemp, j, k)[iDim];
					}
					coord_i[Idx(ni - 2, j, k)][iDim] = MidNodeInter4thLeft1(inter_temp);
					coord_i[Idx(ni - 3, j, k)][iDim] = MidNodeInter4thLeft2(inter_temp);
					coord_i[Idx(ni - 4, j, k)][iDim] = MidNodeInter4th(inter_temp);
				}
			}
		}
		// j+1/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 3; j < nj - 4; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					for (int iDim = 0; iDim < 3; ++iDim)
					{
						for (int iTemp = 0; iTemp < 6; iTemp++)
						{
							inter_temp[iTemp] = node->GetCoord(i, j + iTemp - 2, k)[iDim];
						}
						coord_j[Idx(i, j, k)][iDim] = MidNodeInter6th(inter_temp);
					}
				}
			}
		}
		// j=1/2,3/2,5/2;nj-3/2,nj-5/2,nj-7/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				for (int iDim = 0; iDim < 3; ++iDim)
				{
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = node->GetCoord(i, 1 + iTemp, k)[iDim];
					}
					coord_j[Idx(i, 0, k)][iDim] = MidNodeInter4thRight1(inter_temp);
					coord_j[Idx(i, 1, k)][iDim] = MidNodeInter4thRight2(inter_temp);
					coord_j[Idx(i, 2, k)][iDim] = MidNodeInter4th(inter_temp);
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = node->GetCoord(i, nj - 5 + iTemp, k)[iDim];
					}
					coord_j[Idx(i, nj - 2, k)][iDim] = MidNodeInter4thLeft1(inter_temp);
					coord_j[Idx(i, nj - 3, k)][iDim] = MidNodeInter4thLeft2(inter_temp);
					coord_j[Idx(i, nj - 4, k)][iDim] = MidNodeInter4th(inter_temp);
				}
			}
		}
		// k+1/2
		for (int k = 3; k < nk - 4; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					for (int iDim = 0; iDim < 3; ++iDim)
					{
						for (int iTemp = 0; iTemp < 6; iTemp++)
						{
							inter_temp[iTemp] = node->GetCoord(i, j, k + iTemp - 2)[iDim];
						}
						coord_k[Idx(i, j, k)][iDim] = MidNodeInter6th(inter_temp);
					}
				}
			}
		}
		// k=1/2,3/2,5/2;nk-3/2,nk-5/2,nk-7/2
		for (int j = 1; j < nj - 1; ++j)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				for (int iDim = 0; iDim < 3; ++iDim)
				{
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = node->GetCoord(i, j, 1 + iTemp)[iDim];
					}
					coord_k[Idx(i, j, 0)][iDim] = MidNodeInter4thRight1(inter_temp);
					coord_k[Idx(i, j, 1)][iDim] = MidNodeInter4thRight2(inter_temp);
					coord_k[Idx(i, j, 2)][iDim] = MidNodeInter4th(inter_temp);
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = node->GetCoord(i, j, nk - 5 + iTemp)[iDim];
					}
					coord_k[Idx(i, j, nk - 2)][iDim] = MidNodeInter4thLeft1(inter_temp);
					coord_k[Idx(i, j, nk - 3)][iDim] = MidNodeInter4thLeft2(inter_temp);
					coord_k[Idx(i, j, nk - 4)][iDim] = MidNodeInter4th(inter_temp);
				}
			}
		}
		// 第二步：根据半点坐标计算整点逆变换度量系数
		double diff_temp[6];
		// i direction
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 3; i < ni - 3; ++i)
				{
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_i[Idx(i + iTemp - 3, j, k)][0];
					}
					auto coef_x = coef->GetX(Idx(i, j, k));
					coef_x[0] = NodeDifferece6th(diff_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_i[Idx(i + iTemp - 3, j, k)][1];
					}
					auto coef_y = coef->GetY(Idx(i, j, k));
					coef_y[0] = NodeDifferece6th(diff_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_i[Idx(i + iTemp - 3, j, k)][2];
					}
					auto coef_z = coef->GetZ(Idx(i, j, k));
					coef_z[0] = NodeDifferece6th(diff_temp);
				}
			}
		}
		// j direction
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 3; j < nj - 3; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					auto coef_x = coef->GetX(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_j[Idx(i, j + iTemp - 3, k)][0];
					}
					coef_x[1] = NodeDifferece6th(diff_temp);
					auto coef_y = coef->GetY(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_j[Idx(i, j + iTemp - 3, k)][1];
					}
					coef_y[1] = NodeDifferece6th(diff_temp);
					auto coef_z = coef->GetZ(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_j[Idx(i, j + iTemp - 3, k)][2];
					}
					coef_z[1] = NodeDifferece6th(diff_temp);
				}
			}
		}
		// k direction
		for (int k = 3; k < nk - 3; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					auto coef_x = coef->GetX(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_k[Idx(i, j, k + iTemp - 3)][0];
					}
					coef_x[2] = NodeDifferece6th(diff_temp);
					auto coef_y = coef->GetY(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_k[Idx(i, j, k + iTemp - 3)][1];
					}
					coef_y[2] = NodeDifferece6th(diff_temp);
					auto coef_z = coef->GetZ(Idx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						diff_temp[iTemp] = coord_k[Idx(i, j, k + iTemp - 3)][2];
					}
					coef_z[2] = NodeDifferece6th(diff_temp);
					if (grid->GetDim() == TWO_DIM)
					{
						coef_x[2] = 0.0;
						coef_y[2] = 0.0;
						coef_z[2] = 1.0;
					}
				}
			}
		}
		// 计算边界整点的度量系数
		// i=1,2;ni-3,ni-2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_i[Idx(0 + iTemp, j, k)][0];
				}
				coef->GetX(Idx(1, j, k))[0] = NodeDifferece4thRight(diff_temp);
				coef->GetX(Idx(2, j, k))[0] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_i[Idx(0 + iTemp, j, k)][1];
				}
				coef->GetY(Idx(1, j, k))[0] = NodeDifferece4thRight(diff_temp);
				coef->GetY(Idx(2, j, k))[0] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_i[Idx(0 + iTemp, j, k)][2];
				}
				coef->GetZ(Idx(1, j, k))[0] = NodeDifferece4thRight(diff_temp);
				coef->GetZ(Idx(2, j, k))[0] = NodeDifferece4th(diff_temp);

				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_i[Idx(ni - 2 - iTemp, j, k)][0];
				}
				coef->GetX(Idx(ni - 2, j, k))[0] = NodeDifferece4thLeft(diff_temp);
				coef->GetX(Idx(ni - 3, j, k))[0] = -NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_i[Idx(ni - 2 - iTemp, j, k)][1];
				}
				coef->GetY(Idx(ni - 2, j, k))[0] = NodeDifferece4thLeft(diff_temp);
				coef->GetY(Idx(ni - 3, j, k))[0] = -NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_i[Idx(ni - 2 - iTemp, j, k)][2];
				}
				coef->GetZ(Idx(ni - 2, j, k))[0] = NodeDifferece4thLeft(diff_temp);
				coef->GetZ(Idx(ni - 3, j, k))[0] = -NodeDifferece4th(diff_temp);
			}
		}
		// j=1,2;nj-3,nj-2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_j[Idx(i, 0 + iTemp, k)][0];
				}
				coef->GetX(Idx(i, 1, k))[1] = NodeDifferece4thRight(diff_temp);
				coef->GetX(Idx(i, 2, k))[1] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_j[Idx(i, 0 + iTemp, k)][1];
				}
				coef->GetY(Idx(i, 1, k))[1] = NodeDifferece4thRight(diff_temp);
				coef->GetY(Idx(i, 2, k))[1] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_j[Idx(i, 0 + iTemp, k)][2];
				}
				coef->GetZ(Idx(i, 1, k))[1] = NodeDifferece4thRight(diff_temp);
				coef->GetZ(Idx(i, 2, k))[1] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_j[Idx(i, nj - 2 - iTemp, k)][0];
				}
				coef->GetX(Idx(i, nj - 2, k))[1] = NodeDifferece4thLeft(diff_temp);
				coef->GetX(Idx(i, nj - 3, k))[1] = -NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_j[Idx(i, nj - 2 - iTemp, k)][1];
				}
				coef->GetY(Idx(i, nj - 2, k))[1] = NodeDifferece4thLeft(diff_temp);
				coef->GetY(Idx(i, nj - 3, k))[1] = -NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_j[Idx(i, nj - 2 - iTemp, k)][2];
				}
				coef->GetZ(Idx(i, nj - 2, k))[1] = NodeDifferece4thLeft(diff_temp);
				coef->GetZ(Idx(i, nj - 3, k))[1] = -NodeDifferece4th(diff_temp);
			}
		}
		// k=1,2;nk-3,nk-2
		for (int j = 1; j < nj - 1; ++j)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_k[Idx(i, j, 0 + iTemp)][0];
				}
				coef->GetX(Idx(i, j, 1))[2] = NodeDifferece4thRight(diff_temp);
				coef->GetX(Idx(i, j, 2))[2] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_k[Idx(i, j, 0 + iTemp)][1];
				}
				coef->GetY(Idx(i, j, 1))[2] = NodeDifferece4thRight(diff_temp);
				coef->GetY(Idx(i, j, 2))[2] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_k[Idx(i, j, 0 + iTemp)][2];
				}
				coef->GetZ(Idx(i, j, 1))[2] = NodeDifferece4thRight(diff_temp);
				coef->GetZ(Idx(i, j, 2))[2] = NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_k[Idx(i, j, nk - 2 - iTemp)][0];
				}
				coef->GetX(Idx(i, j, nk - 2))[2] = NodeDifferece4thLeft(diff_temp);
				coef->GetX(Idx(i, j, nk - 3))[2] = -NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_k[Idx(i, j, nk - 2 - iTemp)][1];
				}
				coef->GetY(Idx(i, j, nk - 2))[2] = NodeDifferece4thLeft(diff_temp);
				coef->GetY(Idx(i, j, nk - 3))[2] = -NodeDifferece4th(diff_temp);
				for (int iTemp = 0; iTemp < 6; iTemp++)
				{
					diff_temp[iTemp] = coord_k[Idx(i, j, nk - 2 - iTemp)][2];
				}
				coef->GetZ(Idx(i, j, nk - 2))[2] = NodeDifferece4thLeft(diff_temp);
				coef->GetZ(Idx(i, j, nk - 3))[2] = -NodeDifferece4th(diff_temp);
				if (grid->GetDim() == TWO_DIM)
				{
					coef->GetZ(Idx(i, j, 1))[0] = coef->GetZ(Idx(i, j, 1))[1] = 0.0;
					coef->GetZ(Idx(i, j, 2))[0] = coef->GetZ(Idx(i, j, 2))[1] = 0.0;
					coef->GetZ(Idx(i, j, nk - 2))[0] = coef->GetZ(Idx(i, j, nk - 2))[1] = 0.0;
					coef->GetZ(Idx(i, j, nk - 3))[0] = coef->GetZ(Idx(i, j, nk - 3))[1] = 0.0;
					coef->GetZ(Idx(i, j, 1))[2] = 1.0;
					coef->GetZ(Idx(i, j, 2))[2] = 1.0;
					coef->GetZ(Idx(i, j, nk - 2))[2] = 1.0;
					coef->GetZ(Idx(i, j, nk - 3))[2] = 1.0;
				}
			}
		}

		// 第三步：计算半点处逆变换度量系数和坐标乘积
		//  求出xi方向的临时变量：逆变换度量系数乘以坐标
		struct TempXi
		{
			double y_zeta_z, z_zeta_x, x_zeta_y, y_eta_z, z_eta_x, x_eta_y, z_zeta_y, x_zeta_z, y_zeta_x, z_eta_y, x_eta_z,
				y_eta_x;
		};
		dynamic_array<TempXi> temp_xi(ni * nj * nk);
		// i+1/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 3; i < ni - 4; ++i)
				{
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetY(Idx(i + iTemp - 2, j, k))[2] * node->GetCoord(i + iTemp - 2, j, k)[2];
					}
					temp_xi[Idx(i, j, k)].y_zeta_z = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetZ(Idx(i + iTemp - 2, j, k))[2] * node->GetCoord(i + iTemp - 2, j, k)[0];
					}
					temp_xi[Idx(i, j, k)].z_zeta_x = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetX(Idx(i + iTemp - 2, j, k))[2] * node->GetCoord(i + iTemp - 2, j, k)[1];
					}
					temp_xi[Idx(i, j, k)].x_zeta_y = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetY(Idx(i + iTemp - 2, j, k))[1] * node->GetCoord(i + iTemp - 2, j, k)[2];
					}
					temp_xi[Idx(i, j, k)].y_eta_z = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetZ(Idx(i + iTemp - 2, j, k))[1] * node->GetCoord(i + iTemp - 2, j, k)[0];
					}
					temp_xi[Idx(i, j, k)].z_eta_x = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetX(Idx(i + iTemp - 2, j, k))[1] * node->GetCoord(i + iTemp - 2, j, k)[1];
					}
					temp_xi[Idx(i, j, k)].x_eta_y = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetZ(Idx(i - 2 + iTemp, j, k))[2] * node->GetCoord(i - 2 + iTemp, j, k)[1];
					}
					temp_xi[Idx(i, j, k)].z_zeta_y = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetX(Idx(i - 2 + iTemp, j, k))[2] * node->GetCoord(i - 2 + iTemp, j, k)[2];
					}
					temp_xi[Idx(i, j, k)].x_zeta_z = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetY(Idx(i - 2 + iTemp, j, k))[2] * node->GetCoord(i - 2 + iTemp, j, k)[0];
					}
					temp_xi[Idx(i, j, k)].y_zeta_x = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetZ(Idx(i - 2 + iTemp, j, k))[1] * node->GetCoord(i - 2 + iTemp, j, k)[1];
					}
					temp_xi[Idx(i, j, k)].z_eta_y = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetX(Idx(i - 2 + iTemp, j, k))[1] * node->GetCoord(i - 2 + iTemp, j, k)[2];
					}
					temp_xi[Idx(i, j, k)].x_eta_z = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetY(Idx(i - 2 + iTemp, j, k))[1] * node->GetCoord(i - 2 + iTemp, j, k)[0];
					}
					temp_xi[Idx(i, j, k)].y_eta_x = MidNodeInter6th(inter_temp);
				}
			}
		}
		// i=1/2,3/2,5/2;ni-3/2,ni-5/2,ni-7/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				// i=1/2，3/2，5/2
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(1 + iTemp, j, k))[2] * node->GetCoord(1 + iTemp, j, k)[2];
				}
				temp_xi[Idx(0, j, k)].y_zeta_z = MidNodeInter4thRight1(inter_temp);
				temp_xi[Idx(1, j, k)].y_zeta_z = MidNodeInter4thRight2(inter_temp);
				temp_xi[Idx(2, j, k)].y_zeta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(1 + iTemp, j, k))[2] * node->GetCoord(1 + iTemp, j, k)[0];
				}
				temp_xi[Idx(0, j, k)].z_zeta_x = MidNodeInter4thRight1(inter_temp);
				temp_xi[Idx(1, j, k)].z_zeta_x = MidNodeInter4thRight2(inter_temp);
				temp_xi[Idx(2, j, k)].z_zeta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(1 + iTemp, j, k))[2] * node->GetCoord(1 + iTemp, j, k)[1];
				}
				temp_xi[Idx(0, j, k)].x_zeta_y = MidNodeInter4thRight1(inter_temp);
				temp_xi[Idx(1, j, k)].x_zeta_y = MidNodeInter4thRight2(inter_temp);
				temp_xi[Idx(2, j, k)].x_zeta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(1 + iTemp, j, k))[1] * node->GetCoord(1 + iTemp, j, k)[2];
				}
				temp_xi[Idx(0, j, k)].y_eta_z = MidNodeInter4thRight1(inter_temp);
				temp_xi[Idx(1, j, k)].y_eta_z = MidNodeInter4thRight2(inter_temp);
				temp_xi[Idx(2, j, k)].y_eta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(1 + iTemp, j, k))[1] * node->GetCoord(1 + iTemp, j, k)[0];
				}
				temp_xi[Idx(0, j, k)].z_eta_x = MidNodeInter4thRight1(inter_temp);
				temp_xi[Idx(1, j, k)].z_eta_x = MidNodeInter4thRight2(inter_temp);
				temp_xi[Idx(2, j, k)].z_eta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(1 + iTemp, j, k))[1] * node->GetCoord(1 + iTemp, j, k)[1];
				}
				temp_xi[Idx(0, j, k)].x_eta_y = MidNodeInter4thRight1(inter_temp);
				temp_xi[Idx(1, j, k)].x_eta_y = MidNodeInter4thRight2(inter_temp);
				temp_xi[Idx(2, j, k)].x_eta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(1 + iTemp, j, k))[2] * node->GetCoord(1 + iTemp, j, k)[1];
				}
				temp_xi[Idx(0, j, k)].z_zeta_y = MidNodeInter4thRight1(inter_temp);
				temp_xi[Idx(1, j, k)].z_zeta_y = MidNodeInter4thRight2(inter_temp);
				temp_xi[Idx(2, j, k)].z_zeta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(1 + iTemp, j, k))[2] * node->GetCoord(1 + iTemp, j, k)[2];
				}
				temp_xi[Idx(0, j, k)].x_zeta_z = MidNodeInter4thRight1(inter_temp);
				temp_xi[Idx(1, j, k)].x_zeta_z = MidNodeInter4thRight2(inter_temp);
				temp_xi[Idx(2, j, k)].x_zeta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(1 + iTemp, j, k))[2] * node->GetCoord(1 + iTemp, j, k)[0];
				}
				temp_xi[Idx(0, j, k)].y_zeta_x = MidNodeInter4thRight1(inter_temp);
				temp_xi[Idx(1, j, k)].y_zeta_x = MidNodeInter4thRight2(inter_temp);
				temp_xi[Idx(2, j, k)].y_zeta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(1 + iTemp, j, k))[1] * node->GetCoord(1 + iTemp, j, k)[1];
				}
				temp_xi[Idx(0, j, k)].z_eta_y = MidNodeInter4thRight1(inter_temp);
				temp_xi[Idx(1, j, k)].z_eta_y = MidNodeInter4thRight2(inter_temp);
				temp_xi[Idx(2, j, k)].z_eta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(1 + iTemp, j, k))[1] * node->GetCoord(1 + iTemp, j, k)[2];
				}
				temp_xi[Idx(0, j, k)].x_eta_z = MidNodeInter4thRight1(inter_temp);
				temp_xi[Idx(1, j, k)].x_eta_z = MidNodeInter4thRight2(inter_temp);
				temp_xi[Idx(2, j, k)].x_eta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(1 + iTemp, j, k))[1] * node->GetCoord(1 + iTemp, j, k)[0];
				}
				temp_xi[Idx(0, j, k)].y_eta_x = MidNodeInter4thRight1(inter_temp);
				temp_xi[Idx(1, j, k)].y_eta_x = MidNodeInter4thRight2(inter_temp);
				temp_xi[Idx(2, j, k)].y_eta_x = MidNodeInter4th(inter_temp);
				// i=ni-3/2,ni-5/2,ni-7/2
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(ni - 5 + iTemp, j, k))[2] * node->GetCoord(ni - 5 + iTemp, j, k)[2];
				}
				temp_xi[Idx(ni - 2, j, k)].y_zeta_z = MidNodeInter4thLeft1(inter_temp);
				temp_xi[Idx(ni - 3, j, k)].y_zeta_z = MidNodeInter4thLeft2(inter_temp);
				temp_xi[Idx(ni - 4, j, k)].y_zeta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(ni - 5 + iTemp, j, k))[2] * node->GetCoord(ni - 5 + iTemp, j, k)[0];
				}
				temp_xi[Idx(ni - 2, j, k)].z_zeta_x = MidNodeInter4thLeft1(inter_temp);
				temp_xi[Idx(ni - 3, j, k)].z_zeta_x = MidNodeInter4thLeft2(inter_temp);
				temp_xi[Idx(ni - 4, j, k)].z_zeta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(ni - 5 + iTemp, j, k))[2] * node->GetCoord(ni - 5 + iTemp, j, k)[1];
				}
				temp_xi[Idx(ni - 2, j, k)].x_zeta_y = MidNodeInter4thLeft1(inter_temp);
				temp_xi[Idx(ni - 3, j, k)].x_zeta_y = MidNodeInter4thLeft2(inter_temp);
				temp_xi[Idx(ni - 4, j, k)].x_zeta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(ni - 5 + iTemp, j, k))[1] * node->GetCoord(ni - 5 + iTemp, j, k)[2];
				}
				temp_xi[Idx(ni - 2, j, k)].y_eta_z = MidNodeInter4thLeft1(inter_temp);
				temp_xi[Idx(ni - 3, j, k)].y_eta_z = MidNodeInter4thLeft2(inter_temp);
				temp_xi[Idx(ni - 4, j, k)].y_eta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(ni - 5 + iTemp, j, k))[1] * node->GetCoord(ni - 5 + iTemp, j, k)[0];
				}
				temp_xi[Idx(ni - 2, j, k)].z_eta_x = MidNodeInter4thLeft1(inter_temp);
				temp_xi[Idx(ni - 3, j, k)].z_eta_x = MidNodeInter4thLeft2(inter_temp);
				temp_xi[Idx(ni - 4, j, k)].z_eta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(ni - 5 + iTemp, j, k))[1] * node->GetCoord(ni - 5 + iTemp, j, k)[1];
				}
				temp_xi[Idx(ni - 2, j, k)].x_eta_y = MidNodeInter4thLeft1(inter_temp);
				temp_xi[Idx(ni - 3, j, k)].x_eta_y = MidNodeInter4thLeft2(inter_temp);
				temp_xi[Idx(ni - 4, j, k)].x_eta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(ni - 5 + iTemp, j, k))[2] * node->GetCoord(ni - 5 + iTemp, j, k)[1];
				}
				temp_xi[Idx(ni - 2, j, k)].z_zeta_y = MidNodeInter4thLeft1(inter_temp);
				temp_xi[Idx(ni - 3, j, k)].z_zeta_y = MidNodeInter4thLeft2(inter_temp);
				temp_xi[Idx(ni - 4, j, k)].z_zeta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(ni - 5 + iTemp, j, k))[2] * node->GetCoord(ni - 5 + iTemp, j, k)[2];
				}
				temp_xi[Idx(ni - 2, j, k)].x_zeta_z = MidNodeInter4thLeft1(inter_temp);
				temp_xi[Idx(ni - 3, j, k)].x_zeta_z = MidNodeInter4thLeft2(inter_temp);
				temp_xi[Idx(ni - 4, j, k)].x_zeta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(ni - 5 + iTemp, j, k))[2] * node->GetCoord(ni - 5 + iTemp, j, k)[0];
				}
				temp_xi[Idx(ni - 2, j, k)].y_zeta_x = MidNodeInter4thLeft1(inter_temp);
				temp_xi[Idx(ni - 3, j, k)].y_zeta_x = MidNodeInter4thLeft2(inter_temp);
				temp_xi[Idx(ni - 4, j, k)].y_zeta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(ni - 5 + iTemp, j, k))[1] * node->GetCoord(ni - 5 + iTemp, j, k)[1];
				}
				temp_xi[Idx(ni - 2, j, k)].z_eta_y = MidNodeInter4thLeft1(inter_temp);
				temp_xi[Idx(ni - 3, j, k)].z_eta_y = MidNodeInter4thLeft2(inter_temp);
				temp_xi[Idx(ni - 4, j, k)].z_eta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(ni - 5 + iTemp, j, k))[1] * node->GetCoord(ni - 5 + iTemp, j, k)[2];
				}
				temp_xi[Idx(ni - 2, j, k)].x_eta_z = MidNodeInter4thLeft1(inter_temp);
				temp_xi[Idx(ni - 3, j, k)].x_eta_z = MidNodeInter4thLeft2(inter_temp);
				temp_xi[Idx(ni - 4, j, k)].x_eta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(ni - 5 + iTemp, j, k))[1] * node->GetCoord(ni - 5 + iTemp, j, k)[0];
				}
				temp_xi[Idx(ni - 2, j, k)].y_eta_x = MidNodeInter4thLeft1(inter_temp);
				temp_xi[Idx(ni - 3, j, k)].y_eta_x = MidNodeInter4thLeft2(inter_temp);
				temp_xi[Idx(ni - 4, j, k)].y_eta_x = MidNodeInter4th(inter_temp);
			}
		}
		// 求出eta方向的临时变量：逆变换度量系数乘以坐标
		struct TempEta
		{
			double y_zeta_z, z_zeta_x, x_zeta_y, y_xi_z, z_xi_x, x_xi_y, z_zeta_y, x_zeta_z, y_zeta_x, z_xi_y, x_xi_z,
				y_xi_x;
		};
		dynamic_array<TempEta> temp_eta(ni * nj * nk);
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 3; j < nj - 4; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetY(Idx(i, j + iTemp - 2, k))[2] * node->GetCoord(i, j + iTemp - 2, k)[2];
					}
					temp_eta[Idx(i, j, k)].y_zeta_z = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetZ(Idx(i, j + iTemp - 2, k))[2] * node->GetCoord(i, j + iTemp - 2, k)[0];
					}
					temp_eta[Idx(i, j, k)].z_zeta_x = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetX(Idx(i, j + iTemp - 2, k))[2] * node->GetCoord(i, j + iTemp - 2, k)[1];
					}
					temp_eta[Idx(i, j, k)].x_zeta_y = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetY(Idx(i, j + iTemp - 2, k))[0] * node->GetCoord(i, j + iTemp - 2, k)[2];
					}
					temp_eta[Idx(i, j, k)].y_xi_z = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetZ(Idx(i, j + iTemp - 2, k))[0] * node->GetCoord(i, j + iTemp - 2, k)[0];
					}
					temp_eta[Idx(i, j, k)].z_xi_x = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetX(Idx(i, j + iTemp - 2, k))[0] * node->GetCoord(i, j + iTemp - 2, k)[1];
					}
					temp_eta[Idx(i, j, k)].x_xi_y = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetZ(Idx(i, j - 2 + iTemp, k))[2] * node->GetCoord(i, j - 2 + iTemp, k)[1];
					}
					temp_eta[Idx(i, j, k)].z_zeta_y = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetX(Idx(i, j - 2 + iTemp, k))[2] * node->GetCoord(i, j - 2 + iTemp, k)[2];
					}
					temp_eta[Idx(i, j, k)].x_zeta_z = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetY(Idx(i, j - 2 + iTemp, k))[2] * node->GetCoord(i, j - 2 + iTemp, k)[0];
					}
					temp_eta[Idx(i, j, k)].y_zeta_x = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetZ(Idx(i, j - 2 + iTemp, k))[0] * node->GetCoord(i, j - 2 + iTemp, k)[1];
					}
					temp_eta[Idx(i, j, k)].z_xi_y = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetX(Idx(i, j - 2 + iTemp, k))[0] * node->GetCoord(i, j - 2 + iTemp, k)[2];
					}
					temp_eta[Idx(i, j, k)].x_xi_z = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetY(Idx(i, j - 2 + iTemp, k))[0] * node->GetCoord(i, j - 2 + iTemp, k)[0];
					}
					temp_eta[Idx(i, j, k)].y_xi_x = MidNodeInter6th(inter_temp);
				}
			}
		}
		// j=1/2,3/2,5/2;nj-3/2,nj-5/2,nj-7/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				// j=1/2，3/2，5/2
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, 1 + iTemp, k))[2] * node->GetCoord(i, 1 + iTemp, k)[2];
				}
				temp_eta[Idx(i, 0, k)].y_zeta_z = MidNodeInter4thRight1(inter_temp);
				temp_eta[Idx(i, 1, k)].y_zeta_z = MidNodeInter4thRight2(inter_temp);
				temp_eta[Idx(i, 2, k)].y_zeta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, 1 + iTemp, k))[2] * node->GetCoord(i, 1 + iTemp, k)[0];
				}
				temp_eta[Idx(i, 0, k)].z_zeta_x = MidNodeInter4thRight1(inter_temp);
				temp_eta[Idx(i, 1, k)].z_zeta_x = MidNodeInter4thRight2(inter_temp);
				temp_eta[Idx(i, 2, k)].z_zeta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, 1 + iTemp, k))[2] * node->GetCoord(i, 1 + iTemp, k)[1];
				}
				temp_eta[Idx(i, 0, k)].x_zeta_y = MidNodeInter4thRight1(inter_temp);
				temp_eta[Idx(i, 1, k)].x_zeta_y = MidNodeInter4thRight2(inter_temp);
				temp_eta[Idx(i, 2, k)].x_zeta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, 1 + iTemp, k))[0] * node->GetCoord(i, 1 + iTemp, k)[2];
				}
				temp_eta[Idx(i, 0, k)].y_xi_z = MidNodeInter4thRight1(inter_temp);
				temp_eta[Idx(i, 1, k)].y_xi_z = MidNodeInter4thRight2(inter_temp);
				temp_eta[Idx(i, 2, k)].y_xi_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, 1 + iTemp, k))[0] * node->GetCoord(i, 1 + iTemp, k)[0];
				}
				temp_eta[Idx(i, 0, k)].z_xi_x = MidNodeInter4thRight1(inter_temp);
				temp_eta[Idx(i, 1, k)].z_xi_x = MidNodeInter4thRight2(inter_temp);
				temp_eta[Idx(i, 2, k)].z_xi_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, 1 + iTemp, k))[0] * node->GetCoord(i, 1 + iTemp, k)[1];
				}
				temp_eta[Idx(i, 0, k)].x_xi_y = MidNodeInter4thRight1(inter_temp);
				temp_eta[Idx(i, 1, k)].x_xi_y = MidNodeInter4thRight2(inter_temp);
				temp_eta[Idx(i, 2, k)].x_xi_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, 1 + iTemp, k))[2] * node->GetCoord(i, 1 + iTemp, k)[1];
				}
				temp_eta[Idx(i, 0, k)].z_zeta_y = MidNodeInter4thRight1(inter_temp);
				temp_eta[Idx(i, 1, k)].z_zeta_y = MidNodeInter4thRight2(inter_temp);
				temp_eta[Idx(i, 2, k)].z_zeta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, 1 + iTemp, k))[2] * node->GetCoord(i, 1 + iTemp, k)[2];
				}
				temp_eta[Idx(i, 0, k)].x_zeta_z = MidNodeInter4thRight1(inter_temp);
				temp_eta[Idx(i, 1, k)].x_zeta_z = MidNodeInter4thRight2(inter_temp);
				temp_eta[Idx(i, 2, k)].x_zeta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, 1 + iTemp, k))[2] * node->GetCoord(i, 1 + iTemp, k)[0];
				}
				temp_eta[Idx(i, 0, k)].y_zeta_x = MidNodeInter4thRight1(inter_temp);
				temp_eta[Idx(i, 1, k)].y_zeta_x = MidNodeInter4thRight2(inter_temp);
				temp_eta[Idx(i, 2, k)].y_zeta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, 1 + iTemp, k))[0] * node->GetCoord(i, 1 + iTemp, k)[1];
				}
				temp_eta[Idx(i, 0, k)].z_xi_y = MidNodeInter4thRight1(inter_temp);
				temp_eta[Idx(i, 1, k)].z_xi_y = MidNodeInter4thRight2(inter_temp);
				temp_eta[Idx(i, 2, k)].z_xi_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, 1 + iTemp, k))[0] * node->GetCoord(i, 1 + iTemp, k)[2];
				}
				temp_eta[Idx(i, 0, k)].x_xi_z = MidNodeInter4thRight1(inter_temp);
				temp_eta[Idx(i, 1, k)].x_xi_z = MidNodeInter4thRight2(inter_temp);
				temp_eta[Idx(i, 2, k)].x_xi_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, 1 + iTemp, k))[0] * node->GetCoord(i, 1 + iTemp, k)[0];
				}
				temp_eta[Idx(i, 0, k)].y_xi_x = MidNodeInter4thRight1(inter_temp);
				temp_eta[Idx(i, 1, k)].y_xi_x = MidNodeInter4thRight2(inter_temp);
				temp_eta[Idx(i, 2, k)].y_xi_x = MidNodeInter4th(inter_temp);
				// j=nj-3/2,nj-5/2,nj-7/2
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, nj - 5 + iTemp, k))[2] * node->GetCoord(i, nj - 5 + iTemp, k)[2];
				}
				temp_eta[Idx(i, nj - 2, k)].y_zeta_z = MidNodeInter4thLeft1(inter_temp);
				temp_eta[Idx(i, nj - 3, k)].y_zeta_z = MidNodeInter4thLeft2(inter_temp);
				temp_eta[Idx(i, nj - 4, k)].y_zeta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, nj - 5 + iTemp, k))[2] * node->GetCoord(i, nj - 5 + iTemp, k)[0];
				}
				temp_eta[Idx(i, nj - 2, k)].z_zeta_x = MidNodeInter4thLeft1(inter_temp);
				temp_eta[Idx(i, nj - 3, k)].z_zeta_x = MidNodeInter4thLeft2(inter_temp);
				temp_eta[Idx(i, nj - 4, k)].z_zeta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, nj - 5 + iTemp, k))[2] * node->GetCoord(i, nj - 5 + iTemp, k)[1];
				}
				temp_eta[Idx(i, nj - 2, k)].x_zeta_y = MidNodeInter4thLeft1(inter_temp);
				temp_eta[Idx(i, nj - 3, k)].x_zeta_y = MidNodeInter4thLeft2(inter_temp);
				temp_eta[Idx(i, nj - 4, k)].x_zeta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, nj - 5 + iTemp, k))[0] * node->GetCoord(i, nj - 5 + iTemp, k)[2];
				}
				temp_eta[Idx(i, nj - 2, k)].y_xi_z = MidNodeInter4thLeft1(inter_temp);
				temp_eta[Idx(i, nj - 3, k)].y_xi_z = MidNodeInter4thLeft2(inter_temp);
				temp_eta[Idx(i, nj - 4, k)].y_xi_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, nj - 5 + iTemp, k))[0] * node->GetCoord(i, nj - 5 + iTemp, k)[0];
				}
				temp_eta[Idx(i, nj - 2, k)].z_xi_x = MidNodeInter4thLeft1(inter_temp);
				temp_eta[Idx(i, nj - 3, k)].z_xi_x = MidNodeInter4thLeft2(inter_temp);
				temp_eta[Idx(i, nj - 4, k)].z_xi_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, nj - 5 + iTemp, k))[0] * node->GetCoord(i, nj - 5 + iTemp, k)[1];
				}
				temp_eta[Idx(i, nj - 2, k)].x_xi_y = MidNodeInter4thLeft1(inter_temp);
				temp_eta[Idx(i, nj - 3, k)].x_xi_y = MidNodeInter4thLeft2(inter_temp);
				temp_eta[Idx(i, nj - 4, k)].x_xi_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, nj - 5 + iTemp, k))[2] * node->GetCoord(i, nj - 5 + iTemp, k)[1];
				}
				temp_eta[Idx(i, nj - 2, k)].z_zeta_y = MidNodeInter4thLeft1(inter_temp);
				temp_eta[Idx(i, nj - 3, k)].z_zeta_y = MidNodeInter4thLeft2(inter_temp);
				temp_eta[Idx(i, nj - 4, k)].z_zeta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, nj - 5 + iTemp, k))[2] * node->GetCoord(i, nj - 5 + iTemp, k)[2];
				}
				temp_eta[Idx(i, nj - 2, k)].x_zeta_z = MidNodeInter4thLeft1(inter_temp);
				temp_eta[Idx(i, nj - 3, k)].x_zeta_z = MidNodeInter4thLeft2(inter_temp);
				temp_eta[Idx(i, nj - 4, k)].x_zeta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, nj - 5 + iTemp, k))[2] * node->GetCoord(i, nj - 5 + iTemp, k)[0];
				}
				temp_eta[Idx(i, nj - 2, k)].y_zeta_x = MidNodeInter4thLeft1(inter_temp);
				temp_eta[Idx(i, nj - 3, k)].y_zeta_x = MidNodeInter4thLeft2(inter_temp);
				temp_eta[Idx(i, nj - 4, k)].y_zeta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, nj - 5 + iTemp, k))[0] * node->GetCoord(i, nj - 5 + iTemp, k)[1];
				}
				temp_eta[Idx(i, nj - 2, k)].z_xi_y = MidNodeInter4thLeft1(inter_temp);
				temp_eta[Idx(i, nj - 3, k)].z_xi_y = MidNodeInter4thLeft2(inter_temp);
				temp_eta[Idx(i, nj - 4, k)].z_xi_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, nj - 5 + iTemp, k))[0] * node->GetCoord(i, nj - 5 + iTemp, k)[2];
				}
				temp_eta[Idx(i, nj - 2, k)].x_xi_z = MidNodeInter4thLeft1(inter_temp);
				temp_eta[Idx(i, nj - 3, k)].x_xi_z = MidNodeInter4thLeft2(inter_temp);
				temp_eta[Idx(i, nj - 4, k)].x_xi_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, nj - 5 + iTemp, k))[0] * node->GetCoord(i, nj - 5 + iTemp, k)[0];
				}
				temp_eta[Idx(i, nj - 2, k)].y_xi_x = MidNodeInter4thLeft1(inter_temp);
				temp_eta[Idx(i, nj - 3, k)].y_xi_x = MidNodeInter4thLeft2(inter_temp);
				temp_eta[Idx(i, nj - 4, k)].y_xi_x = MidNodeInter4th(inter_temp);
			}
		}

		// 求出zeta方向的临时变量：逆变换度量系数乘以坐标
		struct TempZeta
		{
			double y_eta_z, z_eta_x, x_eta_y, y_xi_z, z_xi_x, x_xi_y, z_eta_y, x_eta_z, y_eta_x, z_xi_y, x_xi_z, y_xi_x;
		};
		dynamic_array<TempZeta> temp_zeta(ni * nj * nk);
		for (int k = 3; k < nk - 4; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetY(Idx(i, j, k + iTemp - 2))[1] * node->GetCoord(i, j, k + iTemp - 2)[2];
					}
					temp_zeta[Idx(i, j, k)].y_eta_z = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetZ(Idx(i, j, k + iTemp - 2))[1] * node->GetCoord(i, j, k + iTemp - 2)[0];
					}
					temp_zeta[Idx(i, j, k)].z_eta_x = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetX(Idx(i, j, k + iTemp - 2))[1] * node->GetCoord(i, j, k + iTemp - 2)[1];
					}
					temp_zeta[Idx(i, j, k)].x_eta_y = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetY(Idx(i, j, k + iTemp - 2))[0] * node->GetCoord(i, j, k + iTemp - 2)[2];
					}
					temp_zeta[Idx(i, j, k)].y_xi_z = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetZ(Idx(i, j, k + iTemp - 2))[0] * node->GetCoord(i, j, k + iTemp - 2)[0];
					}
					temp_zeta[Idx(i, j, k)].z_xi_x = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetX(Idx(i, j, k + iTemp - 2))[0] * node->GetCoord(i, j, k + iTemp - 2)[1];
					}
					temp_zeta[Idx(i, j, k)].x_xi_y = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetZ(Idx(i, j, k - 2 + iTemp))[1] * node->GetCoord(i, j, k - 2 + iTemp)[1];
					}
					temp_zeta[Idx(i, j, k)].z_eta_y = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetX(Idx(i, j, k - 2 + iTemp))[1] * node->GetCoord(i, j, k - 2 + iTemp)[2];
					}
					temp_zeta[Idx(i, j, k)].x_eta_z = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetY(Idx(i, j, k - 2 + iTemp))[1] * node->GetCoord(i, j, k - 2 + iTemp)[0];
					}
					temp_zeta[Idx(i, j, k)].y_eta_x = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetZ(Idx(i, j, k - 2 + iTemp))[0] * node->GetCoord(i, j, k - 2 + iTemp)[1];
					}
					temp_zeta[Idx(i, j, k)].z_xi_y = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetX(Idx(i, j, k - 2 + iTemp))[0] * node->GetCoord(i, j, k - 2 + iTemp)[2];
					}
					temp_zeta[Idx(i, j, k)].x_xi_z = MidNodeInter6th(inter_temp);
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						inter_temp[iTemp] =
							coef->GetY(Idx(i, j, k - 2 + iTemp))[0] * node->GetCoord(i, j, k - 2 + iTemp)[0];
					}
					temp_zeta[Idx(i, j, k)].y_xi_x = MidNodeInter6th(inter_temp);
				}
			}
		}
		// k=1/2,3/2,5/2;nk-3/2,nk-5/2,nk-7/2
		for (int j = 1; j < nj - 1; ++j)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				// k=1/2，3/2，5/2
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, j, 1 + iTemp))[1] * node->GetCoord(i, j, 1 + iTemp)[2];
				}
				temp_zeta[Idx(i, j, 0)].y_eta_z = MidNodeInter4thRight1(inter_temp);
				temp_zeta[Idx(i, j, 1)].y_eta_z = MidNodeInter4thRight2(inter_temp);
				temp_zeta[Idx(i, j, 2)].y_eta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, j, 1 + iTemp))[1] * node->GetCoord(i, j, 1 + iTemp)[0];
				}
				temp_zeta[Idx(i, j, 0)].z_eta_x = MidNodeInter4thRight1(inter_temp);
				temp_zeta[Idx(i, j, 1)].z_eta_x = MidNodeInter4thRight2(inter_temp);
				temp_zeta[Idx(i, j, 2)].z_eta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, j, 1 + iTemp))[1] * node->GetCoord(i, j, 1 + iTemp)[1];
				}
				temp_zeta[Idx(i, j, 0)].x_eta_y = MidNodeInter4thRight1(inter_temp);
				temp_zeta[Idx(i, j, 1)].x_eta_y = MidNodeInter4thRight2(inter_temp);
				temp_zeta[Idx(i, j, 2)].x_eta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, j, 1 + iTemp))[0] * node->GetCoord(i, j, 1 + iTemp)[2];
				}
				temp_zeta[Idx(i, j, 0)].y_xi_z = MidNodeInter4thRight1(inter_temp);
				temp_zeta[Idx(i, j, 1)].y_xi_z = MidNodeInter4thRight2(inter_temp);
				temp_zeta[Idx(i, j, 2)].y_xi_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, j, 1 + iTemp))[0] * node->GetCoord(i, j, 1 + iTemp)[0];
				}
				temp_zeta[Idx(i, j, 0)].z_xi_x = MidNodeInter4thRight1(inter_temp);
				temp_zeta[Idx(i, j, 1)].z_xi_x = MidNodeInter4thRight2(inter_temp);
				temp_zeta[Idx(i, j, 2)].z_xi_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, j, 1 + iTemp))[0] * node->GetCoord(i, j, 1 + iTemp)[1];
				}
				temp_zeta[Idx(i, j, 0)].x_xi_y = MidNodeInter4thRight1(inter_temp);
				temp_zeta[Idx(i, j, 1)].x_xi_y = MidNodeInter4thRight2(inter_temp);
				temp_zeta[Idx(i, j, 2)].x_xi_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, j, 1 + iTemp))[1] * node->GetCoord(i, j, 1 + iTemp)[1];
				}
				temp_zeta[Idx(i, j, 0)].z_eta_y = MidNodeInter4thRight1(inter_temp);
				temp_zeta[Idx(i, j, 1)].z_eta_y = MidNodeInter4thRight2(inter_temp);
				temp_zeta[Idx(i, j, 2)].z_eta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, j, 1 + iTemp))[1] * node->GetCoord(i, j, 1 + iTemp)[2];
				}
				temp_zeta[Idx(i, j, 0)].x_eta_z = MidNodeInter4thRight1(inter_temp);
				temp_zeta[Idx(i, j, 1)].x_eta_z = MidNodeInter4thRight2(inter_temp);
				temp_zeta[Idx(i, j, 2)].x_eta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, j, 1 + iTemp))[1] * node->GetCoord(i, j, 1 + iTemp)[0];
				}
				temp_zeta[Idx(i, j, 0)].y_eta_x = MidNodeInter4thRight1(inter_temp);
				temp_zeta[Idx(i, j, 1)].y_eta_x = MidNodeInter4thRight2(inter_temp);
				temp_zeta[Idx(i, j, 2)].y_eta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, j, 1 + iTemp))[0] * node->GetCoord(i, j, 1 + iTemp)[1];
				}
				temp_zeta[Idx(i, j, 0)].z_xi_y = MidNodeInter4thRight1(inter_temp);
				temp_zeta[Idx(i, j, 1)].z_xi_y = MidNodeInter4thRight2(inter_temp);
				temp_zeta[Idx(i, j, 2)].z_xi_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, j, 1 + iTemp))[0] * node->GetCoord(i, j, 1 + iTemp)[2];
				}
				temp_zeta[Idx(i, j, 0)].x_xi_z = MidNodeInter4thRight1(inter_temp);
				temp_zeta[Idx(i, j, 1)].x_xi_z = MidNodeInter4thRight2(inter_temp);
				temp_zeta[Idx(i, j, 2)].x_xi_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, j, 1 + iTemp))[0] * node->GetCoord(i, j, 1 + iTemp)[0];
				}
				temp_zeta[Idx(i, j, 0)].y_xi_x = MidNodeInter4thRight1(inter_temp);
				temp_zeta[Idx(i, j, 1)].y_xi_x = MidNodeInter4thRight2(inter_temp);
				temp_zeta[Idx(i, j, 2)].y_xi_x = MidNodeInter4th(inter_temp);
				// k=nk-3/2,nk-5/2,nk-7/2
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, j, nk - 5 + iTemp))[1] * node->GetCoord(i, j, nk - 5 + iTemp)[2];
				}
				temp_zeta[Idx(i, j, nk - 2)].y_eta_z = MidNodeInter4thLeft1(inter_temp);
				temp_zeta[Idx(i, j, nk - 3)].y_eta_z = MidNodeInter4thLeft2(inter_temp);
				temp_zeta[Idx(i, j, nk - 4)].y_eta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, j, nk - 5 + iTemp))[1] * node->GetCoord(i, j, nk - 5 + iTemp)[0];
				}
				temp_zeta[Idx(i, j, nk - 2)].z_eta_x = MidNodeInter4thLeft1(inter_temp);
				temp_zeta[Idx(i, j, nk - 3)].z_eta_x = MidNodeInter4thLeft2(inter_temp);
				temp_zeta[Idx(i, j, nk - 4)].z_eta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, j, nk - 5 + iTemp))[1] * node->GetCoord(i, j, nk - 5 + iTemp)[1];
				}
				temp_zeta[Idx(i, j, nk - 2)].x_eta_y = MidNodeInter4thLeft1(inter_temp);
				temp_zeta[Idx(i, j, nk - 3)].x_eta_y = MidNodeInter4thLeft2(inter_temp);
				temp_zeta[Idx(i, j, nk - 4)].x_eta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, j, nk - 5 + iTemp))[0] * node->GetCoord(i, j, nk - 5 + iTemp)[2];
				}
				temp_zeta[Idx(i, j, nk - 2)].y_xi_z = MidNodeInter4thLeft1(inter_temp);
				temp_zeta[Idx(i, j, nk - 3)].y_xi_z = MidNodeInter4thLeft2(inter_temp);
				temp_zeta[Idx(i, j, nk - 4)].y_xi_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, j, nk - 5 + iTemp))[0] * node->GetCoord(i, j, nk - 5 + iTemp)[0];
				}
				temp_zeta[Idx(i, j, nk - 2)].z_xi_x = MidNodeInter4thLeft1(inter_temp);
				temp_zeta[Idx(i, j, nk - 3)].z_xi_x = MidNodeInter4thLeft2(inter_temp);
				temp_zeta[Idx(i, j, nk - 4)].z_xi_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, j, nk - 5 + iTemp))[0] * node->GetCoord(i, j, nk - 5 + iTemp)[1];
				}
				temp_zeta[Idx(i, j, nk - 2)].x_xi_y = MidNodeInter4thLeft1(inter_temp);
				temp_zeta[Idx(i, j, nk - 3)].x_xi_y = MidNodeInter4thLeft2(inter_temp);
				temp_zeta[Idx(i, j, nk - 4)].x_xi_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, j, nk - 5 + iTemp))[1] * node->GetCoord(i, j, nk - 5 + iTemp)[1];
				}
				temp_zeta[Idx(i, j, nk - 2)].z_eta_y = MidNodeInter4thLeft1(inter_temp);
				temp_zeta[Idx(i, j, nk - 3)].z_eta_y = MidNodeInter4thLeft2(inter_temp);
				temp_zeta[Idx(i, j, nk - 4)].z_eta_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, j, nk - 5 + iTemp))[1] * node->GetCoord(i, j, nk - 5 + iTemp)[2];
				}
				temp_zeta[Idx(i, j, nk - 2)].x_eta_z = MidNodeInter4thLeft1(inter_temp);
				temp_zeta[Idx(i, j, nk - 3)].x_eta_z = MidNodeInter4thLeft2(inter_temp);
				temp_zeta[Idx(i, j, nk - 4)].x_eta_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, j, nk - 5 + iTemp))[1] * node->GetCoord(i, j, nk - 5 + iTemp)[0];
				}
				temp_zeta[Idx(i, j, nk - 2)].y_eta_x = MidNodeInter4thLeft1(inter_temp);
				temp_zeta[Idx(i, j, nk - 3)].y_eta_x = MidNodeInter4thLeft2(inter_temp);
				temp_zeta[Idx(i, j, nk - 4)].y_eta_x = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetZ(Idx(i, j, nk - 5 + iTemp))[0] * node->GetCoord(i, j, nk - 5 + iTemp)[1];
				}
				temp_zeta[Idx(i, j, nk - 2)].z_xi_y = MidNodeInter4thLeft1(inter_temp);
				temp_zeta[Idx(i, j, nk - 3)].z_xi_y = MidNodeInter4thLeft2(inter_temp);
				temp_zeta[Idx(i, j, nk - 4)].z_xi_y = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetX(Idx(i, j, nk - 5 + iTemp))[0] * node->GetCoord(i, j, nk - 5 + iTemp)[2];
				}
				temp_zeta[Idx(i, j, nk - 2)].x_xi_z = MidNodeInter4thLeft1(inter_temp);
				temp_zeta[Idx(i, j, nk - 3)].x_xi_z = MidNodeInter4thLeft2(inter_temp);
				temp_zeta[Idx(i, j, nk - 4)].x_xi_z = MidNodeInter4th(inter_temp);
				for (int iTemp = 0; iTemp < 4; iTemp++)
				{
					inter_temp[iTemp] = coef->GetY(Idx(i, j, nk - 5 + iTemp))[0] * node->GetCoord(i, j, nk - 5 + iTemp)[0];
				}
				temp_zeta[Idx(i, j, nk - 2)].y_xi_x = MidNodeInter4thLeft1(inter_temp);
				temp_zeta[Idx(i, j, nk - 3)].y_xi_x = MidNodeInter4thLeft2(inter_temp);
				temp_zeta[Idx(i, j, nk - 4)].y_xi_x = MidNodeInter4th(inter_temp);
			}
		}
		// 第五步：根据半点坐标和半点逆变换度量系数使用守恒形式计算整点度量系数（CMM1)
		// i direction
		double temp[6][6];
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 3; i < ni - 3; ++i)
				{
					auto coef_xi = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, k));
					auto coef_eta = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, k));
					auto coef_zeta = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, k));
					// i direction
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						int idx = Idx(i - 3 + iTemp, j, k);
						temp[0][iTemp] = 0.5 * (temp_xi[idx].y_zeta_z - temp_xi[idx].z_zeta_y); // y_zeta*z-z_zeta*y
						temp[1][iTemp] = 0.5 * (temp_xi[idx].z_zeta_x - temp_xi[idx].x_zeta_z); // z_zeta*x-x_zeta*z
						temp[2][iTemp] = 0.5 * (temp_xi[idx].x_zeta_y - temp_xi[idx].y_zeta_x); // x_zeta*y-y_zeta*x
						temp[3][iTemp] = 0.5 * (temp_xi[idx].z_eta_y - temp_xi[idx].y_eta_z);   // z_eta*y-y_eta*z
						temp[4][iTemp] = 0.5 * (temp_xi[idx].x_eta_z - temp_xi[idx].z_eta_x);   // x_eta*z-z_eta*x
						temp[5][iTemp] = 0.5 * (temp_xi[idx].y_eta_x - temp_xi[idx].x_eta_y);   // y_eta*x-x_eta*y
					}
					coef_eta[0] += NodeDifferece6th(temp[0]);
					coef_eta[1] += NodeDifferece6th(temp[1]);
					coef_eta[2] += NodeDifferece6th(temp[2]);
					coef_zeta[0] += NodeDifferece6th(temp[3]);
					coef_zeta[1] += NodeDifferece6th(temp[4]);
					coef_zeta[2] += NodeDifferece6th(temp[5]);
				}
			}
		}
		// j direction
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 3; j < nj - 3; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					auto coef_xi = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, k));
					auto coef_eta = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, k));
					auto coef_zeta = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						int idx = Idx(i, j - 3 + iTemp, k);
						temp[0][iTemp] = 0.5 * (temp_eta[idx].z_zeta_y - temp_eta[idx].y_zeta_z); // z_zeta*y-y_zeta*z
						temp[1][iTemp] = 0.5 * (temp_eta[idx].x_zeta_z - temp_eta[idx].z_zeta_x); // x_zeta*z-z_zeta*x
						temp[2][iTemp] = 0.5 * (temp_eta[idx].y_zeta_x - temp_eta[idx].x_zeta_y); //  y_zeta*x-x_zeta*y
						temp[3][iTemp] = 0.5 * (temp_eta[idx].y_xi_z - temp_eta[idx].z_xi_y);     // y_xi*z-z_xi*y
						temp[4][iTemp] = 0.5 * (temp_eta[idx].z_xi_x - temp_eta[idx].x_xi_z);     // z_xi*x-x_xi*z
						temp[5][iTemp] = 0.5 * (temp_eta[idx].x_xi_y - temp_eta[idx].y_xi_x);     // x_xi*y-y_xi*x
					}
					coef_xi[0] += NodeDifferece6th(temp[0]);
					coef_xi[1] += NodeDifferece6th(temp[1]);
					coef_xi[2] += NodeDifferece6th(temp[2]);
					coef_zeta[0] += NodeDifferece6th(temp[3]);
					coef_zeta[1] += NodeDifferece6th(temp[4]);
					coef_zeta[2] += NodeDifferece6th(temp[5]);
				}
			}
		}
		// k direction
		for (int k = 3; k < nk - 3; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					auto coef_xi = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, k));
					auto coef_eta = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, k));
					auto coef_zeta = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, k));
					for (int iTemp = 0; iTemp < 6; iTemp++)
					{
						int idx = Idx(i, j, k - 3 + iTemp);
						temp[0][iTemp] = 0.5 * (temp_zeta[idx].y_eta_z - temp_zeta[idx].z_eta_y); // y_eta*z-z_eta*y
						temp[1][iTemp] = 0.5 * (temp_zeta[idx].z_eta_x - temp_zeta[idx].x_eta_z); // z_eta*x-x_eta*z
						temp[2][iTemp] = 0.5 * (temp_zeta[idx].x_eta_y - temp_zeta[idx].y_eta_x); // x_eta*y-y_eta*x
						temp[3][iTemp] = 0.5 * (temp_zeta[idx].z_xi_y - temp_zeta[idx].y_xi_z);   // z_xi*y-y_xi*z
						temp[4][iTemp] = 0.5 * (temp_zeta[idx].x_xi_z - temp_zeta[idx].z_xi_x);   // x_xi*z-z_xi*x
						temp[5][iTemp] = 0.5 * (temp_zeta[idx].y_xi_x - temp_zeta[idx].x_xi_y);   //  y_xi*x-x_xi*y
					}
					coef_xi[0] += NodeDifferece6th(temp[0]);
					coef_xi[1] += NodeDifferece6th(temp[1]);
					coef_xi[2] += NodeDifferece6th(temp[2]);
					coef_eta[0] += NodeDifferece6th(temp[3]);
					coef_eta[1] += NodeDifferece6th(temp[4]);
					coef_eta[2] += NodeDifferece6th(temp[5]);
				}
			}
		}
		// 第六步：计算边界点度量系数
		// i=1,2;ni-3,ni-2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				auto coef_xi1 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(1, j, k));
				auto coef_eta1 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(1, j, k));
				auto coef_zeta1 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(1, j, k));
				auto coef_xi2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(2, j, k));
				auto coef_eta2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(2, j, k));
				auto coef_zeta2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(2, j, k));
				for (int iTemp = 0; iTemp < 5; iTemp++)
				{
					int idx = Idx(0 + iTemp, j, k);
					temp[0][iTemp] = 0.5 * (temp_xi[idx].y_zeta_z - temp_xi[idx].z_zeta_y); // y_zeta*z-z_zeta*y
					temp[1][iTemp] = 0.5 * (temp_xi[idx].z_zeta_x - temp_xi[idx].x_zeta_z); // z_zeta*x-x_zeta*z
					temp[2][iTemp] = 0.5 * (temp_xi[idx].x_zeta_y - temp_xi[idx].y_zeta_x); // x_zeta*y-y_zeta*x
					temp[3][iTemp] = 0.5 * (temp_xi[idx].z_eta_y - temp_xi[idx].y_eta_z);   // z_eta*y-y_eta*z
					temp[4][iTemp] = 0.5 * (temp_xi[idx].x_eta_z - temp_xi[idx].z_eta_x);   // x_eta*z-z_eta*x
					temp[5][iTemp] = 0.5 * (temp_xi[idx].y_eta_x - temp_xi[idx].x_eta_y);   // y_eta*x-x_eta*y
				}
				coef_eta1[0] += NodeDifferece4thRight(temp[0]);
				coef_eta1[1] += NodeDifferece4thRight(temp[1]);
				coef_eta1[2] += NodeDifferece4thRight(temp[2]);
				coef_zeta1[0] += NodeDifferece4thRight(temp[3]);
				coef_zeta1[1] += NodeDifferece4thRight(temp[4]);
				coef_zeta1[2] += NodeDifferece4thRight(temp[5]);
				coef_eta2[0] += NodeDifferece4th(temp[0]);
				coef_eta2[1] += NodeDifferece4th(temp[1]);
				coef_eta2[2] += NodeDifferece4th(temp[2]);
				coef_zeta2[0] += NodeDifferece4th(temp[3]);
				coef_zeta2[1] += NodeDifferece4th(temp[4]);
				coef_zeta2[2] = NodeDifferece4th(temp[5]);
				auto coef_xi_n2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(ni - 2, j, k));
				auto coef_eta_n2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(ni - 2, j, k));
				auto coef_zeta_n2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(ni - 2, j, k));
				auto coef_xi_n3 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(ni - 3, j, k));
				auto coef_eta_n3 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(ni - 3, j, k));
				auto coef_zeta_n3 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(ni - 3, j, k));
				for (int iTemp = 0; iTemp < 5; iTemp++)
				{
					int idx = Idx(ni - 2 - iTemp, j, k);
					temp[0][iTemp] = 0.5 * (temp_xi[idx].y_zeta_z - temp_xi[idx].z_zeta_y); // y_zeta*z-z_zeta*y
					temp[1][iTemp] = 0.5 * (temp_xi[idx].z_zeta_x - temp_xi[idx].x_zeta_z); // z_zeta*x-x_zeta*z
					temp[2][iTemp] = 0.5 * (temp_xi[idx].x_zeta_y - temp_xi[idx].y_zeta_x); // x_zeta*y-y_zeta*x
					temp[3][iTemp] = 0.5 * (temp_xi[idx].z_eta_y - temp_xi[idx].y_eta_z);   // z_eta*y-y_eta*z
					temp[4][iTemp] = 0.5 * (temp_xi[idx].x_eta_z - temp_xi[idx].z_eta_x);   // x_eta*z-z_eta*x
					temp[5][iTemp] = 0.5 * (temp_xi[idx].y_eta_x - temp_xi[idx].x_eta_y);   // y_eta*x-x_eta*y
				}
				coef_eta_n2[0] += NodeDifferece4thLeft(temp[0]);
				coef_eta_n2[1] += NodeDifferece4thLeft(temp[1]);
				coef_eta_n2[2] += NodeDifferece4thLeft(temp[2]);
				coef_zeta_n2[0] += NodeDifferece4thLeft(temp[3]);
				coef_zeta_n2[1] += NodeDifferece4thLeft(temp[4]);
				coef_zeta_n2[2] += NodeDifferece4thLeft(temp[5]);
				coef_eta_n3[0] += -NodeDifferece4th(temp[0]);
				coef_eta_n3[1] += -NodeDifferece4th(temp[1]);
				coef_eta_n3[2] += -NodeDifferece4th(temp[2]);
				coef_zeta_n3[0] += -NodeDifferece4th(temp[3]);
				coef_zeta_n3[1] += -NodeDifferece4th(temp[4]);
				coef_zeta_n3[2] += -NodeDifferece4th(temp[5]);
			}
		}
		// j=1,2;nj-3,nj-2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				auto coef_xi1 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, 1, k));
				auto coef_eta1 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, 1, k));
				auto coef_zeta1 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, 1, k));
				auto coef_xi2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, 2, k));
				auto coef_eta2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, 2, k));
				auto coef_zeta2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, 2, k));
				for (int iTemp = 0; iTemp < 5; iTemp++)
				{
					int idx = Idx(i, 0 + iTemp, k);
					temp[0][iTemp] = 0.5 * (temp_eta[idx].z_zeta_y - temp_eta[idx].y_zeta_z); // z_zeta*y-y_zeta*z
					temp[1][iTemp] = 0.5 * (temp_eta[idx].x_zeta_z - temp_eta[idx].z_zeta_x); // x_zeta*z-z_zeta*x
					temp[2][iTemp] = 0.5 * (temp_eta[idx].y_zeta_x - temp_eta[idx].x_zeta_y); //  y_zeta*x-x_zeta*y
					temp[3][iTemp] = 0.5 * (temp_eta[idx].y_xi_z - temp_eta[idx].z_xi_y);     // y_xi*z-z_xi*y
					temp[4][iTemp] = 0.5 * (temp_eta[idx].z_xi_x - temp_eta[idx].x_xi_z);     // z_xi*x-x_xi*z
					temp[5][iTemp] = 0.5 * (temp_eta[idx].x_xi_y - temp_eta[idx].y_xi_x);     // x_xi*y-y_xi*x
				}
				coef_xi1[0] += NodeDifferece4thRight(temp[0]);
				coef_xi1[1] += NodeDifferece4thRight(temp[1]);
				coef_xi1[2] += NodeDifferece4thRight(temp[2]);
				coef_zeta1[0] += NodeDifferece4thRight(temp[3]);
				coef_zeta1[1] += NodeDifferece4thRight(temp[4]);
				coef_zeta1[2] += NodeDifferece4thRight(temp[5]);
				coef_xi2[0] += NodeDifferece4th(temp[0]);
				coef_xi2[1] += NodeDifferece4th(temp[1]);
				coef_xi2[2] += NodeDifferece4th(temp[2]);
				coef_zeta2[0] += NodeDifferece4th(temp[3]);
				coef_zeta2[1] += NodeDifferece4th(temp[4]);
				coef_zeta2[2] += NodeDifferece4th(temp[5]);
				auto coef_xi_n2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, nj - 2, k));
				auto coef_eta_n2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, nj - 2, k));
				auto coef_zeta_n2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, nj - 2, k));
				auto coef_xi_n3 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, nj - 3, k));
				auto coef_eta_n3 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, nj - 3, k));
				auto coef_zeta_n3 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, nj - 3, k));
				for (int iTemp = 0; iTemp < 5; iTemp++)
				{
					int idx = Idx(i, nj - 2 - iTemp, k);
					temp[0][iTemp] = 0.5 * (temp_eta[idx].z_zeta_y - temp_eta[idx].y_zeta_z); // z_zeta*y-y_zeta*z
					temp[1][iTemp] = 0.5 * (temp_eta[idx].x_zeta_z - temp_eta[idx].z_zeta_x); // x_zeta*z-z_zeta*x
					temp[2][iTemp] = 0.5 * (temp_eta[idx].y_zeta_x - temp_eta[idx].x_zeta_y); //  y_zeta*x-x_zeta*y
					temp[3][iTemp] = 0.5 * (temp_eta[idx].y_xi_z - temp_eta[idx].z_xi_y);     // y_xi*z-z_xi*y
					temp[4][iTemp] = 0.5 * (temp_eta[idx].z_xi_x - temp_eta[idx].x_xi_z);     // z_xi*x-x_xi*z
					temp[5][iTemp] = 0.5 * (temp_eta[idx].x_xi_y - temp_eta[idx].y_xi_x);     // x_xi*y-y_xi*x
				}
				coef_xi_n2[0] += NodeDifferece4thLeft(temp[0]);
				coef_xi_n2[1] += NodeDifferece4thLeft(temp[1]);
				coef_xi_n2[2] += NodeDifferece4thLeft(temp[2]);
				coef_zeta_n2[0] += NodeDifferece4thLeft(temp[3]);
				coef_zeta_n2[1] += NodeDifferece4thLeft(temp[4]);
				coef_zeta_n2[2] += NodeDifferece4thLeft(temp[5]);
				coef_xi_n3[0] += -NodeDifferece4th(temp[0]);
				coef_xi_n3[1] += -NodeDifferece4th(temp[1]);
				coef_xi_n3[2] += -NodeDifferece4th(temp[2]);
				coef_zeta_n3[0] += -NodeDifferece4th(temp[3]);
				coef_zeta_n3[1] += -NodeDifferece4th(temp[4]);
				coef_zeta_n3[2] += -NodeDifferece4th(temp[5]);
			}
		}
		// k=1,2;nk-3,nk-2
		for (int j = 1; j < nj - 1; ++j)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				auto coef_xi1 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, 1));
				auto coef_eta1 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, 1));
				auto coef_zeta1 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, 1));
				auto coef_xi2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, 2));
				auto coef_eta2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, 2));
				auto coef_zeta2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, 2));
				for (int iTemp = 0; iTemp < 5; iTemp++)
				{
					int idx = Idx(i, j, 0 + iTemp);
					temp[0][iTemp] = 0.5 * (temp_zeta[idx].y_eta_z - temp_zeta[idx].z_eta_y); // y_eta*z-z_eta*y
					temp[1][iTemp] = 0.5 * (temp_zeta[idx].z_eta_x - temp_zeta[idx].x_eta_z); // z_eta*x-x_eta*z
					temp[2][iTemp] = 0.5 * (temp_zeta[idx].x_eta_y - temp_zeta[idx].y_eta_x); // x_eta*y-y_eta*x
					temp[3][iTemp] = 0.5 * (temp_zeta[idx].z_xi_y - temp_zeta[idx].y_xi_z);   // z_xi*y-y_xi*z
					temp[4][iTemp] = 0.5 * (temp_zeta[idx].x_xi_z - temp_zeta[idx].z_xi_x);   // x_xi*z-z_xi*x
					temp[5][iTemp] = 0.5 * (temp_zeta[idx].y_xi_x - temp_zeta[idx].x_xi_y);   //  y_xi*x-x_xi*y
				}
				coef_xi1[0] += NodeDifferece4thRight(temp[0]);
				coef_xi1[1] += NodeDifferece4thRight(temp[1]);
				coef_xi1[2] += NodeDifferece4thRight(temp[2]);
				coef_eta1[0] += NodeDifferece4thRight(temp[3]);
				coef_eta1[1] += NodeDifferece4thRight(temp[4]);
				coef_eta1[2] += NodeDifferece4thRight(temp[5]);
				coef_xi2[0] += NodeDifferece4th(temp[0]);
				coef_xi2[1] += NodeDifferece4th(temp[1]);
				coef_xi2[2] += NodeDifferece4th(temp[2]);
				coef_eta2[0] += NodeDifferece4th(temp[3]);
				coef_eta2[1] += NodeDifferece4th(temp[4]);
				coef_eta2[2] += NodeDifferece4th(temp[5]);
				auto coef_xi_n2 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, nk - 2));
				auto coef_eta_n2 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, nk - 2));
				auto coef_zeta_n2 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, nk - 2));
				auto coef_xi_n3 = GetNodeMetrics()->GetXi(idx_proxy->GetIdx(i, j, nk - 3));
				auto coef_eta_n3 = GetNodeMetrics()->GetEta(idx_proxy->GetIdx(i, j, nk - 3));
				auto coef_zeta_n3 = GetNodeMetrics()->GetZeta(idx_proxy->GetIdx(i, j, nk - 3));
				for (int iTemp = 0; iTemp < 5; iTemp++)
				{
					int idx = Idx(i, j, nk - 2 - iTemp);
					temp[0][iTemp] = 0.5 * (temp_zeta[idx].y_eta_z - temp_zeta[idx].z_eta_y); // y_eta*z-z_eta*y
					temp[1][iTemp] = 0.5 * (temp_zeta[idx].z_eta_x - temp_zeta[idx].x_eta_z); // z_eta*x-x_eta*z
					temp[2][iTemp] = 0.5 * (temp_zeta[idx].x_eta_y - temp_zeta[idx].y_eta_x); // x_eta*y-y_eta*x
					temp[3][iTemp] = 0.5 * (temp_zeta[idx].z_xi_y - temp_zeta[idx].y_xi_z);   // z_xi*y-y_xi*z
					temp[4][iTemp] = 0.5 * (temp_zeta[idx].x_xi_z - temp_zeta[idx].z_xi_x);   // x_xi*z-z_xi*x
					temp[5][iTemp] = 0.5 * (temp_zeta[idx].y_xi_x - temp_zeta[idx].x_xi_y);   //  y_xi*x-x_xi*y
				}
				coef_xi_n2[0] += NodeDifferece4thLeft(temp[0]);
				coef_xi_n2[1] += NodeDifferece4thLeft(temp[1]);
				coef_xi_n2[2] += NodeDifferece4thLeft(temp[2]);
				coef_eta_n2[0] += NodeDifferece4thLeft(temp[3]);
				coef_eta_n2[1] += NodeDifferece4thLeft(temp[4]);
				coef_eta_n2[2] += NodeDifferece4thLeft(temp[5]);
				coef_xi_n3[0] += -NodeDifferece4th(temp[0]);
				coef_xi_n3[1] += -NodeDifferece4th(temp[1]);
				coef_xi_n3[2] += -NodeDifferece4th(temp[2]);
				coef_eta_n3[0] += -NodeDifferece4th(temp[3]);
				coef_eta_n3[1] += -NodeDifferece4th(temp[4]);
				coef_eta_n3[2] += -NodeDifferece4th(temp[5]);
			}
		}
	}
	void NSSolverStruct::CalcMetricsS4()
	{
		// TODO
	}
	void NSSolverStruct::Preprocess()
	{
		NSSolver::Preprocess();
		auto para = GetPara();
		if (para->GetInterSchme() == InterpolationScheme::Grad)
		{
			CalcPrimGradWLS();
			CalcLimiterVK();
		}
	}
	void NSSolverStruct::CalcPrimGradWLS()
	{
		auto grid = GetGrid();
		if (grid->GetDim() == TWO_DIM)
		{
			CalcPrimGradWLS2D();
		}
		else if (grid->GetDim() == THREE_DIM)
		{
			CalcPrimGradWLS3D();
		}
	}
	void NSSolverStruct::CalcPrimGradWLS2D()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto idx_proxy = GetIdxProxy();
		auto data_manager = GetDataManager();
		auto ni = grid->GetNi();
		auto nj = grid->GetNj();
		auto nk = grid->GetNk();
		Matrix2d A, A_inv;
		Eigen::Vector2d b, grad;
		int neighbor_num = 4;
		dynamic_array<int> temp_i, temp_j, temp_k;
		temp_i.resize(neighbor_num);
		temp_j.resize(neighbor_num);
		temp_k.resize(neighbor_num);

		double delta_x, delta_y, delta_z, delta_val;
		double weight;
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					A.setZero();
					int idx = idx_proxy->GetIdx(i, j, k);
					auto coord = node->GetCoord(i, j, k);
					temp_i[0] = i - 1, temp_j[0] = j, temp_k[0] = k;
					temp_i[1] = i + 1, temp_j[1] = j, temp_k[1] = k;
					temp_i[2] = i, temp_j[2] = j - 1, temp_k[2] = k;
					temp_i[3] = i, temp_j[3] = j + 1, temp_k[3] = k;

					for (int iNeigh = 0; iNeigh < neighbor_num; ++iNeigh)
					{
						auto coord_neigh = node->GetCoord(temp_i[iNeigh], temp_j[iNeigh], temp_k[iNeigh]);
						delta_x = coord_neigh[0] - coord[0];
						delta_y = coord_neigh[1] - coord[1];
						weight = 1.0 / sqrt(delta_x * delta_x + delta_y * delta_y);
						A(0, 0) += weight * delta_x * delta_x;
						A(0, 1) += weight * delta_x * delta_y;
						A(1, 0) += weight * delta_y * delta_x;
						A(1, 1) += weight * delta_y * delta_y;
					}
					A_inv = A.inverse();
					for (int idx_eq = 0; idx_eq < 5; ++idx_eq)
					{
						b.setZero();
						for (int iNeigh = 0; iNeigh < neighbor_num; ++iNeigh)
						{
							delta_val = data_manager->GetPrim(
								idx_eq, idx_proxy->GetIdx(temp_i[iNeigh], temp_j[iNeigh], temp_k[iNeigh])) -
								data_manager->GetPrim(idx_eq, idx);
							auto coord_neigh = node->GetCoord(temp_i[iNeigh], temp_j[iNeigh], temp_k[iNeigh]);
							delta_x = coord_neigh[0] - coord[0];
							delta_y = coord_neigh[1] - coord[1];
							weight = 1.0 / sqrt(delta_x * delta_x + delta_y * delta_y);
							b[0] += weight * delta_x * delta_val;
							b[1] += weight * delta_y * delta_val;
						}

						grad = A_inv * b;
						for (size_t iDim = 0; iDim < 2; ++iDim)
						{
							data_manager->SetPrimitiveGrad(idx_eq, iDim, idx, grad(iDim));
						}
						data_manager->SetPrimitiveGrad(idx_eq, 2, idx, 0);
					}
				}
			}
		}
	}
	void NSSolverStruct::CalcPrimGradWLS3D()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto idx_proxy = GetIdxProxy();
		auto data_manager = GetDataManager();
		auto ni = grid->GetNi();
		auto nj = grid->GetNj();
		auto nk = grid->GetNk();
		Matrix3d A, A_inv;
		Eigen::Vector3d b, grad;
		int neighbor_num = 6;
		dynamic_array<int> temp_i, temp_j, temp_k;

		temp_i.resize(neighbor_num);
		temp_j.resize(neighbor_num);
		temp_k.resize(neighbor_num);

		double delta_x, delta_y, delta_z, delta_val;
		double weight;
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					A.setZero();
					int idx = idx_proxy->GetIdx(i, j, k);
					auto coord = node->GetCoord(i, j, k);
					temp_i[0] = i - 1, temp_j[0] = j, temp_k[0] = k;
					temp_i[1] = i + 1, temp_j[1] = j, temp_k[1] = k;
					temp_i[2] = i, temp_j[2] = j - 1, temp_k[2] = k;
					temp_i[3] = i, temp_j[3] = j + 1, temp_k[3] = k;
					temp_i[4] = i, temp_j[4] = j, temp_k[4] = k - 1;
					temp_i[5] = i, temp_j[5] = j, temp_k[5] = k + 1;

					for (int iNeigh = 0; iNeigh < neighbor_num; ++iNeigh)
					{
						auto coord_neigh = node->GetCoord(temp_i[iNeigh], temp_j[iNeigh], temp_k[iNeigh]);
						delta_x = coord_neigh[0] - coord[0];
						delta_y = coord_neigh[1] - coord[1];
						delta_z = coord_neigh[2] - coord[2];
						weight = 1.0 / sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
						A(0, 0) += weight * delta_x * delta_x;
						A(0, 1) += weight * delta_x * delta_y;
						A(0, 2) += weight * delta_x * delta_z;
						A(1, 0) += weight * delta_y * delta_x;
						A(1, 1) += weight * delta_y * delta_y;
						A(1, 2) += weight * delta_y * delta_z;
						A(2, 0) += weight * delta_z * delta_x;
						A(2, 1) += weight * delta_z * delta_y;
						A(2, 2) += weight * delta_z * delta_z;
					}
					A_inv = A.inverse();
					for (int idx_eq = 0; idx_eq < 5; ++idx_eq)
					{
						b.setZero();
						for (int iNeigh = 0; iNeigh < neighbor_num; ++iNeigh)
						{
							delta_val = data_manager->GetPrim(
								idx_eq, idx_proxy->GetIdx(temp_i[iNeigh], temp_j[iNeigh], temp_k[iNeigh])) -
								data_manager->GetPrim(idx_eq, idx);
							auto coord_neigh = node->GetCoord(temp_i[iNeigh], temp_j[iNeigh], temp_k[iNeigh]);
							delta_x = coord_neigh[0] - coord[0];
							delta_y = coord_neigh[1] - coord[1];
							delta_z = coord_neigh[2] - coord[2];
							weight = 1.0 / sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
							b[0] += weight * delta_x * delta_val;
							b[1] += weight * delta_y * delta_val;
							b[2] += weight * delta_z * delta_val;
						}

						grad = A_inv * b;
						for (size_t iDim = 0; iDim < 3; ++iDim)
						{
							data_manager->SetPrimitiveGrad(idx_eq, iDim, idx, grad(iDim));
							// data_manager->SetPrimitiveGrad(idx_eq, iDim, idx, 0.0);
						}
					}
					// Log::info("density grad: {} {} {}", data_manager->GetPrimGrad(0, 0, idx),
					// data_manager->GetPrimGrad(0, 1, idx), data_manager->GetPrimGrad(0, 2, idx));
				}
			}
		}
	}
	void NSSolverStruct::CalcLimiterVK()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto data_manager = GetDataManager();
		auto idx_proxy = GetIdxProxy();
		double eps = 1e-6;
		double vk_coef = 1.0e-5;
		auto ni = grid->GetNi();
		auto nj = grid->GetNj();
		auto nk = grid->GetNk();
		double max_val, min_val;
		int neighbor_num = 4;
		dynamic_array<int> temp_i, temp_j, temp_k;
		if (grid->GetDim() == TWO_DIM)
		{
			temp_i.resize(neighbor_num);
			temp_j.resize(neighbor_num);
			temp_k.resize(neighbor_num);
		}
		else
		{
			neighbor_num = 6;
			temp_i.resize(neighbor_num);
			temp_j.resize(neighbor_num);
			temp_k.resize(neighbor_num);
		}
		for (int idx_eq = 0; idx_eq < 5; idx_eq++)
		{

			for (int k = 1; k < nk - 1; ++k)
			{
				for (int j = 1; j < nj - 1; ++j)
				{
					for (int i = 1; i < ni - 1; ++i)
					{
						int idx = idx_proxy->GetIdx(i, j, k);
						auto coord = node->GetCoord(i, j, k);
						temp_i[0] = i - 1, temp_j[0] = j, temp_k[0] = k;
						temp_i[1] = i + 1, temp_j[1] = j, temp_k[1] = k;
						temp_i[2] = i, temp_j[2] = j - 1, temp_k[2] = k;
						temp_i[3] = i, temp_j[3] = j + 1, temp_k[3] = k;
						if (grid->GetDim() == THREE_DIM)
						{
							temp_i[4] = i, temp_j[4] = j, temp_k[4] = k - 1;
							temp_i[5] = i, temp_j[5] = j, temp_k[5] = k + 1;
						}
						max_val = data_manager->GetPrim(idx_eq, idx);
						min_val = data_manager->GetPrim(idx_eq, idx);
						for (int iNeigh = 0; iNeigh < neighbor_num; ++iNeigh)
						{
							int idx_neigh = idx_proxy->GetIdx(temp_i[iNeigh], temp_j[iNeigh], temp_k[iNeigh]);
							max_val = Max(max_val, data_manager->GetPrim(idx_eq, idx_neigh));
							min_val = Min(min_val, data_manager->GetPrim(idx_eq, idx_neigh));
						}
						eps = vk_coef * (max_val - min_val) + SMALL_NUMBER;
						eps = eps * eps;
						double delta_max = max_val - data_manager->GetPrim(idx_eq, idx);
						double delta_min = min_val - data_manager->GetPrim(idx_eq, idx);
						double temp_coef = LARGE_NUMBER;
						for (int iNeigh = 0; iNeigh < neighbor_num; ++iNeigh)
						{
							int idx_neigh = idx_proxy->GetIdx(temp_i[iNeigh], temp_j[iNeigh], temp_k[iNeigh]);
							auto coord_neigh = node->GetCoord(temp_i[iNeigh], temp_j[iNeigh], temp_k[iNeigh]);
							double delta2 = 0.0;
							for (int iDim = 0; iDim < 3; ++iDim)
							{
								delta2 += (coord_neigh[iDim] - coord[iDim]) * data_manager->GetPrimGrad(idx_eq, iDim, idx);
							}
							delta2 *= 0.5;
							if (delta2 > 0.0)
							{
								temp_coef =
									Min(temp_coef, LimiterVK(max_val - data_manager->GetPrim(idx_eq, idx), delta2, eps));
							}
							else if (delta2 < 0.0)
							{
								temp_coef =
									Min(temp_coef, LimiterVK(min_val - data_manager->GetPrim(idx_eq, idx), delta2, eps));
							}
							else
							{
								temp_coef = Min(temp_coef, 1.0);
							}
						}
						data_manager->SetLimiter(idx_eq, idx, temp_coef);
						// data_manager->SetLimiter(idx_eq, idx, 0.0);
					}
				}
			}
		}
	}
	void NSSolverStruct::CalcTimeStepLocal()
	{
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		int is, ie, js, je, ks, ke;
		grid->GetRange(is, ie, js, je, ks, ke);
		auto para = GetPara();
		double cfl = para->GetCflNumber();
		double gamma = GetGas()->GetGamma();
		auto density = data_manager->GetPrim(ID_DENSITY);
		auto pressure = data_manager->GetPrim(ID_PRESSURE);
		auto velocity_x = data_manager->GetPrim(ID_VELOCITY_X);
		auto velocity_y = data_manager->GetPrim(ID_VELOCITY_Y);
		auto velocity_z = data_manager->GetPrim(ID_VELOCITY_Z);

#pragma omp parallel for
		for (int k = ks; k <= ke; ++k)
		{
			for (int j = js; j <= je; ++j)
			{
				for (int i = is; i <= ie; ++i)
				{
					int idx = m_idx_proxy->GetIdx(i, j, k);
					auto xi = m_node_metrics->GetXi(idx);
					auto eta = m_node_metrics->GetEta(idx);
					auto zeta = m_node_metrics->GetZeta(idx);
					auto jacobi = m_node_metrics->GetJacobian(idx);
					double c = sqrt(gamma * pressure[idx] / density[idx]);
					double norm_xi = sqrt(xi[0] * xi[0] + xi[1] * xi[1] + xi[2] * xi[2]);
					double norm_eta = sqrt(eta[0] * eta[0] + eta[1] * eta[1] + eta[2] * eta[2]);
					double norm_zeta = sqrt(zeta[0] * zeta[0] + zeta[1] * zeta[1] + zeta[2] * zeta[2]);
					double u_xi = velocity_x[idx] * xi[0] + velocity_y[idx] * xi[1] + velocity_z[idx] * xi[2];
					double u_eta = velocity_x[idx] * eta[0] + velocity_y[idx] * eta[1] + velocity_z[idx] * eta[2];
					double u_zeta = velocity_x[idx] * zeta[0] + velocity_y[idx] * zeta[1] + velocity_z[idx] * zeta[2];
					double lamda = abs(u_xi) + abs(u_eta) + c * (norm_xi + norm_eta);
					if (grid->GetDim() == THREE_DIM)
					{
						lamda += abs(u_zeta) + c * norm_zeta;
					}
					lamda = lamda * jacobi;
					data_manager->SetTimeStep(idx, cfl / lamda);
					// data_manager->SetTimeStep(idx, 0.0002);

				}
			}
		}
	}
	void NSSolverStruct::CalcMinTimeStep(double& dt)
	{
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		int is, ie, js, je, ks, ke;
		grid->GetRange(is, ie, js, je, ks, ke);
		auto para = GetPara();
		double cfl = para->GetCflNumber();
		double gamma = GetGas()->GetGamma();
		double min_dt = LARGE_NUMBER;
#pragma omp parallel for reduction(min : min_dt)
		for (int k = ks; k <= ke; ++k)
		{
			for (int j = js; j <= je; ++j)
			{
				for (int i = is; i <= ie; ++i)
				{
					int idx = m_idx_proxy->GetIdx(i, j, k);
					if (data_manager->GetTimeStep(idx) < min_dt)
					{
						min_dt = data_manager->GetTimeStep(idx);
					}
				}
			}
		}
		dt = min_dt;
	}
	void NSSolverStruct::ReduceTimeStep(double& dt)
	{
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		auto ni = grid->GetNi();
		auto nj = grid->GetNj();
		auto nk = grid->GetNk();
		for (int k = 0; k < nk; ++k)
		{
#pragma omp parallel for
			for (int j = 0; j < nj; ++j)
			{
				for (int i = 0; i < ni; ++i)
				{
					int idx = m_idx_proxy->GetIdx(i, j, k);
					data_manager->SetTimeStep(idx, dt);
				}
			}
		}
	}
	void NSSolverStruct::RungeKutta()
	{
		auto grid = GetGrid();
		auto para = GetPara();
		auto data_manager = GetDataManager();
		int is, ie, js, je, ks, ke;
		grid->GetRange(is, ie, js, je, ks, ke);
		int rk_step = para->GetRkStep();
		for (int k = ks; k <= ke; ++k)
		{
#pragma omp parallel for
			for (int j = js; j <= je; ++j)
			{
				for (int i = is; i <= ie; ++i)
				{
					int idx = m_idx_proxy->GetIdx(i, j, k);
					for (int idx_eq = 0; idx_eq < 5; ++idx_eq)
					{
						data_manager->SetConsOld(idx_eq, idx, data_manager->GetCons(idx_eq, idx));
					}
				}
			}
		}
		for (int iStep = 0; iStep < rk_step; ++iStep)
		{
			BoundaryCondition();
			CalcResidual();
			auto& rk_coef = para->GetRkCoef(iStep);
			for (int k = ks; k <= ke; ++k)
			{
#pragma omp parallel for
				for (int j = js; j <= je; ++j)
				{
					for (int i = is; i <= ie; ++i)
					{
						int idx = m_idx_proxy->GetIdx(i, j, k);
						double dt = data_manager->GetTimeStep(idx);
						double jacobi = m_node_metrics->GetJacobian(idx);
						for (int idx_eq = 0; idx_eq < 5; ++idx_eq)
						{
							data_manager->SetCons(idx_eq, idx,
								rk_coef[0] * data_manager->GetConsOld(idx_eq, idx) +
								rk_coef[1] * data_manager->GetCons(idx_eq, idx) +
								rk_coef[2] * dt * jacobi * data_manager->GetResidual(idx_eq, idx));
						}
					}
				}
			}
			UpdateField();
		}
	}
	void NSSolverStruct::Prim2Cons()
	{
		auto gas = GetGas();
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		auto ni = grid->GetNi();
		auto nj = grid->GetNj();
		auto nk = grid->GetNk();
		double prim[5], cons[5];
		for (int k = 0; k < nk; ++k)
		{
#pragma omp parallel for private(prim, cons)
			for (int j = 0; j < nj; ++j)
			{
				for (int i = 0; i < ni; ++i)
				{
					int idx = m_idx_proxy->GetIdx(i, j, k);
					for (int idx_eq = 0; idx_eq < GetPara()->GetEqNum(); ++idx_eq)
					{
						prim[idx_eq] = data_manager->GetPrim(idx_eq, idx);
					}
					gas->Prim2Cons(prim, cons);
					data_manager->SetCons(idx, cons);
				}
			}
		}
	}
	void NSSolverStruct::Cons2Prim()
	{
		auto gas = GetGas();
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		auto ni = grid->GetNi();
		auto nj = grid->GetNj();
		auto nk = grid->GetNk();
		double prim[5], cons[5];
		for (int k = 0; k < nk; ++k)
		{
#pragma omp parallel for private(prim, cons)
			for (int j = 0; j < nj; ++j)
			{
				for (int i = 0; i < ni; ++i)
				{
					int idx = m_idx_proxy->GetIdx(i, j, k);
					for (int idx_eq = 0; idx_eq < GetPara()->GetEqNum(); ++idx_eq)
					{
						cons[idx_eq] = data_manager->GetCons(idx_eq, idx);
					}
					gas->Cons2Prim(cons, prim);
					data_manager->SetPrim(idx, prim);
				}
			}
		}
	}
	void NSSolverStruct::ZeroResidual()
	{
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		auto ni = grid->GetNi();
		auto nj = grid->GetNj();
		auto nk = grid->GetNk();
		double res[5] = { 0.0 };
		for (int k = 0; k < nk; ++k)
		{
#pragma omp parallel for
			for (int j = 0; j < nj; ++j)
			{
				for (int i = 0; i < ni; ++i)
				{
					int idx = m_idx_proxy->GetIdx(i, j, k);
					data_manager->SetResidual(idx, res);
				}
			}
		}
	}

	void NSSolverStruct::InterMidNodePrim_1st()
	{
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		auto ni = grid->GetNi();
		auto nj = grid->GetNj();
		auto nk = grid->GetNk();
		double value[5];
		double value_left[5], value_right[5];
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					int idx = m_idx_proxy->GetIdx(i, j, k);
					int idx_right = m_idx_proxy->GetIdx(i + 1, j, k);
					for (int idx_eq = 0; idx_eq < 5; ++idx_eq)
					{
						data_manager->SetMidNodePrim(idx_eq, 0, idx, data_manager->GetPrim(idx_eq, idx),
							data_manager->GetPrim(idx_eq, idx_right));
					}
					idx_right = m_idx_proxy->GetIdx(i, j + 1, k);
					for (int idx_eq = 0; idx_eq < 5; ++idx_eq)
					{
						data_manager->SetMidNodePrim(idx_eq, 1, idx, data_manager->GetPrim(idx_eq, idx),
							data_manager->GetPrim(idx_eq, idx_right));
					}
					if (grid->GetDim() == THREE_DIM)
					{
						idx_right = m_idx_proxy->GetIdx(i, j, k + 1);
						for (int idx_eq = 0; idx_eq < 5; ++idx_eq)
						{
							data_manager->SetMidNodePrim(idx_eq, 2, idx, data_manager->GetPrim(idx_eq, idx),
								data_manager->GetPrim(idx_eq, idx_right));
						}
					}
				}
			}
		}
	}
	void NSSolverStruct::InterMidNodePrim_Grad()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto data_manager = GetDataManager();
		auto ni = grid->GetNi();
		auto nj = grid->GetNj();
		auto nk = grid->GetNk();
		double value[5];
		double coord_vec[3];
		double value_left[5], value_right[5];
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					int idx = m_idx_proxy->GetIdx(i, j, k);
					int idx_right = m_idx_proxy->GetIdx(i + 1, j, k);
					coord_vec[0] = node->GetCoord(i + 1, j, k)[0] - node->GetCoord(i, j, k)[0];
					coord_vec[1] = node->GetCoord(i + 1, j, k)[1] - node->GetCoord(i, j, k)[1];
					coord_vec[2] = node->GetCoord(i + 1, j, k)[2] - node->GetCoord(i, j, k)[2];
					InterMidNodePrim_Grad(idx, idx_right, coord_vec, value_left, value_right);
					for (int idx_eq = 0; idx_eq < 5; ++idx_eq)
					{
						data_manager->SetMidNodePrim(idx_eq, 0, idx, value_left[idx_eq], value_right[idx_eq]);
					}

					idx_right = m_idx_proxy->GetIdx(i, j + 1, k);
					coord_vec[0] = node->GetCoord(i, j + 1, k)[0] - node->GetCoord(i, j, k)[0];
					coord_vec[1] = node->GetCoord(i, j + 1, k)[1] - node->GetCoord(i, j, k)[1];
					coord_vec[2] = node->GetCoord(i, j + 1, k)[2] - node->GetCoord(i, j, k)[2];
					InterMidNodePrim_Grad(idx, idx_right, coord_vec, value_left, value_right);
					for (int idx_eq = 0; idx_eq < 5; ++idx_eq)
					{
						data_manager->SetMidNodePrim(idx_eq, 1, idx, value_left[idx_eq], value_right[idx_eq]);
					}
					if (grid->GetDim() == THREE_DIM)
					{
						idx_right = m_idx_proxy->GetIdx(i, j, k + 1);
						coord_vec[0] = node->GetCoord(i, j, k + 1)[0] - node->GetCoord(i, j, k)[0];
						coord_vec[1] = node->GetCoord(i, j, k + 1)[1] - node->GetCoord(i, j, k)[1];
						coord_vec[2] = node->GetCoord(i, j, k + 1)[2] - node->GetCoord(i, j, k)[2];
						InterMidNodePrim_Grad(idx, idx_right, coord_vec, value_left, value_right);
						for (int idx_eq = 0; idx_eq < 5; ++idx_eq)
						{
							data_manager->SetMidNodePrim(idx_eq, 2, idx, value_left[idx_eq], value_right[idx_eq]);
						}
					}
				}
			}
		}
	}
	void NSSolverStruct::InterMidNodePrim_MUSCL()
	{
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		int equ_num = GetPara()->GetEqNum();
		auto ni = grid->GetNi();
		auto nj = grid->GetNj();
		auto nk = grid->GetNk();
		// MUSCL整点变量(i-2,i-1,i,i+1,i+2)的值
		double value_temp[5];
		int i_temp[5], j_temp[5], k_temp[5];
		int is, ie, js, je, ks, ke;
		grid->GetRange(is, ie, js, je, ks, ke);
		double left_value, right_value;
		int direction[3][3] = { {1,0,0},{0,1,0},{0,0,1} };
		///@note
		/// 从ks-1开始，到ke结束，因为对第一个计算点ks进行通量差分时，需要使用ks-1/2处的值
		for (int k = ks; k <= ke; ++k)
		{
#pragma omp parallel for private(value_temp, i_temp, j_temp, k_temp, left_value, right_value)
			for (int j = js; j <= je; ++j)
			{
				for (int i = is; i <= ie; ++i)
				{
					for (dimension_type iDim = 0; iDim < grid->GetDim(); ++iDim)
					{
						for (int iTemp = 0; iTemp < 5; ++iTemp)
						{
							i_temp[iTemp] = i + direction[iDim][0] * (iTemp - 2);
							j_temp[iTemp] = j + direction[iDim][1] * (iTemp - 2);
							k_temp[iTemp] = k + direction[iDim][2] * (iTemp - 2);
						}
						auto idx = m_idx_proxy->GetIdx(i, j, k);
						for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
						{
							for (int iTemp = 0; iTemp < 5; ++iTemp)
							{
								value_temp[iTemp] = data_manager->GetPrim(
									idx_eq, m_idx_proxy->GetIdx(i_temp[iTemp], j_temp[iTemp], k_temp[iTemp]));
							}
							InterMidNodePrim_MUSCL(value_temp, left_value, right_value);
							data_manager->SetMidNodePrim(idx_eq, iDim, idx, left_value, right_value);
						}
					}
				}
			}
		}
	}

	void NSSolverStruct::InterMidNodePrim_MUSCLSV()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto data_manager = GetDataManager();
		int equ_num = GetPara()->GetEqNum();
		auto ni = grid->GetNi();
		auto nj = grid->GetNj();
		auto nk = grid->GetNk();
		double** coord_temp = new double* [5];
		for (int i = 0; i < 5; ++i)
		{
			coord_temp[i] = new double[3];
		}
		// MUSCL整点变量(i-2,i-1,i,i+1,i+2)的值
		double value_temp[5];
		int idx_temp[5], i_temp[5], j_temp[5], k_temp[5];
		int is, ie, js, je, ks, ke;
		grid->GetRange(is, ie, js, je, ks, ke);
		double left_value, right_value;
		int direction[3][3] = { {1,0,0},{0,1,0},{0,0,1} };

		///@note
		/// 从ks-1开始，到ke结束，因为对第一个计算点ks进行通量差分时，需要使用ks-1/2处的值
		for (int k = ks; k <= ke; ++k)
		{
			for (int j = js; j <= je; ++j)
			{
				for (int i = is; i <= ie; ++i)
				{
					for (dimension_type iDim = 0; iDim < grid->GetDim(); ++iDim)
					{
						for (int iTemp = 0; iTemp < 5; ++iTemp)
						{
							i_temp[iTemp] = i + direction[iDim][0] * (iTemp - 2);
							j_temp[iTemp] = j + direction[iDim][1] * (iTemp - 2);
							k_temp[iTemp] = k + direction[iDim][2] * (iTemp - 2);
						}
						auto idx = m_idx_proxy->GetIdx(i, j, k);
						for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
						{
							for (int iTemp = 0; iTemp < 5; ++iTemp)
							{
								value_temp[iTemp] = data_manager->GetPrim(
									idx_eq, m_idx_proxy->GetIdx(i_temp[iTemp], j_temp[iTemp], k_temp[iTemp]));
								for (int jDim = 0; jDim < 3; ++jDim)
								{
									coord_temp[iTemp][jDim] = node->GetCoord(i_temp[iTemp], j_temp[iTemp], k_temp[iTemp])[jDim];
								}
							}
							InterMidNodePrim_MUSCLSV(value_temp, coord_temp, left_value, right_value);
							data_manager->SetMidNodePrim(idx_eq, iDim, idx, left_value, right_value);
						}
					}
				}
			}
		}
		for (int i = 0; i < 5; ++i)
		{
			delete[] coord_temp[i];
		}
		delete[] coord_temp;
	}

	void NSSolverStruct::InterGhostMidNode_MUSCL()
	{
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		int equ_num = GetPara()->GetEqNum();
		auto ni = grid->GetNi();
		auto nj = grid->GetNj();
		auto nk = grid->GetNk();
		double value[5];
		int is, ie, js, je, ks, ke;
		grid->GetRange(is, ie, js, je, ks, ke);
		int idx_temp[5];
		double left_value, right_value;
		for (int k = ks; k <= ke; ++k)
		{
#pragma omp parallel for private(value, left_value, right_value, idx_temp)
			for (int j = js; j <= je; ++j)
			{
				// i=is-1处的值
				// 不存在i-2处的值且不参与计算,用i-1处的值代替，不影响计算结果
				idx_temp[0] = m_idx_proxy->GetIdx(is - 3, j, k);
				idx_temp[1] = m_idx_proxy->GetIdx(is - 2, j, k);
				idx_temp[2] = m_idx_proxy->GetIdx(is - 1, j, k);
				idx_temp[3] = m_idx_proxy->GetIdx(is, j, k);
				idx_temp[4] = m_idx_proxy->GetIdx(is + 1, j, k);
				for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
				{
					value[0] = data_manager->GetPrim(idx_eq, idx_temp[0]);
					value[1] = data_manager->GetPrim(idx_eq, idx_temp[1]);
					value[2] = data_manager->GetPrim(idx_eq, idx_temp[2]);
					value[3] = data_manager->GetPrim(idx_eq, idx_temp[3]);
					value[4] = data_manager->GetPrim(idx_eq, idx_temp[4]);
					InterMidNodePrim_MUSCL(value, left_value, right_value);
					data_manager->SetMidNodePrim(idx_eq, 0, idx_temp[2], left_value, right_value);
				}
			}
		}
		for (int k = ks; k <= ke; ++k)
		{
#pragma omp parallel for private(value, left_value, right_value, idx_temp)
			for (int i = is; i <= ie; ++i)
			{
				// j=js-1处的值，(js-1/2)右值会用到，但(js-1/2)左值在计算中不会使用
				// 不存在j-2处的值,用j-1处的值代替，不影响计算结果
				idx_temp[0] = m_idx_proxy->GetIdx(i, js - 3, k);
				idx_temp[1] = m_idx_proxy->GetIdx(i, js - 2, k);
				idx_temp[2] = m_idx_proxy->GetIdx(i, js - 1, k);
				idx_temp[3] = m_idx_proxy->GetIdx(i, js, k);
				idx_temp[4] = m_idx_proxy->GetIdx(i, js + 1, k);
				for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
				{
					value[0] = data_manager->GetPrim(idx_eq, idx_temp[0]);
					value[1] = data_manager->GetPrim(idx_eq, idx_temp[1]);
					value[2] = data_manager->GetPrim(idx_eq, idx_temp[2]);
					value[3] = data_manager->GetPrim(idx_eq, idx_temp[3]);
					value[4] = data_manager->GetPrim(idx_eq, idx_temp[4]);
					InterMidNodePrim_MUSCL(value, left_value, right_value);
					data_manager->SetMidNodePrim(idx_eq, 1, idx_temp[2], left_value, right_value);
				}
			}
		}

		if (grid->GetDim() == THREE_DIM)
		{
#pragma omp parallel for private(value, left_value, right_value, idx_temp)
			for (int j = js; j <= je; ++j)
			{
				for (int i = is; i <= ie; ++i)
				{
					// k=ks-1处的值，(ks-1/2)右值会用到，但(ks-1/2)左值在计算中不会使用
					// 不存在k-2处的值,用k-1处的值代替，不影响计算结果
					idx_temp[0] = m_idx_proxy->GetIdx(i, j, ks - 2);
					idx_temp[1] = m_idx_proxy->GetIdx(i, j, ks - 2);
					idx_temp[2] = m_idx_proxy->GetIdx(i, j, ks - 1);
					idx_temp[3] = m_idx_proxy->GetIdx(i, j, ks);
					idx_temp[4] = m_idx_proxy->GetIdx(i, j, ks + 1);
					for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
					{
						value[0] = data_manager->GetPrim(idx_eq, idx_temp[0]);
						value[1] = data_manager->GetPrim(idx_eq, idx_temp[1]);
						value[2] = data_manager->GetPrim(idx_eq, idx_temp[2]);
						value[3] = data_manager->GetPrim(idx_eq, idx_temp[3]);
						value[4] = data_manager->GetPrim(idx_eq, idx_temp[4]);
						InterMidNodePrim_MUSCL(value, left_value, right_value);
						data_manager->SetMidNodePrim(idx_eq, 2, idx_temp[2], left_value, right_value);
					}
				}
			}
		}
	}

	void NSSolverStruct::CalcPrimGhostMidNodeMUSCL_SV()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto data_manager = GetDataManager();
		int equ_num = GetPara()->GetEqNum();
		auto ni = grid->GetNi();
		auto nj = grid->GetNj();
		auto nk = grid->GetNk();
		double value_temp[5];
		double value_left[5], value_right[5];
		int is, ie, js, je, ks, ke;
		grid->GetRange(is, ie, js, je, ks, ke);
		double** coord_temp = new double* [5];
		for (int i = 0; i < 5; ++i)
		{
			coord_temp[i] = new double[3];
		}
		int idx_temp[5], i_temp[5], j_temp[5], k_temp[5];
		double left_value, right_value;
		int i, j, k;
		for (k = ks; k <= ke; ++k)
		{
			for (j = js; j <= je; ++j)
			{
				// i=is-1处的值
				// 不存在i-2处的值且不参与计算,用i-1处的值代替，不影响计算结果
				i_temp[0] = is - 2, i_temp[1] = is - 2, i_temp[2] = is - 1, i_temp[3] = is, i_temp[4] = is + 1;
				j_temp[0] = j_temp[1] = j_temp[2] = j_temp[3] = j_temp[4] = j;
				k_temp[0] = k_temp[1] = k_temp[2] = k_temp[3] = k_temp[4] = k;

				auto idx = m_idx_proxy->GetIdx(is - 1, j, k);
				for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
				{
					for (int iTemp = 0; iTemp < 5; ++iTemp)
					{
						value_temp[iTemp] =
							data_manager->GetPrim(idx_eq, m_idx_proxy->GetIdx(i_temp[iTemp], j_temp[iTemp], k_temp[iTemp]));
						for (int iDim = 0; iDim < 3; ++iDim)
						{
							coord_temp[iTemp][iDim] = node->GetCoord(i_temp[iTemp], j_temp[iTemp], k_temp[iTemp])[iDim];
						}
					}
					// MidNodeMUSCL(value_temp, left_value, right_value);
					InterMidNodePrim_MUSCLSV(value_temp, coord_temp, left_value, right_value);
					data_manager->SetMidNodePrim(idx_eq, 0, idx, left_value, right_value);
				}
			}
		}

		for (k = ks; k <= ke; ++k)
		{
			for (i = is; i <= ie; ++i)
			{
				// j=js-1处的值，(js-1/2)右值会用到，但(js-1/2)左值在计算中不会使用
				// 不存在j-2处的值,用j-1处的值代替，不影响计算结果
				i_temp[0] = i_temp[1] = i_temp[2] = i_temp[3] = i_temp[4] = i;
				j_temp[0] = js - 2, j_temp[1] = js - 2, j_temp[2] = js - 1, j_temp[3] = js, j_temp[4] = js + 1;
				k_temp[0] = k_temp[1] = k_temp[2] = k_temp[3] = k_temp[4] = k;

				auto idx = m_idx_proxy->GetIdx(i, js - 1, k);
				for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
				{
					for (int iTemp = 0; iTemp < 5; ++iTemp)
					{
						value_temp[iTemp] =
							data_manager->GetPrim(idx_eq, m_idx_proxy->GetIdx(i_temp[iTemp], j_temp[iTemp], k_temp[iTemp]));
						for (int iDim = 0; iDim < 3; ++iDim)
						{
							coord_temp[iTemp][iDim] = node->GetCoord(i_temp[iTemp], j_temp[iTemp], k_temp[iTemp])[iDim];
						}
					}
					// MidNodeMUSCL(value_temp, left_value, right_value);
					InterMidNodePrim_MUSCLSV(value_temp, coord_temp, left_value, right_value);
					data_manager->SetMidNodePrim(idx_eq, 1, idx, left_value, right_value);
				}
			}
		}

		if (grid->GetDim() == THREE_DIM)
		{
			for (j = js; j <= je; ++j)
			{
				for (i = is; i <= ie; ++i)
				{
					// k=ks-1处的值，(ks-1/2)右值会用到，但(ks-1/2)左值在计算中不会使用
					// 不存在k-2处的值,用k-1处的值代替，不影响计算结果
					i_temp[0] = i_temp[1] = i_temp[2] = i_temp[3] = i_temp[4] = i;
					j_temp[0] = j_temp[1] = j_temp[2] = j_temp[3] = j_temp[4] = j;
					k_temp[0] = ks - 2, k_temp[1] = ks - 2, k_temp[2] = ks - 1, k_temp[3] = ks, k_temp[4] = ks + 1;

					auto idx = m_idx_proxy->GetIdx(i, j, ks - 1);
					for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
					{
						for (int iTemp = 0; iTemp < 5; ++iTemp)
						{
							value_temp[iTemp] = data_manager->GetPrim(
								idx_eq, m_idx_proxy->GetIdx(i_temp[iTemp], j_temp[iTemp], k_temp[iTemp]));
							for (int iDim = 0; iDim < 3; ++iDim)
							{
								coord_temp[iTemp][iDim] = node->GetCoord(i_temp[iTemp], j_temp[iTemp], k_temp[iTemp])[iDim];
							}
						}
						// MidNodeMUSCL(value_temp, left_value, right_value);
						InterMidNodePrim_MUSCLSV(value_temp, coord_temp, left_value, right_value);
						data_manager->SetMidNodePrim(idx_eq, 2, idx, left_value, right_value);
					}
				}
			}
		}
		delete[] coord_temp;
	}

	void NSSolverStruct::InterMidNodePrim_WCNS5()
	{
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		int equ_num = GetPara()->GetEqNum();
		auto ni = grid->GetNi();
		auto nj = grid->GetNj();
		auto nk = grid->GetNk();
		// WCNS整点变量(i-2,i-1,i,i+1,i+2)的值
		double value_temp[5];
		int idx_temp[5];
		int is, ie, js, je, ks, ke;
		grid->GetRange(is, ie, js, je, ks, ke);
		double left_value[5]{}, right_value[5]{};
		double min_density = SMALL_NUMBER;
		double max_density = LARGE_NUMBER;
		double min_pressure = SMALL_NUMBER;
		double max_pressure = LARGE_NUMBER;
		int direction[3][3] = { {1,0,0},{0,1,0},{0,0,1} };
		/// 从ks-1开始，到ke结束，因为对第一个计算点ks进行通量差分时，需要使用ks-1/2处的值
		for (int k = ks - 1; k <= ke + 1; ++k)
		{
#pragma omp parallel for private(idx_temp, value_temp, left_value, right_value)
			for (int j = js - 1; j <= je + 1; ++j)
			{
				for (int i = is - 1; i <= ie + 1; ++i)
				{
					for (int iDim = 0; iDim < grid->GetDim(); ++iDim)
					{
						for (int iTemp = 0; iTemp < 5; ++iTemp)
						{
							idx_temp[iTemp] = m_idx_proxy->GetIdx(i + direction[iDim][0] * (iTemp - 2),
								j + direction[iDim][1] * (iTemp - 2), k + direction[iDim][2] * (iTemp - 2));
						}
						double value_base[5];
						value_base[0] = data_manager->GetPrim(0, idx_temp[2]);
						value_base[4] = data_manager->GetPrim(4, idx_temp[2]);
						value_base[1] = value_base[2] = value_base[3] = sqrt(value_base[4] / value_base[0]);
						for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
						{
							for (int iTemp = 0; iTemp < 5; ++iTemp)
							{
								value_temp[iTemp] = data_manager->GetPrim(idx_eq, idx_temp[iTemp]) / value_base[idx_eq];
							}
							InterMidNodePrim_WCNS5(value_temp, left_value[idx_eq], right_value[idx_eq]);
						}
						//if (left_value[0] < min_density || left_value[4] < min_pressure || left_value[0]>max_density || left_value[4]>max_pressure)
						//{
						//	for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
						//	{
						//		left_value[idx_eq] = data_manager->GetPrim(idx_eq, idx_temp[2]);
						//	}
						//}
						//if (right_value[0] < min_density || right_value[4] < min_pressure || right_value[0]>max_density || right_value[4]>max_pressure)
						//{
						//	for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
						//	{
						//		right_value[idx_eq] = data_manager->GetPrim(idx_eq, idx_temp[2]);
						//	}
						//}
						auto idx = m_idx_proxy->GetIdx(i, j, k);
						for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
						{
							data_manager->SetMidNodePrimLeft(idx_eq, iDim, idx_temp[2], left_value[idx_eq] * value_base[idx_eq]);
							data_manager->SetMidNodePrimRight(idx_eq, iDim, idx_temp[1], right_value[idx_eq] * value_base[idx_eq]);
						}
					}
				}
			}
		}
	}

	void NSSolverStruct::InterGhostMidNodePrim_WCNS5()
	{
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		int equ_num = GetPara()->GetEqNum();
		auto ni = grid->GetNi();
		auto nj = grid->GetNj();
		auto nk = grid->GetNk();
		double value_left[5]{}, value_right[5]{};
		int is, ie, js, je, ks, ke;
		grid->GetRange(is, ie, js, je, ks, ke);
		// 左侧的模板点(0,1,2,3,4),右侧的模板点(N-5,N-4,N-3,N-2,N-1)
		int idx_temp_left[5], idx_temp_right[5];
		double left_value, right_value;
		int i, j, k;
#pragma omp parallel for private(j, k, idx_temp_left, idx_temp_right, value_left, value_right, left_value,right_value)
		for (k = ks; k <= ke; ++k)
		{
			for (j = js; j <= je; ++j)
			{
				for (int iTemp = 0; iTemp < 5; ++iTemp)
				{
					idx_temp_left[iTemp] = m_idx_proxy->GetIdx(iTemp, j, k);
					idx_temp_right[iTemp] = m_idx_proxy->GetIdx(ni - 5 + iTemp, j, k);
				}
				for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
				{
					for (int iTemp = 0; iTemp < 5; ++iTemp)
					{
						value_left[iTemp] = data_manager->GetPrim(idx_eq, idx_temp_left[iTemp]);
						value_right[iTemp] = data_manager->GetPrim(idx_eq, idx_temp_right[iTemp]);
					}
					// 1/2
					left_value = right_value =
						(5 * value_left[0] + 15 * value_left[1] - 5 * value_left[2] + value_left[3]) / 16;
					data_manager->SetMidNodePrim(idx_eq, 0, idx_temp_left[0], left_value, right_value);
					// 3/2
					left_value = right_value =
						(-value_left[0] + 9 * value_left[1] + 9 * value_left[2] - value_left[3]) / 16;
					data_manager->SetMidNodePrim(idx_eq, 0, idx_temp_left[1], left_value, right_value);
					// N-5/2
					left_value = right_value =
						(-value_right[4] + 9 * value_right[3] + 9 * value_right[2] - value_right[1]) / 16;
					data_manager->SetMidNodePrim(idx_eq, 0, idx_temp_right[2], left_value, right_value);
					// N-3/2
					left_value = right_value =
						(5 * value_right[4] + 15 * value_right[3] - 5 * value_right[2] + value_right[1]) / 16;
					data_manager->SetMidNodePrim(idx_eq, 0, idx_temp_right[3], left_value, right_value);
				}
			}
		}
#pragma omp parallel for private(i, k, idx_temp_left, idx_temp_right, value_left, value_right, left_value,right_value)
		for (k = ks; k <= ke; ++k)
		{
			for (i = is; i <= ie; ++i)
			{
				// j=js-1处的值，(js-1/2)右值会用到，但(js-1/2)左值在计算中不会使用
				// 不存在j-2处的值,用j-1处的值代替，不影响计算结果
				for (int iTemp = 0; iTemp < 5; ++iTemp)
				{
					idx_temp_left[iTemp] = m_idx_proxy->GetIdx(i, iTemp, k);
					idx_temp_right[iTemp] = m_idx_proxy->GetIdx(i, nj - 5 + iTemp, k);
				}
				for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
				{
					for (int iTemp = 0; iTemp < 5; ++iTemp)
					{
						value_left[iTemp] = data_manager->GetPrim(idx_eq, idx_temp_left[iTemp]);
						value_right[iTemp] = data_manager->GetPrim(idx_eq, idx_temp_right[iTemp]);
					}
					// 1/2
					left_value = right_value =
						(5 * value_left[0] + 15 * value_left[1] - 5 * value_left[2] + value_left[3]) / 16;
					data_manager->SetMidNodePrim(idx_eq, 1, idx_temp_left[0], left_value, right_value);
					// 3/2
					left_value = right_value =
						(-value_left[0] + 9 * value_left[1] + 9 * value_left[2] - value_left[3]) / 16;
					data_manager->SetMidNodePrim(idx_eq, 1, idx_temp_left[1], left_value, right_value);
					// N-5/2
					left_value = right_value =
						(-value_right[4] + 9 * value_right[3] + 9 * value_right[2] - value_right[1]) / 16;
					data_manager->SetMidNodePrim(idx_eq, 1, idx_temp_right[2], left_value, right_value);
					// N-3/2
					left_value = right_value =
						(5 * value_right[4] + 15 * value_right[3] - 5 * value_right[2] + value_right[1]) / 16;
					data_manager->SetMidNodePrim(idx_eq, 1, idx_temp_right[3], left_value, right_value);
				}
			}
		}

		if (grid->GetDim() == THREE_DIM)
		{
#pragma omp parallel for private(i, j, idx_temp_left, idx_temp_right, value_left, value_right, left_value, right_value)
			for (j = js; j <= je; ++j)
			{
				for (i = is; i <= ie; ++i)
				{
					// k=ks-1处的值，(ks-1/2)右值会用到，但(ks-1/2)左值在计算中不会使用
					// 不存在k-2处的值,用k-1处的值代替，不影响计算结果
					for (int iTemp = 0; iTemp < 5; ++iTemp)
					{
						idx_temp_left[iTemp] = m_idx_proxy->GetIdx(i, j, iTemp);
						idx_temp_right[iTemp] = m_idx_proxy->GetIdx(i, j, nk - 5 + iTemp);
					}
					for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
					{
						for (int iTemp = 0; iTemp < 5; ++iTemp)
						{
							value_left[iTemp] = data_manager->GetPrim(idx_eq, idx_temp_left[iTemp]);
							value_right[iTemp] = data_manager->GetPrim(idx_eq, idx_temp_right[iTemp]);
						}
						// 1/2
						left_value = right_value =
							(5 * value_left[0] + 15 * value_left[1] - 5 * value_left[2] + value_left[3]) / 16;
						data_manager->SetMidNodePrim(idx_eq, 2, idx_temp_left[0], left_value, right_value);
						// 3/2
						left_value = right_value =
							(-value_left[0] + 9 * value_left[1] + 9 * value_left[2] - value_left[3]) / 16;
						data_manager->SetMidNodePrim(idx_eq, 2, idx_temp_left[1], left_value, right_value);
						// N-5/2
						left_value = right_value =
							(-value_right[4] + 9 * value_right[3] + 9 * value_right[2] - value_right[1]) / 16;
						data_manager->SetMidNodePrim(idx_eq, 2, idx_temp_right[2], left_value, right_value);
						// N-3/2
						left_value = right_value =
							(5 * value_right[4] + 15 * value_right[3] - 5 * value_right[2] + value_right[1]) / 16;
						data_manager->SetMidNodePrim(idx_eq, 2, idx_temp_right[3], left_value, right_value);
					}
				}
			}
		}
	}

	void NSSolverStruct::InterMidNodePrim_Grad(int idx_left, int idx_right, double* coord_vec, double* value_left,
		double* value_right)
	{
		int equ_num = GetPara()->GetEqNum();
		auto data_manager = GetDataManager();
		for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
		{
			value_left[idx_eq] = data_manager->GetPrim(idx_eq, idx_left) +
				0.0 * 0.5 * data_manager->GetLimiter(idx_eq, idx_left) *
				(coord_vec[0] * data_manager->GetPrimGrad(idx_eq, 0, idx_left) +
					coord_vec[1] * data_manager->GetPrimGrad(idx_eq, 1, idx_left) +
					coord_vec[2] * data_manager->GetPrimGrad(idx_eq, 2, idx_left));
			value_right[idx_eq] = data_manager->GetPrim(idx_eq, idx_right) -
				0.0 * 0.5 * data_manager->GetLimiter(idx_eq, idx_right) *
				(coord_vec[0] * data_manager->GetPrimGrad(idx_eq, 0, idx_right) +
					coord_vec[1] * data_manager->GetPrimGrad(idx_eq, 1, idx_right) +
					coord_vec[2] * data_manager->GetPrimGrad(idx_eq, 2, idx_right));
		}
	}
	void NSSolverStruct::InterMidNodePrim_1st(int index_left, int index_right, double* value_rec_left,
		double* value_rec_right)
	{
		int equ_num = GetPara()->GetEqNum();
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
		{
			value_rec_left[idx_eq] = data_manager->GetPrim(idx_eq, index_left);
			value_rec_right[idx_eq] = data_manager->GetPrim(idx_eq, index_right);
		}
	}

	void NSSolverStruct::InterMidNodePrim_MUSCL(const double* value, double& value_left, double& value_right)
	{
		double k = -1.0;
		double delta_plus = value[3] - value[2];
		double delta_minus = value[2] - value[1];
		double delta = LimiterVanLeer(delta_plus, delta_minus);
		value_left = value[2] + 0.125 * ((1.0 - k) * delta + (1.0 + k) * delta);
		//value_left = value[2] + 0.25 * ((1.0 - k) * delta_minus + (1.0 + k) * delta_plus);
		delta_plus = value[4] - value[3];
		delta_minus = value[3] - value[2];
		delta = LimiterVanLeer(delta_plus, delta_minus);
		value_right = value[3] - 0.125 * ((1.0 - k) * delta + (1.0 + k) * delta);
		//value_right = value[3] - 0.25 * ((1.0 - k) * delta_plus + (1.0 + k) * delta_minus);
		if (std::isnan(value_left) || std::isnan(value_right))
		{
			Log::error("MUSCL interpolation failed!");
			Log::error("value[0] = {}, value[1] = {}, value[2] = {}, value[3] = {}, "
				"value[4] = {}",
				value[0], value[1], value[2], value[3], value[4]);
			exit(0);
		}
	}

	void NSSolverStruct::InterMidNodePrim_MUSCLSV(const double* value, double** coord, double& value_left,
		double& value_right)
	{
		double k = -1.0;
		double l =
			sqrt(pow(coord[3][0] - coord[2][0], 2) + pow(coord[3][1] - coord[2][1], 2) + pow(coord[3][2] - coord[2][2], 2));
		double delta_plus = value[3] - value[2];
		double l_plus =
			sqrt(pow(coord[3][0] - coord[2][0], 2) + pow(coord[3][1] - coord[2][1], 2) + pow(coord[3][2] - coord[2][2], 2));
		double delta_minus = value[2] - value[1];
		double l_minus =
			sqrt(pow(coord[2][0] - coord[1][0], 2) + pow(coord[2][1] - coord[1][1], 2) + pow(coord[2][2] - coord[1][2], 2));
		l_plus = l / l_plus;
		l_minus = l / l_minus;
		// double delta = LimiterVanLeer(delta_plus, delta_minus);
		// value_left = value[2] + 0.125 * ((1.0 - k) * delta + (1.0 + k) * delta);
		value_left = value[2] + 0.25 * ((1.0 - k) * l_minus * delta_minus + (1.0 + k) * l_plus * delta_plus);
		delta_plus = value[4] - value[3];
		l_plus =
			sqrt(pow(coord[4][0] - coord[3][0], 2) + pow(coord[4][1] - coord[3][1], 2) + pow(coord[4][2] - coord[3][2], 2));
		delta_minus = value[3] - value[2];
		l_minus =
			sqrt(pow(coord[3][0] - coord[2][0], 2) + pow(coord[3][1] - coord[2][1], 2) + pow(coord[3][2] - coord[2][2], 2));
		l_plus = l / l_plus;
		l_minus = l / l_minus;
		// delta = LimiterVanLeer(delta_plus, delta_minus);
		// value_right = value[3] - 0.125 * ((1.0 - k) * delta + (1.0 + k) * delta);
		value_right = value[3] - 0.25 * ((1.0 - k) * l_plus * delta_plus + (1.0 + k) * l_minus * delta_minus);
		if (std::isnan(value_left) || std::isnan(value_right))
		{
			Log::error("MUSCL interpolation failed!");
			Log::error("value[0] = {}, value[1] = {}, value[2] = {}, value[3] = {}, "
				"value[4] = {}",
				value[0], value[1], value[2], value[3], value[4]);
			exit(0);
		}
	}

	// 参考《计算空气动力学》pp.84
	void NSSolverStruct::InterMidNodePrim_WCNS5(double* value, double& value_left, double& value_right)
	{
		//double ref_value = value[2];
		//if (abs(ref_value) > SMALL_NUMBER)
		//{
		//	for (int i = 0; i < 5; ++i)
		//	{
		//		value[i] = value[i] / ref_value;
		//	}
		//}
		//1st order derivative
		double f[3];
		f[0] = 0.5 * (value[0] - 4.0 * value[1] + 3.0 * value[2]);
		f[1] = 0.5 * (value[3] - value[1]);
		f[2] = 0.5 * (-3.0 * value[2] + 4.0 * value[3] - value[4]);
		//2nd order derivative
		double s[3];
		s[0] = (value[0] - 2.0 * value[1] + value[2]);
		s[1] = (value[1] - 2.0 * value[2] + value[3]);
		s[2] = (value[2] - 2.0 * value[3] + value[4]);
		//weight
		double c_left[3], c_right[3];
		c_left[0] = 1.0 / 16.0;
		c_left[1] = 10.0 / 16.0;
		c_left[2] = 5.0 / 16.0;
		c_right[0] = 5.0 / 16.0;
		c_right[1] = 10.0 / 16.0;
		c_right[2] = 1.0 / 16.0;
		double eps = 1e-6;
		double IS[3];
		double beta_left[3], beta_right[3];
		for (int i = 0; i < 3; ++i)
		{
			IS[i] = f[i] * f[i] + s[i] * s[i] + eps;
			beta_left[i] = c_left[i] / (IS[i] * IS[i]);
			beta_right[i] = c_right[i] / (IS[i] * IS[i]);
		}
		double omega_left[3], omega_right[3];
		for (int i = 0; i < 3; ++i)
		{
			omega_left[i] = beta_left[i] / (beta_left[0] + beta_left[1] + beta_left[2]);
			omega_right[i] = beta_right[i] / (beta_right[0] + beta_right[1] + beta_right[2]);
		}
		double f_left, f_right;
		f_left = omega_left[0] * f[0] + omega_left[1] * f[1] + omega_left[2] * f[2];
		f_right = omega_right[0] * f[0] + omega_right[1] * f[1] + omega_right[2] * f[2];
		double s_left, s_right;
		s_left = omega_left[0] * s[0] + omega_left[1] * s[1] + omega_left[2] * s[2];
		s_right = omega_right[0] * s[0] + omega_right[1] * s[1] + omega_right[2] * s[2];
		value_left = value[2] + f_left / 2.0 + s_left / 8.0;
		value_right = value[2] - f_right / 2.0 + s_right / 8.0;
		//if (abs(ref_value) > SMALL_NUMBER)
		//{
		//	value_left = value_left * ref_value;
		//	value_right = value_right * ref_value;
		//}
	}
	void NSSolverStruct::FluxDifference()
	{
		auto para = GetPara();
		auto flux_difference_scheme = para->GetFluxDifferenceScheme();
		if (flux_difference_scheme == FluxDifferenceScheme::SecondOrder)
		{
			FluxDifference2nd();
		}
		else if (flux_difference_scheme == FluxDifferenceScheme::FourthOrder)
		{
			FluxDifference4th();
		}
		else if (flux_difference_scheme == FluxDifferenceScheme::SixthOrder)
		{
			FluxDifference6th();
		}
		else
		{
			Log::error("Flux difference scheme not found!");
			FluxDifference2nd();
		}
	}
	void NSSolverStruct::CalcConvectionResidual()
	{
		auto para = GetPara();
		auto inter_scheme = para->GetInterSchme();
		if (inter_scheme == InterpolationScheme::FirstOrder)
		{
			CalcConvectionRes_1st();
		}
		else if (inter_scheme == InterpolationScheme::MUSCL)
		{
			CalcConvectionRes_MUSCL();
		}
		else if (inter_scheme == InterpolationScheme::MUSCL_SV)
		{
			CalcConvectionRes_MUSCLSV();
		}
		else if (inter_scheme == InterpolationScheme::WCNS5)
		{
			CalcConvectionRes_WCNS5();
		}
		else if (inter_scheme == InterpolationScheme::Grad)
		{
			CalcConvectionRes_Grad();
		}
		else
		{
			Log::error("Interpolation scheme not found!");
			CalcConvectionRes_1st();
		}
	}
	void NSSolverStruct::CalcConvectionRes_1st()
	{
		InterMidNodePrim_1st();
		InterGhostMidNodePrim_WCNS5();
		FluxDifference2nd();
	}
	void NSSolverStruct::CalcConvectionRes_Grad()
	{
		InterMidNodePrim_Grad();
		FluxDifference();
	}
	void NSSolverStruct::CalcConvectionRes_MUSCL()
	{
		InterMidNodePrim_MUSCL();
		InterGhostMidNode_MUSCL();
		FluxDifference();
	}
	void NSSolverStruct::CalcConvectionRes_MUSCLSV()
	{
		InterMidNodePrim_MUSCLSV();
		CalcPrimGhostMidNodeMUSCL_SV();
		FluxDifference();
	}
	void NSSolverStruct::CalcConvectionRes_WCNS5()
	{
		int firstOrderSteps = GlobalData::GetInt("firstOrderSteps");
		int currentIter = GlobalData::GetInt("currentIter");
		if (currentIter < firstOrderSteps)
		{
			CalcConvectionRes_1st();
		}
		else
		{
			InterMidNodePrim_WCNS5();
			InterGhostMidNodePrim_WCNS5();
			FluxDifference();
		}
	}
	void NSSolverStruct::CalcViscousResidual()
	{
		// todo
	}
	void NSSolverStruct::CalcViscousFlux()
	{
		// todo
	}
	void NSSolverStruct::CalcViscousFluxGrad()
	{
		// todo
	}
	void NSSolverStruct::CalcSourceResidual()
	{
		// todo
	}
	void NSSolverStruct::BoundaryCondition()
	{
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		auto bound_map = grid->GetBoundMap();
		auto ni = grid->GetNi();
		auto nj = grid->GetNj();
		auto nk = grid->GetNk();
		for (auto& boundary : bound_map->GetBoundMap())
		{
			auto& bound_name = boundary.first;
			auto& bound = boundary.second;
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
			else
			{
				Log::error("Unsupport boundary: {}", bound_name);
				exit(0);
			}
		}
	}
	void NSSolverStruct::BCInflow(dynamic_array<BoundStruct>& bound)
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto data_manager = GetDataManager();
		int ghost_size = grid->GetGhostLevel();
		double prim_far[5] = { GetPara()->GetInflowDensity(), GetPara()->GetInflowVelocityX(),
							  GetPara()->GetInflowVelocityY(), GetPara()->GetInflowVelocityZ(),
							  GetPara()->GetInflowPressure() };
		//double prim_far[5] = { 8.0,8.25 * cos(30.0 / 180.0 * PI),-8.25 * sin(30.0 / 180.0 * PI),0.0,116.5 };
		//double prim_far[5] = { 6.4,3.125,0,0,18.5 };
		double cons_far[5];
		GetGas()->Prim2Cons(prim_far, cons_far);
#pragma omp parallel for 
		for (size_t iBound = 0; iBound < bound.size(); ++iBound)
		{
			index_type i_bound, j_bound, k_bound;
			index_type i_ghost, j_ghost, k_ghost;
			bound[iBound].GetIdx(i_bound, j_bound, k_bound);
			auto bound_direction = bound[iBound].GetDirectionSrc();
			for (size_t iGhost = 1; iGhost <= ghost_size; ++iGhost)
			{
				i_ghost = i_bound + iGhost * bound_direction[0];
				j_ghost = j_bound + iGhost * bound_direction[1];
				k_ghost = k_bound + iGhost * bound_direction[2];
				index_type idx_ghost = m_idx_proxy->GetIdx(i_ghost, j_ghost, k_ghost);
				//double y = node->GetCoord(i_ghost, j_ghost, k_ghost)[1];
				//prim_far[0] = 1 + 0.1 * y*y;
				//GetGas()->Prim2Cons(prim_far, cons_far);
				data_manager->SetPrim(idx_ghost, prim_far);
				data_manager->SetCons(idx_ghost, cons_far);
			}
		}
	}

	void NSSolverStruct::BCOutflow(dynamic_array<BoundStruct>& bound)
	{
		//BCDoubleMach(bound);
		//return;
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		size_t ghost_size = grid->GetGhostLevel();
		index_type eq_num = GetPara()->GetEqNum();
		dynamic_array<double> prim_vals(eq_num), cons_vals(eq_num);
#pragma omp parallel for private( prim_vals, cons_vals)
		for (size_t iBound = 0; iBound < bound.size(); ++iBound)
		{
			index_type i_bound, j_bound, k_bound;
			index_type i_ghost, j_ghost, k_ghost;
			bound[iBound].GetIdx(i_bound, j_bound, k_bound);
			index_type idx_bound = m_idx_proxy->GetIdx(i_bound, j_bound, k_bound);
			auto bound_direction = bound[iBound].GetDirectionSrc();
			// 预先获取边界的原始变量和守恒变量
			for (size_t idx_eq = 0; idx_eq < eq_num; ++idx_eq)
			{
				prim_vals[idx_eq] = data_manager->GetPrim(idx_eq, idx_bound);
			}
			GetGas()->Prim2Cons(prim_vals.data(), cons_vals.data());
			for (size_t iGhost = 1; iGhost <= ghost_size; ++iGhost)
			{
				i_ghost = i_bound + iGhost * bound_direction[0];
				j_ghost = j_bound + iGhost * bound_direction[1];
				k_ghost = k_bound + iGhost * bound_direction[2];
				index_type idx_ghost = m_idx_proxy->GetIdx(i_ghost, j_ghost, k_ghost);
				data_manager->SetPrim(idx_ghost, prim_vals.data());
				data_manager->SetCons(idx_ghost, cons_vals.data());
			}
		}
	}

	void NSSolverStruct::BCWall(dynamic_array<BoundStruct>& bound)
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto data_manager = GetDataManager();
		int ghost_size = grid->GetGhostLevel();
#pragma omp for
		for (size_t iBound = 0; iBound < bound.size(); ++iBound)
		{
			index_type i_bound, j_bound, k_bound;
			index_type i_ref, j_ref, k_ref;
			index_type i_ghost, j_ghost, k_ghost;
			auto norm_bnd = bound[iBound].GetNorm();
			bound[iBound].GetIdx(i_bound, j_bound, k_bound);
			index_type idx_bound = m_idx_proxy->GetIdx(i_bound, j_bound, k_bound);
			auto bound_coord = node->GetCoord(i_bound, j_bound, k_bound);
			for (size_t iGhost = 1; iGhost <= ghost_size; ++iGhost)
			{
				i_ghost = i_bound + iGhost * bound[iBound].GetDirectionSrc()[0];
				j_ghost = j_bound + iGhost * bound[iBound].GetDirectionSrc()[1];
				k_ghost = k_bound + iGhost * bound[iBound].GetDirectionSrc()[2];
				index_type idx_ghost = m_idx_proxy->GetIdx(i_ghost, j_ghost, k_ghost);
				i_ref = i_bound - iGhost * bound[iBound].GetDirectionSrc()[0];
				j_ref = j_bound - iGhost * bound[iBound].GetDirectionSrc()[1];
				k_ref = k_bound - iGhost * bound[iBound].GetDirectionSrc()[2];
				index_type idx_ref = m_idx_proxy->GetIdx(i_ref, j_ref, k_ref);
				double prim_ghost[5];
				for (int idx_eq = 0; idx_eq < data_manager->GetEqNum(); ++idx_eq)
				{
					prim_ghost[idx_eq] = data_manager->GetPrim(idx_eq, idx_ref);
				}
				double vel_ref[3] = { prim_ghost[1], prim_ghost[2], prim_ghost[3] };
				double vn_ref = vel_ref[0] * norm_bnd[0] + vel_ref[1] * norm_bnd[1] + vel_ref[2] * norm_bnd[2];
				prim_ghost[1] = vel_ref[0] - 2.0 * vn_ref * norm_bnd[0];
				prim_ghost[2] = vel_ref[1] - 2.0 * vn_ref * norm_bnd[1];
				prim_ghost[3] = vel_ref[2] - 2.0 * vn_ref * norm_bnd[2];
				double cons_ghost[5];
				GetGas()->Prim2Cons(prim_ghost, cons_ghost);
				data_manager->SetPrim(idx_ghost, prim_ghost);
				data_manager->SetCons(idx_ghost, cons_ghost);
			}
		}
	}

	void NSSolverStruct::BCFarfield(dynamic_array<BoundStruct>& bound)
	{
		auto data_manager = GetDataManager();
#pragma omp parallel for
		for (size_t iBound = 0; iBound < bound.size(); ++iBound)
		{
			index_type i_bound, j_bound, k_bound;
			index_type i_ghost, j_ghost, k_ghost;
			bound[iBound].GetIdx(i_bound, j_bound, k_bound);
			int idx_bound = m_idx_proxy->GetIdx(i_bound, j_bound, k_bound);
			i_ghost = i_bound + bound[iBound].GetDirectionSrc()[0];
			j_ghost = j_bound + bound[iBound].GetDirectionSrc()[1];
			k_ghost = k_bound + bound[iBound].GetDirectionSrc()[2];

			int idx_ref = m_idx_proxy->GetIdx(i_ghost, j_ghost, k_ghost);
			int idx_bnd = m_idx_proxy->GetIdx(i_bound, j_bound, k_bound);
			auto norm_bnd = bound[iBound].GetNorm();
			double prim_in[5] = { data_manager->GetPrim(0, idx_ref), data_manager->GetPrim(1, idx_ref),
								 data_manager->GetPrim(2, idx_ref), data_manager->GetPrim(3, idx_ref),
								 data_manager->GetPrim(4, idx_ref) };
			double prim_bnd[5] = { data_manager->GetPrim(0, idx_bnd), data_manager->GetPrim(1, idx_bnd),
								  data_manager->GetPrim(2, idx_bnd), data_manager->GetPrim(3, idx_bnd),
								  data_manager->GetPrim(4, idx_bnd) };
			auto para = GetPara();
			auto gas = GetGas();
			double gamma = gas->GetGamma();
			double prim_far[5] = { para->GetInflowDensity(), para->GetInflowVelocityX(), para->GetInflowVelocityY(),
								  para->GetInflowVelocityZ(), para->GetInflowPressure() };
			double vel_in[3] = { prim_in[1], prim_in[2], prim_in[3] };
			// 速度在边界法向上的投影
			double vn_in = vel_in[0] * norm_bnd[0] + vel_in[1] * norm_bnd[1] + vel_in[2] * norm_bnd[2];
			double vel_far[3] = { para->GetInflowVelocityX(), para->GetInflowVelocityY(), para->GetInflowVelocityZ() };
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
						prim_bnd[iVal] = prim_in[iVal];
					}
				}
				else // 超声速入口
				{
					for (int iVal = 0; iVal < 5; ++iVal)
					{
						prim_bnd[iVal] = prim_far[iVal];
					}
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
					prim_bnd[0] = pow(c_bound * c_bound / (entropy * gamma), 1.0 / gamma1);
					prim_bnd[1] = prim_far[1] + norm_bnd[0] * (vn_bound - vn_far);
					prim_bnd[2] = prim_far[2] + norm_bnd[1] * (vn_bound - vn_far);
					prim_bnd[3] = prim_far[3] + norm_bnd[2] * (vn_bound - vn_far);
					prim_bnd[4] = c_bound * c_bound * prim_bnd[0] / gamma;
				}
				else // 亚声速出口
				{
					double entropy = prim_in[4] / pow(prim_in[0], gamma);
					prim_bnd[0] = pow(c_bound * c_bound / (entropy * gamma), 1.0 / gamma1);
					prim_bnd[1] = prim_in[1] + norm_bnd[0] * (vn_bound - vn_in);
					prim_bnd[2] = prim_in[2] + norm_bnd[1] * (vn_bound - vn_in);
					prim_bnd[3] = prim_in[3] + norm_bnd[2] * (vn_bound - vn_in);
					prim_bnd[4] = c_bound * c_bound * prim_bnd[0] / gamma;
				}
			}
			double cons_bound[5];
			GetGas()->Prim2Cons(prim_bnd, cons_bound);
			auto grid = GetGrid();
			int ghost_size = grid->GetGhostLevel();
			for (size_t iGhost = 1; iGhost <= ghost_size; ++iGhost)
			{
				i_ghost = i_bound + iGhost * bound[iBound].GetDirectionSrc()[0];
				j_ghost = j_bound + iGhost * bound[iBound].GetDirectionSrc()[1];
				k_ghost = k_bound + iGhost * bound[iBound].GetDirectionSrc()[2];
				int idx_ghost = m_idx_proxy->GetIdx(i_ghost, j_ghost, k_ghost);
				data_manager->SetPrim(idx_ghost, prim_bnd);
				data_manager->SetCons(idx_ghost, cons_bound);
			}
		}
	}


	void NSSolverStruct::BCSymmetry(dynamic_array<BoundStruct>& bound)
	{
		BCWall(bound);
	}

	void NSSolverStruct::BCVortex(dynamic_array<BoundStruct>& bound)
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto para = GetPara();
		auto data_manager = GetDataManager();
		int ghost_size = grid->GetGhostLevel();
		double beta = 5.0;
		double r2;
		double prim[5];
		double gamma = 1.4;
#pragma omp for private( r2, prim)
		for (size_t iBound = 0; iBound < bound.size(); ++iBound)
		{
			index_type i_bound, j_bound, k_bound;
			index_type i_ghost, j_ghost, k_ghost;
			bound[iBound].GetIdx(i_bound, j_bound, k_bound);
			index_type idx_bound = m_idx_proxy->GetIdx(i_bound, j_bound, k_bound);
			for (size_t iGhost = 1; iGhost <= ghost_size; ++iGhost)
			{
				i_ghost = i_bound + iGhost * bound[iBound].GetDirectionSrc()[0];
				j_ghost = j_bound + iGhost * bound[iBound].GetDirectionSrc()[1];
				k_ghost = k_bound + iGhost * bound[iBound].GetDirectionSrc()[2];
				index_type idx_ghost = m_idx_proxy->GetIdx(i_ghost, j_ghost, k_ghost);
				auto coord = node->GetCoord(i_ghost, j_ghost, k_ghost);
				r2 = coord[0] * coord[0] + coord[1] * coord[1];
				prim[0] = pow(1.0 - (gamma - 1.0) * beta * beta * exp(1.0 - r2) / (8.0 * gamma * PI * PI), 1.0 / (gamma - 1.0));
				prim[4] = pow(prim[0], gamma);
				prim[1] = beta * exp(0.5 * (1.0 - r2)) / (2.0 * PI) * (-coord[1]);
				prim[2] = beta * exp(0.5 * (1.0 - r2)) / (2.0 * PI) * (coord[0]);
				prim[3] = 0.0;
				data_manager->SetPrim(idx_ghost, prim);
			}
		}
	}

	void NSSolverStruct::BCDoubleMach(dynamic_array<BoundStruct>& bound)
	{
		double current_time = GlobalData::GetDouble("currentTime");
		double shock_velocity = 10.0 / sin(PI / 3.0);
		double shock_x = 1.0 / 6.0 + 1.0 / tan(PI / 3.0) + shock_velocity * current_time;
		double prim_left[5] = { 8.0,8.25 * cos(30.0 / 180.0 * PI),-8.25 * sin(30.0 / 180.0 * PI),0.0,116.5 };
		double prim_right[5] = { 1.4,0.0,0.0,0.0,1.0 };
		double cons_left[5], cons_right[5];
		GetGas()->Prim2Cons(prim_left, cons_left);
		GetGas()->Prim2Cons(prim_right, cons_right);
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		int ghost_size = grid->GetGhostLevel();
#pragma omp for
		for (size_t iBound = 0; iBound < bound.size(); ++iBound)
		{
			index_type i_bound, j_bound, k_bound;
			bound[iBound].GetIdx(i_bound, j_bound, k_bound);
			index_type idx_bound = m_idx_proxy->GetIdx(i_bound, j_bound, k_bound);
			auto bound_direction = bound[iBound].GetDirectionSrc();
			for (size_t iGhost = 1; iGhost <= ghost_size; ++iGhost)
			{
				index_type i_ghost = i_bound + iGhost * bound_direction[0];
				index_type j_ghost = j_bound + iGhost * bound_direction[1];
				index_type k_ghost = k_bound + iGhost * bound_direction[2];
				index_type idx_ghost = m_idx_proxy->GetIdx(i_ghost, j_ghost, k_ghost);
				double x = grid->GetNode()->GetCoord(i_ghost, j_ghost, k_ghost)[0];
				if (x <= shock_x)
				{
					data_manager->SetPrim(idx_ghost, prim_left);
					data_manager->SetCons(idx_ghost, cons_left);
				}
				else
				{
					data_manager->SetPrim(idx_ghost, prim_right);
					data_manager->SetCons(idx_ghost, cons_right);
				}
			}
		}
	}

	void NSSolverStruct::CheckPrimtive()
	{

	}
	void NSSolverStruct::CheckResidual()
	{
	}
	void NSSolverStruct::FixPrimtive()
	{
	}
	void NSSolverStruct::BackupField(std::string& back_folder)
	{
	}

	GridStruct* NSSolverStruct::GetGrid()
	{
		return m_grid.get();
	}
	Metric* NSSolverStruct::GetMidMetricsI()
	{
		return m_metrics_mid_i;
	}

	Metric* NSSolverStruct::GetMidMetricsJ()
	{
		return m_metrics_mid_j;
	}

	Metric* NSSolverStruct::GetMidMetricsK()
	{
		return m_metrics_mid_k;
	}

	double NSSolverStruct::NodeDifferece2nd(double* mid_data)
	{
		return 0.5 * (mid_data[1] - mid_data[0]);
	}

	double NSSolverStruct::NodeDifferece4th(double* mid_data)
	{
		return 27.0 / 24.0 * (mid_data[2] - mid_data[1]) - 1.0 / 24.0 * (mid_data[3] - mid_data[0]);
	}

	double NSSolverStruct::NodeDifferece4thLeft(double* mid_data)
	{
		return 1.0 / 24.0 *
			(22.0 * mid_data[0] - 17.0 * mid_data[1] - 9.0 * mid_data[2] + 5.0 * mid_data[3] - 1.0 * mid_data[4]);
	}

	double NSSolverStruct::NodeDifferece4thRight(double* mid_data)
	{
		return -1.0 / 24.0 *
			(22.0 * mid_data[0] - 17.0 * mid_data[1] - 9.0 * mid_data[2] + 5.0 * mid_data[3] - 1.0 * mid_data[4]);
	}

	double NSSolverStruct::NodeDifferece6th(double* mid_data)
	{
		return 75.0 / 64.0 * (mid_data[3] - mid_data[2]) - 25.0 / 384.0 * (mid_data[4] - mid_data[1]) +
			3.0 / 640.0 * (mid_data[5] - mid_data[0]);
	}

	double NSSolverStruct::MidNodeInter2nd(double* node_data)
	{
		return 0.5 * (node_data[0] + node_data[1]);
	}

	double NSSolverStruct::MidNodeInter2ndLeft(double* node_data)
	{
		return 1.5 * node_data[0] - 0.5 * node_data[1];
	}

	double NSSolverStruct::MidNodeInter2ndRight(double* node_data)
	{
		return -0.5 * node_data[0] + 1.5 * node_data[1];
	}

	double NSSolverStruct::MidNodeInter4th(double* node_data)
	{
		return 1.0 / 16.0 * (-1.0 * node_data[0] + 9.0 * node_data[1] + 9.0 * node_data[2] - 1.0 * node_data[3]);
	}

	double NSSolverStruct::MidNodeInter4thLeft1(double* node_data)
	{
		return 1.0 / 16.0 * (-5.0 * node_data[0] + 21.0 * node_data[1] - 35.0 * node_data[2] + 35.0 * node_data[3]);
	}

	double NSSolverStruct::MidNodeInter4thLeft2(double* node_data)
	{
		return 1.0 / 16.0 * (1.0 * node_data[0] - 5.0 * node_data[1] + 15.0 * node_data[2] + 5.0 * node_data[3]);
	}

	double NSSolverStruct::MidNodeInter4thRight1(double* node_data)
	{
		return 1.0 / 16.0 * (35.0 * node_data[0] - 35.0 * node_data[1] + 21.0 * node_data[2] - 5.0 * node_data[3]);
	}

	double NSSolverStruct::MidNodeInter4thRight2(double* node_data)
	{
		return 1.0 / 16.0 * (5.0 * node_data[0] + 15.0 * node_data[1] - 5.0 * node_data[2] + 1.0 * node_data[3]);
	}

	double NSSolverStruct::MidNodeInter6th(double* node_data)
	{
		return 75.0 / 128.0 * (node_data[2] + node_data[3]) - 25.0 / 256.0 * (node_data[1] + node_data[4]) +
			3.0 / 256.0 * (node_data[0] + node_data[5]);
	}

} // namespace  zaran