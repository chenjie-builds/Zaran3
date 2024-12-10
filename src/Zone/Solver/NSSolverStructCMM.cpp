#include "NSSolverStructCMM.h"
#include "Log.h"
#include "NSSolverStruct.h"
#include "Metric.h"

namespace zaran
{
	NSSolverStructCMM::NSSolverStructCMM(Id index, string name, std::shared_ptr<FlowSolverParamStruct> para, std::shared_ptr < GridStruct> grid, std::shared_ptr < DataManagerNSStruct>data_manager)
		: NSSolverStruct(index, name, para, grid, data_manager)
	{
	}
	NSSolverStructCMM::~NSSolverStructCMM()
	{
	}

	void NSSolverStructCMM::CalcCoordTransCoef()
	{
		NSSolverStruct::CalcCoordTransCoef();
		CalcMidNodeMetrics();
	}
	void NSSolverStructCMM::CalcMidNodeMetrics()
	{
		auto para = GetPara();
		auto flux_diff_scheme = para->GetFluxDifferenceScheme();
		if (flux_diff_scheme == FluxDifferenceScheme::SecondOrder)
		{
			CalcMidNodeMetrics2nd();
		}
		else if (flux_diff_scheme == FluxDifferenceScheme::SixthOrder)
		{
			CalcMidNodeMetrics6th();
		}
		else
		{
			CalcMidNodeMetrics2nd();
		}
	}
	void NSSolverStructCMM::CalcMidNodeMetrics2nd()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto coef = GetNodeMetrics();
		Metric* coef_mid[3];
		coef_mid[0] = GetMidMetricsI();
		coef_mid[1] = GetMidMetricsJ();
		coef_mid[2] = GetMidMetricsK();
		auto idx_proxy = GetIdxProxy();
		int ni = grid->GetNi();
		int nj = grid->GetNj();
		int nk = grid->GetNk();
		int idx_temp[3];
		// 用于i,j,k方向寻找模板的系数
		int idx_direct[3][3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };
		// 内部半点：根据整点度量系数计算半点度量系数
		for (int k = 1; k < nk - 2; ++k)
		{
			for (int j = 1; j < nj - 2; ++j)
			{
				for (int i = 1; i < ni - 2; ++i)
				{
					int idx = idx_proxy->GetIdx(i, j, k);
					for (int iDim = 0; iDim < 3; iDim++)
					{
						auto coef_xi = coef_mid[iDim]->GetXi(idx);
						auto coef_eta = coef_mid[iDim]->GetEta(idx);
						auto coef_zeta = coef_mid[iDim]->GetZeta(idx);
						idx_temp[0] = idx_proxy->GetIdx(i, j, k);
						idx_temp[1] = idx_proxy->GetIdx(i + 1 * idx_direct[iDim][0], j + 1 * idx_direct[iDim][1],
							k + 1 * idx_direct[iDim][2]);
						for (int jDim = 0; jDim < 3; ++jDim)
						{
							coef_xi[jDim] = 0.5 * (coef->GetXi(idx_temp[0])[jDim] + coef->GetXi(idx_temp[1])[jDim]);
							coef_eta[jDim] = 0.5 * (coef->GetEta(idx_temp[0])[jDim] + coef->GetEta(idx_temp[1])[jDim]);
							coef_zeta[jDim] = 0.5 * (coef->GetZeta(idx_temp[0])[jDim] + coef->GetZeta(idx_temp[1])[jDim]);
						}
					}
				}
			}
		}
		// 边界半点：根据整点度量系数计算半点度量系数
		// i=1/2，i=(ni-1)-1/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				auto coef_xi_left = coef_mid[0]->GetXi(idx_proxy->GetIdx(0, j, k));           // xi[1/2]
				auto coef_eta_left = coef_mid[0]->GetEta(idx_proxy->GetIdx(0, j, k));         // eta[1/2]
				auto coef_zeta_left = coef_mid[0]->GetZeta(idx_proxy->GetIdx(0, j, k));       // zeta[1/2]
				auto coef_xi_right = coef_mid[0]->GetXi(idx_proxy->GetIdx(ni - 2, j, k));     // xi[(ni-1)-1/2]
				auto coef_eta_right = coef_mid[0]->GetEta(idx_proxy->GetIdx(ni - 2, j, k));   // eta[(ni-1)-1/2]
				auto coef_zeta_right = coef_mid[0]->GetZeta(idx_proxy->GetIdx(ni - 2, j, k)); // zeta[(ni-1)-1/2]
				for (int iDim = 0; iDim < 3; ++iDim)
				{
					coef_xi_left[iDim] = 0.5 * (3.0 * coef->GetXi(idx_proxy->GetIdx(1, j, k))[iDim] -
						coef->GetXi(idx_proxy->GetIdx(2, j, k))[iDim]);
					coef_eta_left[iDim] = 0.5 * (3.0 * coef->GetEta(idx_proxy->GetIdx(1, j, k))[iDim] -
						coef->GetEta(idx_proxy->GetIdx(2, j, k))[iDim]);
					coef_zeta_left[iDim] = 0.5 * (3.0 * coef->GetZeta(idx_proxy->GetIdx(1, j, k))[iDim] -
						coef->GetZeta(idx_proxy->GetIdx(2, j, k))[iDim]);
					coef_xi_right[iDim] = 0.5 * (3.0 * coef->GetXi(idx_proxy->GetIdx(ni - 2, j, k))[iDim] -
						coef->GetXi(idx_proxy->GetIdx(ni - 3, j, k))[iDim]);
					coef_eta_right[iDim] = 0.5 * (3.0 * coef->GetEta(idx_proxy->GetIdx(ni - 2, j, k))[iDim] -
						coef->GetEta(idx_proxy->GetIdx(ni - 3, j, k))[iDim]);
					coef_zeta_right[iDim] = 0.5 * (3.0 * coef->GetZeta(idx_proxy->GetIdx(ni - 2, j, k))[iDim] -
						coef->GetZeta(idx_proxy->GetIdx(ni - 3, j, k))[iDim]);
				}
			}
		}
		// j=1/2,j=(nj-1)-1/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				auto coef_xi_left = coef_mid[1]->GetXi(idx_proxy->GetIdx(i, 0, k));           // xi[1/2]
				auto coef_eta_left = coef_mid[1]->GetEta(idx_proxy->GetIdx(i, 0, k));         // eta[1/2]
				auto coef_zeta_left = coef_mid[1]->GetZeta(idx_proxy->GetIdx(i, 0, k));       // zeta[1/2]
				auto coef_xi_right = coef_mid[1]->GetXi(idx_proxy->GetIdx(i, nj - 2, k));     // xi[(nj-1)-1/2]
				auto coef_eta_right = coef_mid[1]->GetEta(idx_proxy->GetIdx(i, nj - 2, k));   // eta[(nj-1)-1/2]
				auto coef_zeta_right = coef_mid[1]->GetZeta(idx_proxy->GetIdx(i, nj - 2, k)); // zeta[(nj-1)-1/2]
				for (int iDim = 0; iDim < 3; ++iDim)
				{
					coef_xi_left[iDim] = 0.5 * (3.0 * coef->GetXi(idx_proxy->GetIdx(i, 1, k))[iDim] -
						coef->GetXi(idx_proxy->GetIdx(i, 2, k))[iDim]);
					coef_eta_left[iDim] = 0.5 * (3.0 * coef->GetEta(idx_proxy->GetIdx(i, 1, k))[iDim] -
						coef->GetEta(idx_proxy->GetIdx(i, 2, k))[iDim]);
					coef_zeta_left[iDim] = 0.5 * (3.0 * coef->GetZeta(idx_proxy->GetIdx(i, 1, k))[iDim] -
						coef->GetZeta(idx_proxy->GetIdx(i, 2, k))[iDim]);
					coef_xi_right[iDim] = 0.5 * (3.0 * coef->GetXi(idx_proxy->GetIdx(i, nj - 2, k))[iDim] -
						coef->GetXi(idx_proxy->GetIdx(i, nj - 3, k))[iDim]);
					coef_eta_right[iDim] = 0.5 * (3.0 * coef->GetEta(idx_proxy->GetIdx(i, nj - 2, k))[iDim] -
						coef->GetEta(idx_proxy->GetIdx(i, nj - 3, k))[iDim]);
					coef_zeta_right[iDim] = 0.5 * (3.0 * coef->GetZeta(idx_proxy->GetIdx(i, nj - 2, k))[iDim] -
						coef->GetZeta(idx_proxy->GetIdx(i, nj - 3, k))[iDim]);
				}
			}
		}
		// k=1/2,k=(nk-1)-1/2
		for (int j = 1; j < nj - 1; ++j)
		{
			for (int i = 1; i < ni - 1; ++i)
			{
				auto coef_xi_left = coef_mid[2]->GetXi(idx_proxy->GetIdx(i, j, 0));           // xi[1/2]
				auto coef_eta_left = coef_mid[2]->GetEta(idx_proxy->GetIdx(i, j, 0));         // eta[1/2]
				auto coef_zeta_left = coef_mid[2]->GetZeta(idx_proxy->GetIdx(i, j, 0));       // zeta[1/2]
				auto coef_xi_right = coef_mid[2]->GetXi(idx_proxy->GetIdx(i, j, nk - 2));     // xi[(nk-1)-1/2]
				auto coef_eta_right = coef_mid[2]->GetEta(idx_proxy->GetIdx(i, j, nk - 2));   // eta[(nk-1)-1/2]
				auto coef_zeta_right = coef_mid[2]->GetZeta(idx_proxy->GetIdx(i, j, nk - 2)); // zeta[(nk-1)-1/2]
				for (int iDim = 0; iDim < 3; ++iDim)
				{
					coef_xi_left[iDim] = 0.5 * (3.0 * coef->GetXi(idx_proxy->GetIdx(i, j, 1))[iDim] -
						coef->GetXi(idx_proxy->GetIdx(i, j, 2))[iDim]);
					coef_eta_left[iDim] = 0.5 * (3.0 * coef->GetEta(idx_proxy->GetIdx(i, j, 1))[iDim] -
						coef->GetEta(idx_proxy->GetIdx(i, j, 2))[iDim]);
					coef_zeta_left[iDim] = 0.5 * (3.0 * coef->GetZeta(idx_proxy->GetIdx(i, j, 1))[iDim] -
						coef->GetZeta(idx_proxy->GetIdx(i, j, 2))[iDim]);
					coef_xi_right[iDim] = 0.5 * (3.0 * coef->GetXi(idx_proxy->GetIdx(i, j, nk - 2))[iDim] -
						coef->GetXi(idx_proxy->GetIdx(i, j, nk - 3))[iDim]);
					coef_eta_right[iDim] = 0.5 * (3.0 * coef->GetEta(idx_proxy->GetIdx(i, j, nk - 2))[iDim] -
						coef->GetEta(idx_proxy->GetIdx(i, j, nk - 3))[iDim]);
					coef_zeta_right[iDim] = 0.5 * (3.0 * coef->GetZeta(idx_proxy->GetIdx(i, j, nk - 2))[iDim] -
						coef->GetZeta(idx_proxy->GetIdx(i, j, nk - 3))[iDim]);
				}
			}
		}
	}
	void NSSolverStructCMM::CalcMidNodeMetrics6th()
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
		int idx_temp[6];
		// 内部半点：根据整点度量系数计算半点度量系数
		// 第三步：根据整点逆变换度量系数计算半点逆变换度量系数
		double inter_temp[6];
		// i+1/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 3; i < ni - 3; ++i)
				{
					for (int iDim = 0; iDim < 3; iDim++)
					{
						for (int iTemp = 0; iTemp < 6; iTemp++)
						{
							inter_temp[iTemp] = coef->GetXi(Idx(i + iTemp - 2, j, k))[iDim];
						}
						coef_mid[0]->GetXi(Idx(i, j, k))[iDim] = MidNodeInter6th(inter_temp);
						for (int iTemp = 0; iTemp < 6; iTemp++)
						{
							inter_temp[iTemp] = coef->GetEta(Idx(i + iTemp - 2, j, k))[iDim];
						}
						coef_mid[0]->GetEta(Idx(i, j, k))[iDim] = MidNodeInter6th(inter_temp);
						for (int iTemp = 0; iTemp < 6; iTemp++)
						{
							inter_temp[iTemp] = coef->GetZeta(Idx(i + iTemp - 2, j, k))[iDim];
						}
						coef_mid[0]->GetZeta(Idx(i, j, k))[iDim] = MidNodeInter6th(inter_temp);
					}
				}
			}
		}
		// j+1/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 3; j < nj - 3; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					for (int iDim = 0; iDim < 3; iDim++)
					{
						for (int iTemp = 0; iTemp < 6; iTemp++)
						{
							inter_temp[iTemp] = coef->GetXi(Idx(i, j + iTemp - 2, k))[iDim];
						}
						coef_mid[1]->GetXi(Idx(i, j, k))[iDim] = MidNodeInter6th(inter_temp);
						for (int iTemp = 0; iTemp < 6; iTemp++)
						{
							inter_temp[iTemp] = coef->GetEta(Idx(i, j + iTemp - 2, k))[iDim];
						}
						coef_mid[1]->GetEta(Idx(i, j, k))[iDim] = MidNodeInter6th(inter_temp);
						for (int iTemp = 0; iTemp < 6; iTemp++)
						{
							inter_temp[iTemp] = coef->GetZeta(Idx(i, j + iTemp - 2, k))[iDim];
						}
						coef_mid[1]->GetZeta(Idx(i, j, k))[iDim] = MidNodeInter6th(inter_temp);
					}
				}
			}
		}
		// k+1/2
		for (int k = 3; k < nk - 3; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int i = 1; i < ni - 1; ++i)
				{
					for (int iDim = 0; iDim < 3; iDim++)
					{
						for (int iTemp = 0; iTemp < 6; iTemp++)
						{
							inter_temp[iTemp] = coef->GetXi(Idx(i, j, k + iTemp - 2))[iDim];
						}
						coef_mid[2]->GetXi(Idx(i, j, k))[iDim] = MidNodeInter6th(inter_temp);
						for (int iTemp = 0; iTemp < 6; iTemp++)
						{
							inter_temp[iTemp] = coef->GetEta(Idx(i, j, k + iTemp - 2))[iDim];
						}
						coef_mid[2]->GetEta(Idx(i, j, k))[iDim] = MidNodeInter6th(inter_temp);
						for (int iTemp = 0; iTemp < 6; iTemp++)
						{
							inter_temp[iTemp] = coef->GetZeta(Idx(i, j, k + iTemp - 2))[iDim];
						}
						coef_mid[2]->GetZeta(Idx(i, j, k))[iDim] = MidNodeInter6th(inter_temp);
					}
				}
			}
		}
		// 第四步：计算边界半点的逆变换度量系数(二阶偏置插值)
		// i=1/2,3/2,5/2;ni-3/2,ni-5/2,ni-7/2
		for (int k = 1; k < nk - 1; ++k)
		{
			for (int j = 1; j < nj - 1; ++j)
			{
				for (int iDim = 0; iDim < 3; ++iDim)
				{
					// i=1/2，3/2，5/2
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = coef->GetXi(Idx(1 + iTemp, j, k))[iDim];
					}
					coef_mid[0]->GetXi(Idx(0, j, k))[iDim] = MidNodeInter4thRight1(inter_temp);
					coef_mid[0]->GetXi(Idx(1, j, k))[iDim] = MidNodeInter4thRight2(inter_temp);
					coef_mid[0]->GetXi(Idx(2, j, k))[iDim] = MidNodeInter4th(inter_temp);
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = coef->GetEta(Idx(1 + iTemp, j, k))[iDim];
					}
					coef_mid[0]->GetEta(Idx(0, j, k))[iDim] = MidNodeInter4thRight1(inter_temp);
					coef_mid[0]->GetEta(Idx(1, j, k))[iDim] = MidNodeInter4thRight2(inter_temp);
					coef_mid[0]->GetEta(Idx(2, j, k))[iDim] = MidNodeInter4th(inter_temp);
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = coef->GetZeta(Idx(1 + iTemp, j, k))[iDim];
					}
					coef_mid[0]->GetZeta(Idx(0, j, k))[iDim] = MidNodeInter4thRight1(inter_temp);
					coef_mid[0]->GetZeta(Idx(1, j, k))[iDim] = MidNodeInter4thRight2(inter_temp);
					coef_mid[0]->GetZeta(Idx(2, j, k))[iDim] = MidNodeInter4th(inter_temp);
					// i=ni-3/2,ni-5/2,ni-7/2
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = coef->GetXi(Idx(ni - 5 + iTemp, j, k))[iDim];
					}
					coef_mid[0]->GetXi(Idx(ni - 2, j, k))[iDim] = MidNodeInter4thLeft1(inter_temp);
					coef_mid[0]->GetXi(Idx(ni - 3, j, k))[iDim] = MidNodeInter4thLeft2(inter_temp);
					coef_mid[0]->GetXi(Idx(ni - 4, j, k))[iDim] = MidNodeInter4th(inter_temp);
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = coef->GetEta(Idx(ni - 5 + iTemp, j, k))[iDim];
					}
					coef_mid[0]->GetEta(Idx(ni - 2, j, k))[iDim] = MidNodeInter4thLeft1(inter_temp);
					coef_mid[0]->GetEta(Idx(ni - 3, j, k))[iDim] = MidNodeInter4thLeft2(inter_temp);
					coef_mid[0]->GetEta(Idx(ni - 4, j, k))[iDim] = MidNodeInter4th(inter_temp);
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = coef->GetZeta(Idx(ni - 5 + iTemp, j, k))[iDim];
					}
					coef_mid[0]->GetZeta(Idx(ni - 2, j, k))[iDim] = MidNodeInter4thLeft1(inter_temp);
					coef_mid[0]->GetZeta(Idx(ni - 3, j, k))[iDim] = MidNodeInter4thLeft2(inter_temp);
					coef_mid[0]->GetZeta(Idx(ni - 4, j, k))[iDim] = MidNodeInter4th(inter_temp);
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
					// j=1/2，3/2，5/2
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = coef->GetXi(Idx(i, 1 + iTemp, k))[iDim];
					}
					coef_mid[1]->GetXi(Idx(i, 0, k))[iDim] = MidNodeInter4thRight1(inter_temp);
					coef_mid[1]->GetXi(Idx(i, 1, k))[iDim] = MidNodeInter4thRight2(inter_temp);
					coef_mid[1]->GetXi(Idx(i, 2, k))[iDim] = MidNodeInter4th(inter_temp);
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = coef->GetEta(Idx(i, 1 + iTemp, k))[iDim];
					}
					coef_mid[1]->GetEta(Idx(i, 0, k))[iDim] = MidNodeInter4thRight1(inter_temp);
					coef_mid[1]->GetEta(Idx(i, 1, k))[iDim] = MidNodeInter4thRight2(inter_temp);
					coef_mid[1]->GetEta(Idx(i, 2, k))[iDim] = MidNodeInter4th(inter_temp);
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = coef->GetZeta(Idx(i, 1 + iTemp, k))[iDim];
					}
					coef_mid[1]->GetZeta(Idx(i, 0, k))[iDim] = MidNodeInter4thRight1(inter_temp);
					coef_mid[1]->GetZeta(Idx(i, 1, k))[iDim] = MidNodeInter4thRight2(inter_temp);
					coef_mid[1]->GetZeta(Idx(i, 2, k))[iDim] = MidNodeInter4th(inter_temp);
					// j=nj-3/2,nj-5/2,nj-7/2
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = coef->GetXi(Idx(i, nj - 2 - iTemp, k))[iDim];
					}
					coef_mid[1]->GetXi(Idx(i, nj - 2, k))[iDim] = MidNodeInter4thLeft1(inter_temp);
					coef_mid[1]->GetXi(Idx(i, nj - 3, k))[iDim] = MidNodeInter4thLeft2(inter_temp);
					coef_mid[1]->GetXi(Idx(i, nj - 4, k))[iDim] = MidNodeInter4th(inter_temp);
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = coef->GetEta(Idx(i, nj - 5 + iTemp, k))[iDim];
					}
					coef_mid[1]->GetEta(Idx(i, nj - 2, k))[iDim] = MidNodeInter4thLeft1(inter_temp);
					coef_mid[1]->GetEta(Idx(i, nj - 3, k))[iDim] = MidNodeInter4thLeft2(inter_temp);
					coef_mid[1]->GetEta(Idx(i, nj - 4, k))[iDim] = MidNodeInter4th(inter_temp);
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = coef->GetZeta(Idx(i, nj - 5 + iTemp, k))[iDim];
					}
					coef_mid[1]->GetZeta(Idx(i, nj - 2, k))[iDim] = MidNodeInter4thLeft1(inter_temp);
					coef_mid[1]->GetZeta(Idx(i, nj - 3, k))[iDim] = MidNodeInter4thLeft2(inter_temp);
					coef_mid[1]->GetZeta(Idx(i, nj - 4, k))[iDim] = MidNodeInter4th(inter_temp);
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
					// k=1/2，3/2，5/2
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = coef->GetXi(Idx(i, j, 1 + iTemp))[iDim];
					}
					coef_mid[2]->GetXi(Idx(i, j, 0))[iDim] = MidNodeInter4thRight1(inter_temp);
					coef_mid[2]->GetXi(Idx(i, j, 1))[iDim] = MidNodeInter4thRight2(inter_temp);
					coef_mid[2]->GetXi(Idx(i, j, 2))[iDim] = MidNodeInter4th(inter_temp);
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = coef->GetEta(Idx(i, j, 1 + iTemp))[iDim];
					}
					coef_mid[2]->GetEta(Idx(i, j, 0))[iDim] = MidNodeInter4thRight1(inter_temp);
					coef_mid[2]->GetEta(Idx(i, j, 1))[iDim] = MidNodeInter4thRight2(inter_temp);
					coef_mid[2]->GetEta(Idx(i, j, 2))[iDim] = MidNodeInter4th(inter_temp);
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = coef->GetZeta(Idx(i, j, 1 + iTemp))[iDim];
					}
					coef_mid[2]->GetZeta(Idx(i, j, 0))[iDim] = MidNodeInter4thRight1(inter_temp);
					coef_mid[2]->GetZeta(Idx(i, j, 1))[iDim] = MidNodeInter4thRight2(inter_temp);
					coef_mid[2]->GetZeta(Idx(i, j, 2))[iDim] = MidNodeInter4th(inter_temp);
					// k=nk-3/2,nk-5/2,nk-7/2
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = coef->GetXi(Idx(i, j, nk - 5 + iTemp))[iDim];
					}
					coef_mid[2]->GetXi(Idx(i, j, nk - 2))[iDim] = MidNodeInter4thLeft1(inter_temp);
					coef_mid[2]->GetXi(Idx(i, j, nk - 3))[iDim] = MidNodeInter4thLeft2(inter_temp);
					coef_mid[2]->GetXi(Idx(i, j, nk - 4))[iDim] = MidNodeInter4th(inter_temp);
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = coef->GetEta(Idx(i, j, nk - 5 + iTemp))[iDim];
					}
					coef_mid[2]->GetEta(Idx(i, j, nk - 2))[iDim] = MidNodeInter4thLeft1(inter_temp);
					coef_mid[2]->GetEta(Idx(i, j, nk - 3))[iDim] = MidNodeInter4thLeft2(inter_temp);
					coef_mid[2]->GetEta(Idx(i, j, nk - 4))[iDim] = MidNodeInter4th(inter_temp);
					for (int iTemp = 0; iTemp < 4; iTemp++)
					{
						inter_temp[iTemp] = coef->GetZeta(Idx(i, j, nk - 5 + iTemp))[iDim];
					}
					coef_mid[2]->GetZeta(Idx(i, j, nk - 2))[iDim] = MidNodeInter4thLeft1(inter_temp);
					coef_mid[2]->GetZeta(Idx(i, j, nk - 3))[iDim] = MidNodeInter4thLeft2(inter_temp);
					coef_mid[2]->GetZeta(Idx(i, j, nk - 4))[iDim] = MidNodeInter4th(inter_temp);
				}
			}
		}
	}

	void NSSolverStructCMM::CalcConvectionRes_1st()
	{
		InterMidNodePrim_1st();
		InterGhostMidNode_MUSCL();
		CalcConvectionFluxMidNode();
		FluxDifference2nd();
	}

	void NSSolverStructCMM::CalcConvectionRes_MUSCL()
	{
		InterMidNodePrim_MUSCL();
		InterGhostMidNode_MUSCL();
		CalcConvectionFluxMidNode();
		FluxDifference();
	}

	void NSSolverStructCMM::CalcConvectionRes_WCNS5()
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
			CalcConvectionFluxMidNode();
			FluxDifference();
		}
	}

	void NSSolverStructCMM::CalcConvectionFluxMidNode()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto data_manager = GetDataManager();
		auto node_metrics = GetNodeMetrics();
		auto idx_proxy = GetIdxProxy();
		auto gas = GetGas();
		auto equ_num = GetPara()->GetEqNum();
		int ni, nj, nk;
		ni = grid->GetNi();
		nj = grid->GetNj();
		nk = grid->GetNk();
		RiemannSolverPara riemann_para;
		for (int k = 0; k < nk - 1; ++k)
		{
#pragma omp parallel for private(riemann_para)
			for (int j = 0; j < nj - 1; ++j)
			{
				for (int i = 0; i < ni - 1; ++i)
				{
					riemann_para.gamma_left = riemann_para.gamma_right = gas->GetGamma();
					int idx = idx_proxy->GetIdx(i, j, k);
					double jacobi = node_metrics->GetJacobian(idx);
					// i direction
					riemann_para.norm[0] = GetMidMetricsI()->GetXi(idx)[0];
					riemann_para.norm[1] = GetMidMetricsI()->GetXi(idx)[1];
					riemann_para.norm[2] = GetMidMetricsI()->GetXi(idx)[2];
					riemann_para.nt = GetMidMetricsI()->GetXi(idx)[3];
					for (int iEqu = 0; iEqu < equ_num; ++iEqu)
					{
						riemann_para.prim_left[iEqu] = data_manager->GetMidNodePrimLeft(iEqu, 0, idx);
						riemann_para.prim_right[iEqu] = data_manager->GetMidNodePrimRight(iEqu, 0, idx);
					}
					//if (abs(riemann_para.prim_left[0] -riemann_para.prim_right[0])>1e-10)
					//{
					//	double y = 0.5 * (node->GetCoord(i, j, k)[1] + node->GetCoord(i + 1, j, k)[1]);
					//	double density_exact = 1 + 0.1 * y;
					//	Log::error("i direction: left and right density are not equal");
					//	Log::error("(i,j,k)=({},{},{}) left and right density are not equal", i, j, k);
					//	Log::error("density left= {}, right= {}, exact= {}", riemann_para.prim_left[0], riemann_para.prim_right[0], density_exact);
					//}
					m_riemann_solver->Solver(riemann_para);
					for (int iEqu = 0; iEqu < equ_num; ++iEqu)
					{
						data_manager->SetMidNodeFlux(iEqu, 0, idx, riemann_para.flux[iEqu]);
					}
					// j direction
					riemann_para.norm[0] = GetMidMetricsJ()->GetEta(idx)[0];
					riemann_para.norm[1] = GetMidMetricsJ()->GetEta(idx)[1];
					riemann_para.norm[2] = GetMidMetricsJ()->GetEta(idx)[2];
					riemann_para.nt = GetMidMetricsJ()->GetEta(idx)[3];
					for (int iEqu = 0; iEqu < equ_num; ++iEqu)
					{
						riemann_para.prim_left[iEqu] = data_manager->GetMidNodePrimLeft(iEqu, 1, idx);
						riemann_para.prim_right[iEqu] = data_manager->GetMidNodePrimRight(iEqu, 1, idx);
					}
					//if (abs(riemann_para.prim_left[0] - riemann_para.prim_right[0]) > 1e-10)
					//{
					//	double y = 0.5 * (node->GetCoord(i, j, k)[1] + node->GetCoord(i , j+1, k)[1]);
					//	double density_exact = 1 + 0.1 * y;
					//	Log::error("j direction: left and right density are not equal");
					//	Log::error("(i,j,k)=({},{},{}) left and right density are not equal", i, j, k);
					//	Log::error("density left= {}, right= {}, exact= {}", riemann_para.prim_left[0], riemann_para.prim_right[0], density_exact);
					//}
					m_riemann_solver->Solver(riemann_para);
					for (int iEqu = 0; iEqu < equ_num; ++iEqu)
					{
						data_manager->SetMidNodeFlux(iEqu, 1, idx, riemann_para.flux[iEqu]);
					}

					// k direction
					if (grid->GetDim() == 3)
					{
						riemann_para.norm[0] = GetMidMetricsK()->GetZeta(idx)[0];
						riemann_para.norm[1] = GetMidMetricsK()->GetZeta(idx)[1];
						riemann_para.norm[2] = GetMidMetricsK()->GetZeta(idx)[2];
						riemann_para.nt = GetMidMetricsK()->GetZeta(idx)[3];
						for (int iEqu = 0; iEqu < equ_num; ++iEqu)
						{
							riemann_para.prim_left[iEqu] = data_manager->GetMidNodePrimLeft(iEqu, 2, idx);
							riemann_para.prim_right[iEqu] = data_manager->GetMidNodePrimRight(iEqu, 2, idx);
						}
						m_riemann_solver->Solver(riemann_para);
						for (int iEqu = 0; iEqu < equ_num; ++iEqu)
						{
							data_manager->SetMidNodeFlux(iEqu, 2, idx, riemann_para.flux[iEqu]);
						}
					}
				}
			}
		}
	}
	void NSSolverStructCMM::FluxDifference2nd()
	{
		auto grid = GetGrid();
		auto para = GetPara();
		auto data_manager = GetDataManager();
		auto node_metrics = GetNodeMetrics();
		auto idx_proxy = GetIdxProxy();
		auto gas = GetGas();
		auto equ_num = para->GetEqNum();
		int is, ie, js, je, ks, ke;
		grid->GetRange(is, ie, js, je, ks, ke);
		// 差分模板的编号
		int idx_temp[3];
		// 残差的临时变量
		double res_tmp[5];
		double direction[3][3] = { {1,0,0},{0,1,0},{0,0,1} };
		for (int k = ks; k <= ke; ++k)
		{
#pragma omp parallel for private( idx_temp, res_tmp)
			for (int j = js; j <= je; ++j)
			{
				for (int i = is; i <= ie; ++i)
				{
					int idx = idx_proxy->GetIdx(i, j, k);
					for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
					{
						res_tmp[idx_eq] = data_manager->GetResidual(idx_eq, idx);
					}
					for (int dim = 0; dim < grid->GetDim(); ++dim)
					{
						for (int iTemp = 0; iTemp < 3; iTemp++)
						{
							idx_temp[iTemp] = idx_proxy->GetIdx(
								i + (iTemp - 1) * direction[dim][0],
								j + (iTemp - 1) * direction[dim][1],
								k + (iTemp - 1) * direction[dim][2]);
						}
						for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
						{
							res_tmp[idx_eq] -=
								data_manager->GetMidNodeFlux(idx_eq, dim, idx_temp[1]) - data_manager->GetMidNodeFlux(idx_eq, dim, idx_temp[0]);
						}
					}
					data_manager->SetResidual(idx, res_tmp);
				}
			}
		}
	}
	void NSSolverStructCMM::FluxDifference6th()
	{
		auto grid = GetGrid();
		auto para = GetPara();
		auto data_manager = GetDataManager();
		auto node_metrics = GetNodeMetrics();
		auto idx_proxy = GetIdxProxy();
		auto gas = GetGas();
		auto equ_num = para->GetEqNum();
		int is, ie, js, je, ks, ke;
		grid->GetRange(is, ie, js, je, ks, ke);
		// 差分模板的编号
		int idx_temp[6];
		double res_tmp[5];
		double coef1 = 75.0 / 64.0;
		double coef2 = -25.0 / 384.0;
		double coef3 = 3.0 / 640.0;
		double direction[3][3] = { {1,0,0},{0,1,0},{0,0,1} };
		for (int k = ks; k <= ke; ++k)
		{
#pragma omp parallel for private(  idx_temp, res_tmp)
			for (int j = js; j <= je; ++j)
			{
				for (int i = is; i <= ie; ++i)
				{
					int idx = idx_proxy->GetIdx(i, j, k);
					for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
					{
						res_tmp[idx_eq] = data_manager->GetResidual(idx_eq, idx);
					}
					for (int dim = 0; dim < grid->GetDim(); ++dim)
					{
						for (int iTemp = 0; iTemp < 6; iTemp++)
						{
							idx_temp[iTemp] = idx_proxy->GetIdx(
								i + (iTemp - 3) * direction[dim][0],
								j + (iTemp - 3) * direction[dim][1],
								k + (iTemp - 3) * direction[dim][2]);
						}
						for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
						{
							res_tmp[idx_eq] -=
								coef1 * (data_manager->GetMidNodeFlux(idx_eq, dim, idx_temp[3]) - data_manager->GetMidNodeFlux(idx_eq, dim, idx_temp[2])) +
								coef2 * (data_manager->GetMidNodeFlux(idx_eq, dim, idx_temp[4]) - data_manager->GetMidNodeFlux(idx_eq, dim, idx_temp[1])) +
								coef3 * (data_manager->GetMidNodeFlux(idx_eq, dim, idx_temp[5]) - data_manager->GetMidNodeFlux(idx_eq, dim, idx_temp[0]));
						}
					}
					data_manager->SetResidual(idx, res_tmp);
				}
			}
		}
	}

} // namespace zaran