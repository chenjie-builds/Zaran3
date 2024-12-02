#include "NSSolverStructDEER.h"
namespace zaran
{
	NSSolverStructDEER::NSSolverStructDEER(int index, string name, FlowSolverPara* para, GridStruct* grid,
		DataManagerNSStruct* data_manager)
		: NSSolverStruct(index, name, para, grid, data_manager)
	{
	}
	NSSolverStructDEER::~NSSolverStructDEER()
	{
	}
	void NSSolverStructDEER::CalcFluxDifference2ndOrder()
	{
		auto para = GetPara();
		if (para->GetInterSchme() == InterpolationScheme::Grad)
		{
			FlueDifference2nd_MoveMidNode();
			return;
		}
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		auto node_metrics = GetNodeMetrics();
		auto idx_proxy = GetIdxProxy();
		auto gas = GetGas();
		auto equ_num = para->GetEqNum();
		int is, ie, js, je, ks, ke;
		grid->GetRange(is, ie, js, je, ks, ke);
		RiemannSolverPara riemann_para[2];
		double direction[3][3] = { {1,0,0},{0,1,0},{0,0,1} };
		double res_tmp[5];
		for (int k = ks; k <= ke; ++k)
		{
#pragma omp parallel for private( res_tmp, riemann_para)
			for (int j = js; j <= je; ++j)
			{
				for (int i = is; i <= ie; ++i)
				{
					int idx = idx_proxy->GetIdx(i, j, k);
					for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
					{
						res_tmp[idx_eq] = data_manager->GetResidual(idx_eq, idx);
					}
					for (int dim = 0; dim < grid->GetDim(); dim++)
					{
						for (int iTemp = 0; iTemp < 2; iTemp++)
						{
							riemann_para[iTemp].gamma_left = riemann_para[iTemp].gamma_right = gas->GetGamma();
							riemann_para[iTemp].norm(0) = node_metrics->GetMetrics(dim, idx)[0];
							riemann_para[iTemp].norm(1) = node_metrics->GetMetrics(dim, idx)[1];
							riemann_para[iTemp].norm(2) = node_metrics->GetMetrics(dim, idx)[2];
							riemann_para[iTemp].nt = node_metrics->GetMetrics(dim, idx)[3];
							int idx_temp = idx_proxy->GetIdx(
								i + (iTemp - 1) * direction[dim][0],
								j + (iTemp - 1) * direction[dim][1],
								k + (iTemp - 1) * direction[dim][2]);
							for (int iVal = 0; iVal < equ_num; ++iVal)
							{
								riemann_para[iTemp].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, dim, idx_temp);
								riemann_para[iTemp].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, dim, idx_temp);
							}
							m_riemann_solver->Solver(riemann_para[iTemp]);
						}
						for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
						{
							res_tmp[idx_eq] -=
								riemann_para[1].flux[idx_eq] - riemann_para[0].flux[idx_eq];
						}
					}
					data_manager->SetResidual(idx, res_tmp);

				}
			}
		}
	}

	void NSSolverStructDEER::FlueDifference2nd_MoveMidNode()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto para = GetPara();
		auto data_manager = GetDataManager();
		auto node_metrics = GetNodeMetrics();
		auto idx_proxy = GetIdxProxy();
		auto gas = GetGas();
		auto equ_num = para->GetEqNum();
		int is, ie, js, je, ks, ke;
		grid->GetRange(is, ie, js, je, ks, ke);
		RiemannSolverPara riemann_para[6];
		for (int i = 0; i < 6; ++i)
		{
			riemann_para[i].gamma_left = riemann_para[i].gamma_right = gas->GetGamma();
		}
		// 当前节点的编号
		int idx;
		// 差分模板的编号
		int idx_temp[5];
		// 差分模板的值
		double value[5];
		double res_tmp[5];
		double mid_coord_left[3], mid_coord_right[3]; // i-1/2,i+1/2
		double mid_coord[3];                          // i-1/2,i+1/2的中点
		double move_vector[3];                        // mid_coord到i的向量用于移动mid_coord_left和mid_coord_right
#pragma omp parallel for private(idx, idx_temp, value, res_tmp, mid_coord_left, mid_coord_right, mid_coord,            \
                                     move_vector, riemann_para)
		for (int k = ks; k <= ke; ++k)
		{
			for (int j = js; j <= je; ++j)
			{
				for (int i = is; i <= ie; ++i)
				{
					for (int iVal = 0; iVal < equ_num; ++iVal)
					{
						res_tmp[iVal] = data_manager->GetResidual(iVal, idx_proxy->GetIdx(i, j, k));
					}
					// i direction
					auto left_coord = node->GetCoord(i - 1, j, k);
					auto coord = node->GetCoord(i, j, k);
					auto right_coord = node->GetCoord(i + 1, j, k);
					for (int iDim = 0; iDim < 3; ++iDim)
					{
						mid_coord_left[iDim] = (left_coord[iDim] + coord[iDim]) / 2;
						mid_coord_right[iDim] = (right_coord[iDim] + coord[iDim]) / 2;
						mid_coord[iDim] = (mid_coord_left[iDim] + mid_coord_right[iDim]) / 2;
						move_vector[iDim] = coord[iDim] - mid_coord[iDim];
						mid_coord_left[iDim] += move_vector[iDim];
						mid_coord_right[iDim] += move_vector[iDim];
					}
					idx = idx_proxy->GetIdx(i, j, k);
					riemann_para[0].norm(0) = node_metrics->GetXi(idx)[0];
					riemann_para[0].norm(1) = node_metrics->GetXi(idx)[1];
					riemann_para[0].norm(2) = node_metrics->GetXi(idx)[2];
					riemann_para[0].nt = node_metrics->GetXi(idx)[3];
					riemann_para[1].norm(0) = node_metrics->GetXi(idx)[0];
					riemann_para[1].norm(1) = node_metrics->GetXi(idx)[1];
					riemann_para[1].norm(2) = node_metrics->GetXi(idx)[2];
					riemann_para[1].nt = node_metrics->GetXi(idx)[3];
					MidNodeGrad(idx_proxy->GetIdx(i, j, k), idx_proxy->GetIdx(i + 1, j, k), coord, mid_coord_right,
						right_coord, riemann_para[0].prim_left.data(), riemann_para[0].prim_right.data());
					MidNodeGrad(idx_proxy->GetIdx(i - 1, j, k), idx_proxy->GetIdx(i, j, k), left_coord, mid_coord_left,
						coord, riemann_para[1].prim_left.data(), riemann_para[1].prim_right.data());
					m_riemann_solver->Solver(riemann_para[0]);
					m_riemann_solver->Solver(riemann_para[1]);
					for (int iVal = 0; iVal < equ_num; ++iVal)
					{
						res_tmp[iVal] -= (riemann_para[0].flux[iVal] - riemann_para[1].flux[iVal]);
					}
					// j direction
					left_coord = node->GetCoord(i, j - 1, k);
					coord = node->GetCoord(i, j, k);
					right_coord = node->GetCoord(i, j + 1, k);
					for (int iDim = 0; iDim < 3; ++iDim)
					{
						mid_coord_left[iDim] = (left_coord[iDim] + coord[iDim]) / 2;
						mid_coord_right[iDim] = (right_coord[iDim] + coord[iDim]) / 2;
						mid_coord[iDim] = (mid_coord_left[iDim] + mid_coord_right[iDim]) / 2;
						move_vector[iDim] = coord[iDim] - mid_coord[iDim];
						mid_coord_left[iDim] += move_vector[iDim];
						mid_coord_right[iDim] += move_vector[iDim];
					}
					riemann_para[2].norm(0) = node_metrics->GetEta(idx)[0];
					riemann_para[2].norm(1) = node_metrics->GetEta(idx)[1];
					riemann_para[2].norm(2) = node_metrics->GetEta(idx)[2];
					riemann_para[2].nt = node_metrics->GetEta(idx)[3];
					riemann_para[3].norm(0) = node_metrics->GetEta(idx)[0];
					riemann_para[3].norm(1) = node_metrics->GetEta(idx)[1];
					riemann_para[3].norm(2) = node_metrics->GetEta(idx)[2];
					riemann_para[3].nt = node_metrics->GetEta(idx)[3];
					MidNodeGrad(idx_proxy->GetIdx(i, j, k), idx_proxy->GetIdx(i, j + 1, k), coord, mid_coord_right,
						right_coord, riemann_para[2].prim_left.data(), riemann_para[2].prim_right.data());
					MidNodeGrad(idx_proxy->GetIdx(i, j - 1, k), idx_proxy->GetIdx(i, j, k), left_coord, mid_coord_left,
						coord, riemann_para[3].prim_left.data(), riemann_para[3].prim_right.data());
					m_riemann_solver->Solver(riemann_para[2]);
					m_riemann_solver->Solver(riemann_para[3]);
					for (int iVal = 0; iVal < equ_num; ++iVal)
					{
						res_tmp[iVal] -= (riemann_para[2].flux[iVal] - riemann_para[3].flux[iVal]);
					}
					// k direction
					if (grid->GetDim() == 3)
					{
						left_coord = node->GetCoord(i, j, k - 1);
						coord = node->GetCoord(i, j, k);
						right_coord = node->GetCoord(i, j, k + 1);
						for (int iDim = 0; iDim < 3; ++iDim)
						{
							mid_coord_left[iDim] = (left_coord[iDim] + coord[iDim]) / 2;
							mid_coord_right[iDim] = (right_coord[iDim] + coord[iDim]) / 2;
							mid_coord[iDim] = (mid_coord_left[iDim] + mid_coord_right[iDim]) / 2;
							move_vector[iDim] = coord[iDim] - mid_coord[iDim];
							mid_coord_left[iDim] += move_vector[iDim];
							mid_coord_right[iDim] += move_vector[iDim];
						}
						riemann_para[4].norm(0) = node_metrics->GetZeta(idx)[0];
						riemann_para[4].norm(1) = node_metrics->GetZeta(idx)[1];
						riemann_para[4].norm(2) = node_metrics->GetZeta(idx)[2];
						riemann_para[4].nt = node_metrics->GetZeta(idx)[3];
						riemann_para[5].norm(0) = node_metrics->GetZeta(idx)[0];
						riemann_para[5].norm(1) = node_metrics->GetZeta(idx)[1];
						riemann_para[5].norm(2) = node_metrics->GetZeta(idx)[2];
						riemann_para[5].nt = node_metrics->GetZeta(idx)[3];
						MidNodeGrad(idx_proxy->GetIdx(i, j, k), idx_proxy->GetIdx(i, j, k + 1), coord, mid_coord_right,
							right_coord, riemann_para[4].prim_left.data(), riemann_para[4].prim_right.data());
						MidNodeGrad(idx_proxy->GetIdx(i, j, k - 1), idx_proxy->GetIdx(i, j, k), left_coord, mid_coord_left,
							coord, riemann_para[5].prim_left.data(), riemann_para[5].prim_right.data());
						m_riemann_solver->Solver(riemann_para[4]);
						m_riemann_solver->Solver(riemann_para[5]);
						for (int iVal = 0; iVal < equ_num; ++iVal)
						{
							res_tmp[iVal] -= (riemann_para[4].flux[iVal] - riemann_para[5].flux[iVal]);
						}
					}
					for (int iVar = 0; iVar < equ_num; ++iVar)
					{
						data_manager->SetResidual(iVar, idx, res_tmp[iVar]);
					}
				}
			}
		}
	}
	void NSSolverStructDEER::CalcFluxDifference4thOrder()
	{
		auto grid = GetGrid();
		auto para = GetPara();
		auto node = grid->GetNode();
		auto data_manager = GetDataManager();
		auto idx_proxy = GetIdxProxy();
		auto node_metrics = GetNodeMetrics();
		auto gas = GetGas();
		auto equ_num = para->GetEqNum();
		int is, ie, js, je, ks, ke;
		grid->GetRange(is, ie, js, je, ks, ke);
		RiemannSolverPara riemann_para[6];
		// 差分模板的值
		double res_tmp[5];
		double coef1 = 27.0 / 24.0;
		double coef2 = -1.0 / 24.0;
		double direction[3][3] = { {1,0,0},{0,1,0},{0,0,1} };
		for (int k = ks; k <= ke; ++k)
		{
#pragma omp parallel for private(res_tmp, riemann_para)
			for (int j = js; j <= je; ++j)
			{
				for (int i = is; i <= ie; ++i)
				{
					int idx = idx_proxy->GetIdx(i, j, k);
					for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
					{
						res_tmp[idx_eq] = data_manager->GetResidual(idx_eq, idx);
					}
					for (int dim = 0; dim < grid->GetDim(); dim++)
					{
						for (int iTemp = 0; iTemp < 4; iTemp++)
						{
							riemann_para[iTemp].gamma_left = riemann_para[iTemp].gamma_right = gas->GetGamma();
							riemann_para[iTemp].norm(0) = node_metrics->GetMetrics(dim, idx)[0];
							riemann_para[iTemp].norm(1) = node_metrics->GetMetrics(dim, idx)[1];
							riemann_para[iTemp].norm(2) = node_metrics->GetMetrics(dim, idx)[2];
							riemann_para[iTemp].nt = node_metrics->GetMetrics(dim, idx)[3];
							int idx_temp = idx_proxy->GetIdx(
								i + (iTemp - 2) * direction[dim][0],
								j + (iTemp - 2) * direction[dim][1],
								k + (iTemp - 2) * direction[dim][2]);
							for (int iVal = 0; iVal < equ_num; ++iVal)
							{
								riemann_para[iTemp].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, dim, idx_temp);
								riemann_para[iTemp].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, dim, idx_temp);
							}
							m_riemann_solver->Solver(riemann_para[iTemp]);
						}
						for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
						{
							res_tmp[idx_eq] -=
								coef1 * (riemann_para[2].flux[idx_eq] - riemann_para[1].flux[idx_eq]) +
								coef2 * (riemann_para[3].flux[idx_eq] - riemann_para[0].flux[idx_eq]);
						}
					}
					data_manager->SetResidual(idx, res_tmp);
				}
			}
		}
	}

	void NSSolverStructDEER::CalcFluxDifference6thOrder()
	{
		auto grid = GetGrid();
		auto para = GetPara();
		auto node = grid->GetNode();
		auto data_manager = GetDataManager();
		auto idx_proxy = GetIdxProxy();
		auto node_metrics = GetNodeMetrics();
		auto gas = GetGas();
		auto equ_num = para->GetEqNum();
		int is, ie, js, je, ks, ke;
		grid->GetRange(is, ie, js, je, ks, ke);
		RiemannSolverPara riemann_para[6];
		// 差分模板的值
		double res_tmp[5];
		double coef1 = 75.0 / 64.0;
		double coef2 = -25.0 / 384.0;
		double coef3 = 3.0 / 640.0;
		double direction[3][3] = { {1,0,0},{0,1,0},{0,0,1} };
		for (int k = ks; k <= ke; ++k)
		{
#pragma omp parallel for private(  res_tmp, riemann_para)
			for (int j = js; j <= je; ++j)
			{
				for (int i = is; i <= ie; ++i)
				{
					int idx = idx_proxy->GetIdx(i, j, k);
					for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
					{
						res_tmp[idx_eq] = data_manager->GetResidual(idx_eq, idx);
					}
					for (int dim = 0; dim < grid->GetDim(); dim++)
					{
						for (int iTemp = 0; iTemp < 6; iTemp++)
						{
							riemann_para[iTemp].gamma_left = riemann_para[iTemp].gamma_right = gas->GetGamma();
							riemann_para[iTemp].norm(0) = node_metrics->GetMetrics(dim, idx)[0];
							riemann_para[iTemp].norm(1) = node_metrics->GetMetrics(dim, idx)[1];
							riemann_para[iTemp].norm(2) = node_metrics->GetMetrics(dim, idx)[2];
							riemann_para[iTemp].nt = node_metrics->GetMetrics(dim, idx)[3];
							int idx_temp = idx_proxy->GetIdx(
								i + (iTemp - 3) * direction[dim][0],
								j + (iTemp - 3) * direction[dim][1],
								k + (iTemp - 3) * direction[dim][2]);
							for (int iVal = 0; iVal < equ_num; ++iVal)
							{
								riemann_para[iTemp].prim_left(iVal) = data_manager->GetMidNodePrimLeft(iVal, dim, idx_temp);
								riemann_para[iTemp].prim_right(iVal) = data_manager->GetMidNodePrimRight(iVal, dim, idx_temp);
							}
							m_riemann_solver->Solver(riemann_para[iTemp]);
						}
						for (int idx_eq = 0; idx_eq < equ_num; ++idx_eq)
						{
							res_tmp[idx_eq] -=
								coef1 * (riemann_para[3].flux[idx_eq] - riemann_para[2].flux[idx_eq]) +
								coef2 * (riemann_para[4].flux[idx_eq] - riemann_para[1].flux[idx_eq]) +
								coef3 * (riemann_para[5].flux[idx_eq] - riemann_para[0].flux[idx_eq]);
						}
					}
					data_manager->SetResidual(idx, res_tmp);
				}
			}
		}
	}
	void NSSolverStructDEER::MidNodeGrad(int idx_left, int idx_right, const double* lef_coord, const double* mid_coord,
		const double* right_coord, double* value_left, double* value_right)
	{
		int equ_num = GetPara()->GetEqNum();
		auto data_manager = GetDataManager();

		for (int iEqu = 0; iEqu < equ_num; ++iEqu)
		{
			value_left[iEqu] = data_manager->GetPrim(iEqu, idx_left);
			value_right[iEqu] = data_manager->GetPrim(iEqu, idx_right);

			for (int iDim = 0; iDim < 3; ++iDim)
			{
				double left_grad =
					data_manager->GetPrimGrad(iEqu, iDim, idx_left) * data_manager->GetLimiter(iEqu, idx_left);
				double right_grad =
					data_manager->GetPrimGrad(iEqu, iDim, idx_right) * data_manager->GetLimiter(iEqu, idx_right);

				value_left[iEqu] += (mid_coord[iDim] - lef_coord[iDim]) * left_grad;
				value_right[iEqu] += (mid_coord[iDim] - right_coord[iDim]) * right_grad;
			}
		}
	}
} // namespace zaran