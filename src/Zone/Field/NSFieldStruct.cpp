#include "NSFieldStruct.h"
#include "DataManagerNSStruct.h"
#include "FlowSolverStructPara.h"
#include "NSSolverStructCMM.h"
#include "NSSolverStructDEER.h"
#include "Log.h"
namespace zaran
{
	NSFieldStruct::NSFieldStruct(shared_ptr<GridStruct> grid)
		: FieldNS(grid, FieldType::NS_Structured)
	{
	}
	NSFieldStruct::~NSFieldStruct()
	{
	}

	void NSFieldStruct::CalcResidual()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto para = GetSolverPara();
		IdProxyStruct& idx_proxy = grid->GetIdxProxy();
		int equ_num = para->GetEqNum();
		auto data_manager = GetDataManager();
		int is, ie, js, je, ks, ke;
		grid->GetRange(is, ie, js, je, ks, ke);
		int total_node_num = (ie - is + 1) * (je - js + 1) * (ke - ks + 1);
		double norm_inf = -LARGE_NUMBER;
		double norm_l2 = 0;
		int norm_inf_node = -1;
		double norm_inf_coord[3];
		for (int iEqu = 0; iEqu < equ_num; iEqu++)
		{
			norm_inf = -LARGE_NUMBER;
			norm_l2 = 0;
			auto res = GetDataManager()->GetResidual(iEqu);
#pragma omp parallel for reduction(max : norm_inf) reduction(+ : norm_l2)
			for (int k = ks; k <= ke; k++)
			{
				for (int j = js; j <= je; j++)
				{
					for (int i = is; i <= ie; i++)
					{
						int idx =idx_proxy(i, j, k);
						double res_val = abs(res[idx]);
						auto coord = node->GetCoord(i, j, k);
						//if (iEqu == 0)
						//{
						//	double gamma = 1.4;
						//	double beta = 5.0;
						//	double xc = 5.0;
						//	double yc = 5.0;
						//	double r2 = (coord[0] - xc) * (coord[0] - xc) + (coord[1] - yc) * (coord[1] - yc);
						//	double density_num = data_manager->GetPrim(ID_DENSITY, idx);
						//	double density_exact = pow(1.0 - (gamma - 1.0) * beta * beta * exp(1.0 - r2) / (8.0 * gamma * PI * PI), 1.0 / (gamma - 1.0));
						//	res_val = abs(density_num - density_exact) / density_exact;
						//}
						//if (iEqu == 0)
						//{
						//	double density_num = data_manager->GetPrim(ID_DENSITY, idx);
						//	double density_exact = 1 - 0.005 * coord[0] + 0.01 * coord[1];
						//	res_val = abs(density_num - density_exact) / density_exact;
						//}
						if (res_val > norm_inf)
						{
							norm_inf = res_val;
							norm_inf_node = idx;
							for (dimension_type iDim = 0; iDim < grid->GetDim(); iDim++)
							{
								norm_inf_coord[iDim] = coord[iDim];
							}
						}
						norm_l2 += res_val * res_val;
					}
				}
			}
			norm_l2 = sqrt(norm_l2 / total_node_num);
			m_res_info->SetInfNorm(iEqu, norm_inf);
			m_res_info->SetL2Norm(iEqu, norm_l2);
			m_res_info->SetInfNormCoord(iEqu, norm_inf_coord);
			m_res_info->SetInfNormIdx(iEqu, norm_inf_node);
		}
	}
	void NSFieldStruct::Allocate()
	{
		FieldNS::Allocate();
	}

	shared_ptr<zaran::GridStruct> NSFieldStruct::GetGrid()
	{
		return std::static_pointer_cast<zaran::GridStruct>(Field::GetGrid());
	}

	shared_ptr<zaran::FlowSolverParamStruct> NSFieldStruct::GetSolverPara()
	{
		return std::static_pointer_cast<zaran::FlowSolverParamStruct>(Field::GetSolverPara());
	}

	shared_ptr<zaran::NSSolverStruct> NSFieldStruct::GetSolver()
	{
		return std::static_pointer_cast<zaran::NSSolverStruct>(Field::GetSolver());
	}

	shared_ptr<zaran::DataManagerNSStruct> NSFieldStruct::GetDataManager()
	{
		return std::static_pointer_cast<zaran::DataManagerNSStruct>(Field::GetDataManager());
	}

	void NSFieldStruct::AllocateResInfo()
	{
		m_res_info = make_shared<ResInfo>(GetSolverPara()->GetEqNum());
	}
	void NSFieldStruct::AllocateSolver()
	{
		auto para = GetSolverPara();
		if (para->GetMidMetricsScheme() == MidMetricsScheme::CMM)
		{
			m_solver = make_shared<NSSolverStructCMM>(GetIdx(), "NS_Struct", GetSolverPara(), GetGrid(), GetDataManager());
		}
		else if (para->GetMidMetricsScheme() == MidMetricsScheme::DEER)
		{
			m_solver = make_shared<NSSolverStructDEER>(GetIdx(), "NS_Struct", GetSolverPara(), GetGrid(), GetDataManager());
		}
		else
		{
			Log::error("Invalid metrics scheme");
			throw std::runtime_error("Invalid metrics scheme");
		}
	}
	void NSFieldStruct::AllocateDataManager()
	{
		int ni = GetGrid()->GetNi();
		int nj = GetGrid()->GetNj();
		int nk = GetGrid()->GetNk();
		m_data_manager = make_shared<DataManagerNSStruct>(GetData(), ni, nj, nk);
		m_data_manager->CreateData();
		m_data_manager->RegisterData();
	}
	void NSFieldStruct::AllocateSolverPara()
	{
		m_solver_para = make_shared<FlowSolverParamStruct>();
		GetSolverPara()->Init();
	}
} // namespace zaran