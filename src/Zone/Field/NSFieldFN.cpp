#include "NSFieldFN.h"

namespace zaran
{
	NSFieldFNFDM::NSFieldFNFDM(shared_ptr<GridFN> grid)
		: FieldNS(grid, FieldType::NS_FlexibleNode)
	{

	}

	NSFieldFNFDM::~NSFieldFNFDM()
	{

	}

	shared_ptr<zaran::GridFN> NSFieldFNFDM::GetGrid()
	{
		return std::static_pointer_cast<GridFN>(FieldNS::GetGrid());
	}

	shared_ptr<zaran::NSSolverFNFDM> NSFieldFNFDM::GetSolver()
	{
		return std::static_pointer_cast<NSSolverFNFDM>(FieldNS::GetSolver());
	}

	void NSFieldFNFDM::CalcResidual()
	{
		auto grid = GetGrid();
		auto& node = grid->GetNode();
		auto para = GetSolverPara();
		int equ_num = para->GetEqNum();
		int inner_node_num = grid->GetInnerNodeNum();
		auto inner_node = grid->GetInnerNode();
		double norm_inf = -LARGE_NUMBER;
		double norm_L2 = 0;
		int norm_inf_node = 0;
		double norm_inf_coord[3];
		for (int iEqu = 0; iEqu < equ_num; iEqu++)
		{
			norm_inf = -LARGE_NUMBER;
			norm_L2 = 0;
			norm_inf_node = -1;
			auto res = GetDataManager()->GetResidual(iEqu);
#pragma omp parallel for reduction(max:norm_inf) reduction(+:norm_L2)
			for (int idx = 0; idx < inner_node_num; idx++)
			{
				int iNode = inner_node[idx];
				if (abs(res[iNode]) > norm_inf)
				{
					norm_inf = abs(res[iNode]);
					norm_inf_node = iNode;
					for (dimension_type iDim = 0; iDim < grid->GetDim(); iDim++)
					{
						norm_inf_coord[iDim] = node.GetCoord(iNode)[iDim];
					}
				}
				norm_L2 += res[iNode] * res[iNode];
			}
			norm_L2 = sqrt(norm_L2 / inner_node_num);
			m_res_info->SetInfNorm(iEqu, norm_inf);
			m_res_info->SetL2Norm(iEqu, norm_L2);
			m_res_info->SetInfNormCoord(iEqu, norm_inf_coord);
		}
	}

	void NSFieldFNFDM::AllocateSolver()
	{
		m_solver = make_shared<NSSolverFNFDM>(1, "NS_FNFDM", GetSolverPara(), GetGrid(), GetDataManager());
	}
	void NSFieldFNFDM::AllocateDataManager()
	{
		m_data_manager = make_shared<DataManagerNS>(GetData(), GetGrid()->GetTotalNodeNum());
		m_data_manager->CreateData();
		m_data_manager->RegisterData();
	}
} // namespace zaran