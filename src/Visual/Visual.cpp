#include "Visual.h"
#include "GlobalData.h"
#include <fstream>
#include <string>
namespace zaran
{
	void Visual::WriteTecplot(std::shared_ptr<FieldSolver>& solver)
	{
		int step = GlobalData::GetInt("step");
		std::string filename = "result/" + std::to_string(step) + ".dat";
		std::ofstream fout(filename);
		fout << "variables=x,y,z,rho,u,v,w,p\n";
		auto& grid = solver->GetGrid();
		int nInnerNum = grid->GetInnerNodeNum();
		int nBoundNum = grid->GetBoundNodeNum();
		int nTotalNum = grid->GetTotalNodeNum();
		auto& nodeTopo = grid->GetNodeTopoInfo();
		auto& cellTopo = grid->GetCellTopoInfo();
		auto& data = solver->GetFieldData();
		auto& rho = data->GetData("rho");
		auto& u = data->GetData("u");
		auto& v = data->GetData("v");
		auto& w = data->GetData("w");
		auto& p = data->GetData("p");
		fout << "ZONE T= " << grid->GetName() << std::endl;
		fout << "N=" << nInnerNum + nBoundNum << ", E= " << cellTopo.size() << ", F=FEPOINT, ET=Brick" << std::endl;
		fout << "solutiontime= " << GlobalData::GetDouble("globalTime") << std::endl;
		for (int iNode = 0; iNode < nInnerNum + nBoundNum; ++iNode)
		{
			auto& currentCoord = nodeTopo[iNode].GetCoordinate();
			fout << currentCoord(0) << "  " << currentCoord(1) << "  " << currentCoord(2) << "  ";
			fout << rho[iNode] << "  " << u[iNode] << "  " << v[iNode] << "  " << w[iNode] << "  " << p[iNode] << "   ";
			fout << std::endl;
		}
		int nCell = cellTopo.size();
		for (int iCell = 0; iCell < nCell; ++iCell)
		{
			auto& cell2node = cellTopo[iCell].GetNode();
			for (int iNode = 0; iNode < cell2node.size(); ++iNode)
			{
				fout << cell2node[iNode] + 1 << "  ";
			}
			fout << std::endl;
		}
	}


	void Visual::WriteVTK(std::shared_ptr<FieldSolver>& solver)
	{
		// TO DO
	}

}