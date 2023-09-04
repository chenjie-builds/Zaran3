#include "Visual.h"
#include "GlobalData.h"
#include <fstream>
#include <string>
using namespace zaran;
void Visual::WriteTecplot(Ptr<FieldSolver>& solver)
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
void zaran::Visual::WriteTecplotPoint(Ptr<FieldSolver>& solver)
{
	int step = GlobalData::GetInt("step");
	std::string filename = "result/" + std::to_string(step) + ".dat";
	std::ofstream fout(filename);
	fout << "variables=x,y,z,rho,u,v,w,p\n";
	auto& grid = solver->GetGrid();
	auto& nodeTopo = grid->GetNodeTopoInfo();
	auto& cellTopo = grid->GetCellTopoInfo();
	auto& data = solver->GetFieldData();
	auto& rho = data->GetData("rho");
	auto& u = data->GetData("u");
	auto& v = data->GetData("v");
	auto& w = data->GetData("w");
	auto& p = data->GetData("p");
	//fout << "solutiontime= " << GlobalData::GetDouble("globalTime") << std::endl;
	for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
	{
		auto& currentCoord = nodeTopo[iNode].GetCoordinate();
		fout << currentCoord(0) << "  " << currentCoord(1) << "  " << currentCoord(2) << "  ";
		fout << rho[iNode] << "  " << u[iNode] << "  " << v[iNode] << "  " << w[iNode] << "  " << p[iNode];
		fout << std::endl;
	}
	fout.close();

	auto& boundMap = grid->GetBoundaryMap();
	//fout.open("result/inlet.dat");
	//fout << "variables=x,y,z\n";
	//auto& boundNode = boundMap->GetBoundary("inlet");
	//for (int iBound = 0; iBound < boundNode.size(); ++iBound)
	//{
	//	auto& boundIndex = boundNode[iBound].GetIndex();
	//	auto& currentNode = nodeTopo[boundIndex];
	//	auto& currentCoord = currentNode.GetCoordinate();
	//	fout << currentCoord(0) << "  " << currentCoord(1) << "  " << currentCoord(2) << "  ";
	//	fout << std::endl;
	//}
	//fout.close();
	//fout.open("result/outlet.dat");
	//fout << "variables=x,y,z\n";
	//boundNode = boundMap->GetBoundary("outlet");
	//for (int iBound = 0; iBound < boundNode.size(); ++iBound)
	//{
	//	auto& boundIndex = boundNode[iBound].GetIndex();
	//	auto& currentNode = nodeTopo[boundIndex];
	//	auto& currentCoord = currentNode.GetCoordinate();
	//	fout << currentCoord(0) << "  " << currentCoord(1) << "  " << currentCoord(2) << "  ";
	//	fout << std::endl;
	//}
	//fout.close();
	fout.open("result/wall_" + std::to_string(step) + ".dat");
	fout << "variables=x,y,z,rho,u,v,w,p\n";
	auto& boundNode = boundMap->GetBoundary("slipWall");
	//std::ofstream fout1("result/innerNode.dat");
	//fout1 << "variables=x,y,z\n";
	for (int iBound = 0; iBound < boundNode.size(); ++iBound)
	{
		auto& boundIndex = boundNode[iBound].GetIndex();
		//auto& innerIndex = boundNode[iBound].GetInnerNodeIndex();
		auto& currentNode = nodeTopo[boundIndex];
		auto& currentCoord = currentNode.GetCoordinate();
		fout << currentCoord(0) << "  " << currentCoord(1) << "  " << currentCoord(2)<<"  ";
		fout << rho[boundIndex] << "  " << u[boundIndex] << "  " << v[boundIndex] << "  " << w[boundIndex] << "  " << p[boundIndex];
		fout << std::endl;
		//auto& innerNode = nodeTopo[innerIndex];
		//auto& innerNodeCoord = innerNode.GetCoordinate();
		//fout1 << innerNodeCoord(0) << "  " << innerNodeCoord(1) << "  " << innerNodeCoord(2);
		//fout1 << std::endl;
	}
	fout.close();
	//fout1.close();



}

void Visual::WriteVTK(Ptr<FieldSolver>& solver)
{
	// TO DO
}

