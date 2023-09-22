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
	auto& nodeTopo = grid->GetNodeTopo();
	auto& nodeCoord = nodeTopo->GetCoordinate();
	auto& cellTopo = grid->GetCellTopo();
	auto& cell2node = cellTopo->GetNodeIndex();
	auto& data = solver->GetFieldData();
	auto& rho = data->GetData("rho");
	auto& u = data->GetData("u");
	auto& v = data->GetData("v");
	auto& w = data->GetData("w");
	auto& p = data->GetData("p");
	fout << "ZONE T= grid_" << grid->GetName() << std::endl;
	fout << "N=" << grid->GetTotalNodeNum() << ", E= " << cell2node.size() << ", F=FEPOINT, ET=Brick" << std::endl;
	fout << "solutiontime= " << GlobalData::GetDouble("globalTime") << std::endl;
	for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
	{
		auto& currentCoord = nodeCoord[iNode];
		fout << currentCoord(0) << "  " << currentCoord(1) << "  " << currentCoord(2) << "  ";
		fout << rho[iNode] << "  " << u[iNode] << "  " << v[iNode] << "  " << w[iNode] << "  " << p[iNode];
		fout << std::endl;
	}
	int nCell = cell2node.size();
	for (int iCell = 0; iCell < nCell; ++iCell)
	{
		for (int iNode = 0; iNode < cell2node[iCell].size(); ++iNode)
		{
			fout << cell2node[iCell][iNode] + 1 << "  ";
		}
		fout << std::endl;
	}
	fout.close();
}
void zaran::Visual::WriteTecplot2D(Ptr<FieldSolver>& solver)
{
	int step = GlobalData::GetInt("step");
	std::string filename = "result/" + std::to_string(step) + ".dat";
	std::ofstream fout(filename);
	fout << "variables=x,y,rho,u,v,p\n";
	auto& grid = solver->GetGrid();
	int nInnerNum = grid->GetInnerNodeNum();
	int nBoundNum = grid->GetBoundNodeNum();
	int nTotalNum = grid->GetTotalNodeNum();
	auto& nodeTopo = grid->GetNodeTopo();
	auto& nodeCoord = nodeTopo->GetCoordinate();
	auto& cellTopo = grid->GetCellTopo();
	auto& cell2node = cellTopo->GetNodeIndex();
	auto& data = solver->GetFieldData();
	auto& rho = data->GetData("rho");
	auto& u = data->GetData("u");
	auto& v = data->GetData("v");
	auto& w = data->GetData("w");
	auto& p = data->GetData("p");
	fout << "ZONE T= grid_" << grid->GetName() << std::endl;
	fout << "N=" << grid->GetTotalNodeNum() << ", E= " << cell2node.size() << ", F=FEPOINT, ET=QUADRILATERAL" << std::endl;
	fout << "solutiontime= " << GlobalData::GetDouble("globalTime") << std::endl;
	for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
	{
		auto& currentCoord = nodeCoord[iNode];
		fout << currentCoord(0) << "  " << currentCoord(1) << "  ";
		fout << rho[iNode] << "  " << u[iNode] << "  " << v[iNode] << "  " << p[iNode];
		fout << std::endl;
	}
	int nCell = cell2node.size();
	for (int iCell = 0; iCell < nCell; ++iCell)
	{
		for (int iNode = 0; iNode < cell2node[iCell].size(); ++iNode)
		{
			fout << cell2node[iCell][iNode] + 1 << "  ";
		}
		fout << std::endl;
	}
	fout.close();
}
void zaran::Visual::WriteTecplotPoint(Ptr<FieldSolver>& solver)
{
	int step = GlobalData::GetInt("step");
	std::string filename = "result/" + std::to_string(step) + ".dat";
	std::ofstream fout(filename);
	fout << "variables=x,y,z,rho,u,v,w,p\n";
	auto& grid = solver->GetGrid();
	auto& nodeTopo = grid->GetNodeTopo();
	auto& nodeCoord = nodeTopo->GetCoordinate();
	auto& data = solver->GetFieldData();
	auto& rho = data->GetData("rho");
	auto& u = data->GetData("u");
	auto& v = data->GetData("v");
	auto& w = data->GetData("w");
	auto& p = data->GetData("p");
	//fout << "solutiontime= " << GlobalData::GetDouble("globalTime") << std::endl;
	for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
	{
		auto& currentCoord = nodeCoord[iNode];
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



	//fout.open("result/wall_" + std::to_string(step) + ".dat");
	//fout << "variables=x,y,z,rho,u,v,w,p\n";
	//auto& wallBoundNode = boundMap->GetBoundary("slipWall");
	//for (int iBound = 0; iBound < wallBoundNode.size(); ++iBound)
	//{
	//	auto& boundIndex = wallBoundNode[iBound].GetIndex();
	//	auto& currentNode = nodeTopo[boundIndex];
	//	auto& currentCoord = currentNode.GetCoordinate();
	//	fout << currentCoord(0) << "  " << currentCoord(1) << "  " << currentCoord(2) << "  ";
	//	fout << rho[boundIndex] << "  " << u[boundIndex] << "  " << v[boundIndex] << "  " << w[boundIndex] << "  " << p[boundIndex];
	//	fout << std::endl;
	//}
	//fout.close();
	//fout.open("result/wallInnerNode_" + std::to_string(step) + ".dat");
	//fout << "variables=x,y,z,rho,u,v,w,p\n";
	//for (int iBound = 0; iBound < wallBoundNode.size(); ++iBound)
	//{
	//	auto& innerIndex = wallBoundNode[iBound].GetInnerNodeIndex();
	//	auto& innerNode = nodeTopo[innerIndex];
	//	auto& innerNodeCoord = innerNode.GetCoordinate();
	//	fout << innerNodeCoord(0) << "  " << innerNodeCoord(1) << "  " << innerNodeCoord(2);
	//	fout << "  " << rho[innerIndex] << "  " << u[innerIndex] << "  " << v[innerIndex] << "  " << w[innerIndex] << "  " << p[innerIndex];
	//	fout << std::endl;
	//}
	//fout.close();

	//fout.open("result/outlet_" + std::to_string(step) + ".dat");
	//fout << "variables=x,y,z,rho,u,v,w,p\n";
	//auto& outletBoundNode = boundMap->GetBoundary("outlet");
	//for (int iBound = 0; iBound < outletBoundNode.size(); ++iBound)
	//{
	//	auto& boundIndex = outletBoundNode[iBound].GetIndex();
	//	auto& currentNode = nodeTopo[boundIndex];
	//	auto& currentCoord = currentNode.GetCoordinate();
	//	fout << currentCoord(0) << "  " << currentCoord(1) << "  " << currentCoord(2) << "  ";
	//	fout << rho[boundIndex] << "  " << u[boundIndex] << "  " << v[boundIndex] << "  " << w[boundIndex] << "  " << p[boundIndex];
	//	fout << std::endl;
	//}
	//fout.close();
	//fout.open("result/outletInnerNode_" + std::to_string(step) + ".dat");
	//fout << "variables=x,y,z,rho,u,v,w,p\n";
	//for (int iBound = 0; iBound < outletBoundNode.size(); ++iBound)
	//{
	//	auto& innerIndex = outletBoundNode[iBound].GetInnerNodeIndex();
	//	auto& innerNode = nodeTopo[innerIndex];
	//	auto& innerNodeCoord = innerNode.GetCoordinate();
	//	fout << innerNodeCoord(0) << "  " << innerNodeCoord(1) << "  " << innerNodeCoord(2);
	//	fout << "  " << rho[innerIndex] << "  " << u[innerIndex] << "  " << v[innerIndex] << "  " << w[innerIndex] << "  " << p[innerIndex];
	//	fout << std::endl;
	//}
	//fout.close();

	//fout.open("result/inner_" + std::to_string(step) + ".dat");
	//fout << "variables=x,y,z,rho,u,v,w,p\n";
	//for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
	//{
	//	auto&currentNode= nodeTopo[iNode];
	//	auto& currentCoord = nodeTopo[iNode].GetCoordinate();
	//	if(currentNode.GetType()!=NodeType::inner)
	//		continue;
	//	fout << currentCoord(0) << "  " << currentCoord(1) << "  " << currentCoord(2) << "  ";
	//	fout << rho[iNode] << "  " << u[iNode] << "  " << v[iNode] << "  " << w[iNode] << "  " << p[iNode];
	//	fout << std::endl;
	//}
	//fout.close();

	//fout.open("result/unset_" + std::to_string(step) + ".dat");
	//fout << "variables=x,y,z,rho,u,v,w,p\n";
	//for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
	//{
	//	auto&currentNode = nodeTopo[iNode];
	//	auto& currentCoord = nodeTopo[iNode].GetCoordinate();
	//	if (currentNode.GetType() != NodeType::undefined)
	//		continue;
	//	fout << currentCoord(0) << "  " << currentCoord(1) << "  " << currentCoord(2) << "  ";
	//	fout << rho[iNode] << "  " << u[iNode] << "  " << v[iNode] << "  " << w[iNode] << "  " << p[iNode];
	//	fout << std::endl;
	//}


}

void Visual::WriteVTK(Ptr<FieldSolver>& solver)
{
	// TO DO
}

