#include "Visual.h"
#include "GlobalData.h"
#include "Grid_Zaran_3D.h"
#include <fstream>
#include <string>
using namespace zaran;
void Visual::WriteTecplot(Ptr<FieldSolver>& solver)
{
	int step = GlobalData::GetInt("step");
	std::string filename = "result/" + std::to_string(step) + ".dat";
	std::ofstream fout(filename);
	fout << "variables=x,y,z,rho,u,v,w,p,jacobi\n";
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
	auto& jacobi = data->GetData("coordTransJ");
	fout << "ZONE T= grid_" << grid->GetName() << std::endl;
	fout << "N=" << grid->GetTotalNodeNum() << ", E= " << cell2node.size() << ", F=FEPOINT, ET=Brick" << std::endl;
	fout << "solutiontime= " << GlobalData::GetDouble("globalTime") << std::endl;
	for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
	{
		auto& currentCoord = nodeCoord[iNode];
		fout << currentCoord(0) << "  " << currentCoord(1) << "  " << currentCoord(2) << "  ";
		fout << rho[iNode] << "  " << u[iNode] << "  " << v[iNode] << "  " << w[iNode] << "  " << p[iNode] << "  " << jacobi[iNode];
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
	fout << "variables=x,y,rho,u,v,w,p,J,lim0,lim1,lim2,lim3,lim4,rhoGradX,rhoGradY,uGradX,uGradY,vGradX,vGradY,pGradX,pGradY\n";
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
	auto& jacobi = data->GetData("coordTransJ");
	auto& rhoGradX = data->GetData("rhoGradX");
	auto& rhoGradY = data->GetData("rhoGradY");
	auto& uGradX = data->GetData("uGradX");
	auto& uGradY = data->GetData("uGradY");
	auto& vGradX = data->GetData("vGradX");
	auto& vGradY = data->GetData("vGradY");
	auto& pGradX = data->GetData("pGradX");
	auto& pGradY = data->GetData("pGradY");
	auto& limiterCoef0 = data->GetData("limiterCoef0");
	auto& limiterCoef1 = data->GetData("limiterCoef1");
	auto& limiterCoef2 = data->GetData("limiterCoef2");
	auto& limiterCoef3 = data->GetData("limiterCoef3");
	auto& limiterCoef4 = data->GetData("limiterCoef4");
	fout << "ZONE T= grid_" << grid->GetName() << std::endl;
	fout << "N=" << grid->GetTotalNodeNum() << ", E= " << cell2node.size() << ", F=FEPOINT, ET=QUADRILATERAL" << std::endl;
	fout << "solutiontime= " << GlobalData::GetDouble("globalTime") << std::endl;
	for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
	{
		auto& currentCoord = nodeCoord[iNode];
		fout << currentCoord(0) << "  " << currentCoord(1) << "  ";
		fout << rho[iNode] << "  " << u[iNode] << "  " << v[iNode] << "  " << w[iNode] << "  " << p[iNode];
		fout << "  " << jacobi[iNode];
		fout << "  " << limiterCoef0[iNode] << "  " << limiterCoef1[iNode] << "  " << limiterCoef2[iNode] << "  " << limiterCoef3[iNode] << "  " << limiterCoef4[iNode];
		fout << "  " << rhoGradX[iNode] << "  " << rhoGradY[iNode];
		fout << "  " << uGradX[iNode] << "  " << uGradY[iNode];
		fout << "  " << vGradX[iNode] << "  " << vGradY[iNode];
		fout << "  " << pGradX[iNode] << "  " << pGradY[iNode];
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

void zaran::Visual::WriteTecplotZaran3D(Ptr<FieldSolver>& solver)
{
	int step = GlobalData::GetInt("step");
	std::string filename = "result/" + std::to_string(step) + ".dat";
	std::ofstream fout(filename);
	fout << "variables=x,y,z,label\n";
	Ptr<Grid_Zaran_3D>& grid = std::static_pointer_cast<Grid_Zaran_3D>(solver->GetGrid());
	auto cellTopo = std::static_pointer_cast<CellTopoInfoZaran>(grid->GetCellTopo());
	auto cell_type = cellTopo->GetType();
	int ni, nj, nk;
	grid->GetNodeNum(ni, nj, nk);
	double x_min, x_max, y_min, y_max, z_min, z_max;
	grid->GetBox(x_min, x_max, y_min, y_max, z_min, z_max);
	double dx = (x_max - x_min) / (ni - 1);
	double dy = (y_max - y_min) / (nj - 1);
	double dz = (z_max - z_min) / (nk - 1);
	fout << "ZONE T= grid_" << grid->GetName() << std::endl;
	fout << "I=" << ni << ", J=" << nj << ", K=" << nk << ", DATAPACKING=BLOCK, VARLOCATION=([4]CELLCENTERED)" << std::endl;
	fout << "solutiontime= " << GlobalData::GetDouble("globalTime") << std::endl;
	int count = 0;
	for (int k = 0; k < nk - 1; ++k)
	{
		for (int j = 0; j < nj - 1; ++j)
		{
			for (int i = 0; i < ni - 1; ++i)
			{
				fout << x_min + i * dx + 0.5 * dx << " ";
				if (++count % 10 == 0)
				{
					count = 0;
					fout << std::endl;
				}
			}
		}
	}
	for (int k = 0; k < nk - 1; ++k)
	{
		for (int j = 0; j < nj - 1; ++j)
		{
			for (int i = 0; i < ni - 1; ++i)
			{
				fout << y_min + j * dy + 0.5 * dy << " ";
				if (++count % 10 == 0)
				{
					count = 0;
					fout << std::endl;
				}
			}
		}
	}
	for (int k = 0; k < nk - 1; ++k)
	{
		for (int j = 0; j < nj - 1; ++j)
		{
			for (int i = 0; i < ni - 1; ++i)
			{
				fout << z_min + k * dz + 0.5 * dz << " ";
				if (++count % 10 == 0)
				{
					count = 0;
					fout << std::endl;
				}
			}
		}
	}
	for (int k = 0; k < nk - 1; ++k)
	{
		for (int j = 0; j < nj - 1; ++j)
		{
			for (int i = 0; i < ni - 1; ++i)
			{
				fout << int(cell_type[grid->GetCellIndex(i, j, k)]) << " ";
				if (++count % 10 == 0)
				{
					count = 0;
					fout << std::endl;
				}
			}
		}
	}
	fout.close();

}

