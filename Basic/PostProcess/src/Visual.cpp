#include "Visual.h"
#include "GlobalData.h"
#include "Grid_Zaran_3D.h"
#include<TECIO.h>
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
	fout << "variables=x,y,z,rho,u,v,w,p\n";
	auto& data = solver->GetFieldData();
	auto& rho = data->GetData("rho");
	auto& u = data->GetData("u");
	auto& v = data->GetData("v");
	auto& w = data->GetData("w");
	auto& p = data->GetData("p");
	auto& jacobi = data->GetData("coordTransJ");
	Ptr<Grid_Zaran_3D>& grid = std::static_pointer_cast<Grid_Zaran_3D>(solver->GetGrid());
	auto& cellTopo = grid->GetCellTopo();
	auto& cell_type = cellTopo->GetType();
	auto& cell_center = cellTopo->GetCenterCoord();
	auto& node_topo = grid->GetNodeTopo();
	auto& node_coord = node_topo->GetCoordinate();
	int is, ie, js, je, ks, ke;
	grid->GetRange(is, ie, js, je, ks, ke);
	fout << "ZONE T= grid_" << grid->GetName() << std::endl;
	fout << "I=" << ie - is << ", J=" << je - js << ", K=" << ke - ks << ", DATAPACKING=BLOCK, VARLOCATION=([4-8]=CELLCENTERED)" << std::endl;
	fout << "solutiontime= " << GlobalData::GetDouble("globalTime") << std::endl;
	auto CellIndex = [&](int i, int j, int k) {return grid->GetCellIndex(i, j, k); };
	int iCell = 0;
	int count = 0;
	for (int k = ks; k < ke; ++k)
	{
		for (int j = js; j < je; ++j)
		{
			for (int i = is; i < ie; ++i)
			{
				fout << node_coord[grid->GetNodeIndex(i, j, k)].x() << " ";
				if (++count % 10 == 0)
				{
					count = 0;
					fout << std::endl;
				}
			}
		}
	}
	for (int k = ks; k < ke; ++k)
	{
		for (int j = js; j < je; ++j)
		{
			for (int i = is; i < ie; ++i)
			{
				fout << node_coord[grid->GetNodeIndex(i, j, k)].y() << " ";
				if (++count % 10 == 0)
				{
					count = 0;
					fout << std::endl;
				}
			}
		}
	}
	for (int k = ks; k < ke; ++k)
	{
		for (int j = js; j < je; ++j)
		{
			for (int i = is; i < ie; ++i)
			{
				fout << node_coord[grid->GetNodeIndex(i, j, k)].z() << " ";
				if (++count % 10 == 0)
				{
					count = 0;
					fout << std::endl;
				}
			}
		}
	}
	for (int k = ks; k < ke - 1; ++k)
	{
		for (int j = js; j < je - 1; ++j)
		{
			for (int i = is; i < ie - 1; ++i)
			{
				iCell = CellIndex(i, j, k);
				fout << rho[iCell] << " ";
				if (++count % 10 == 0)
				{
					count = 0;
					fout << std::endl;
				}
			}
		}
	}
	for (int k = ks; k < ke - 1; ++k)
	{
		for (int j = js; j < je - 1; ++j)
		{
			for (int i = is; i < ie - 1; ++i)
			{
				iCell = CellIndex(i, j, k);
				fout << u[iCell] << " ";
				if (++count % 10 == 0)
				{
					count = 0;
					fout << std::endl;
				}
			}
		}
	}
	for (int k = ks; k < ke - 1; ++k)
	{
		for (int j = js; j < je - 1; ++j)
		{
			for (int i = is; i < ie - 1; ++i)
			{
				iCell = CellIndex(i, j, k);
				fout << v[iCell] << " ";
				if (++count % 10 == 0)
				{
					count = 0;
					fout << std::endl;
				}
			}
		}
	}
	for (int k = ks; k < ke - 1; ++k)
	{
		for (int j = js; j < je - 1; ++j)
		{
			for (int i = is; i < ie - 1; ++i)
			{
				iCell = CellIndex(i, j, k);
				fout << w[iCell] << " ";
				if (++count % 10 == 0)
				{
					count = 0;
					fout << std::endl;
				}
			}
		}
	}
	for (int k = ks; k < ke - 1; ++k)
	{
		for (int j = js; j < je - 1; ++j)
		{
			for (int i = is; i < ie - 1; ++i)
			{
				iCell = CellIndex(i, j, k);
				fout << jacobi[iCell] << " ";
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

void zaran::Visual::WriteTecplotZaran3DBinary(Ptr<FieldSolver>& solver)
{
	auto& data = solver->GetFieldData();
	auto& rho = data->GetData("rho");
	auto& u = data->GetData("u");
	auto& v = data->GetData("v");
	auto& w = data->GetData("w");
	auto& p = data->GetData("p");
	auto& jacobi = data->GetData("coordTransJ");
	Ptr<Grid_Zaran_3D>& grid = std::static_pointer_cast<Grid_Zaran_3D>(solver->GetGrid());
	auto& cellTopo = grid->GetCellTopo();
	auto& cell_type = cellTopo->GetType();
	auto& cell_center = cellTopo->GetCenterCoord();
	auto& node_topo = grid->GetNodeTopo();
	auto& node_coord = node_topo->GetCoordinate();

	int is, ie, js, je, ks, ke;
	grid->GetRange(is, ie, js, je, ks, ke);
	int n_node = (ie - is) * (je - js) * (ke - ks);
	int n_cell = (ie - is - 1) * (je - js - 1) * (ke - ks - 1);
	DArray x(n_node), y(n_node), z(n_node), density(n_cell), x_vel(n_cell), y_vel(n_cell), z_vel(n_cell), pressure(n_cell);
	int step = GlobalData::GetInt("step");
	std::string file_name = "result/" + std::to_string(step) + ".plt";
	string zone_name = "grid_" + grid->GetName();
	string var_name = "x y z rho u v w p";
	double solution_time = GlobalData::GetDouble("globalTime");
	int node_index = 0;
	for (int k = ks; k < ke; ++k)
	{
		for (int j = js; j < je; ++j)
		{
			for (int i = is; i < ie; ++i)
			{
				node_index=(k - ks) * (je - js) * (ie - is) + (j - js) * (ie - is) + (i - is);
				x[node_index] = node_coord[grid->GetNodeIndex(i, j, k)].x();
				y[node_index] = node_coord[grid->GetNodeIndex(i, j, k)].y();
				z[node_index] = node_coord[grid->GetNodeIndex(i, j, k)].z();
			}
		}
	}
	int cell_index = 0;
	for (int k = ks; k < ke - 1; ++k)
	{
		for (int j = js; j < je - 1; ++j)
		{
			for (int i = is; i < ie - 1; ++i)
			{
				cell_index = (k - ks) * (je - js - 1) * (ie - is - 1) + (j - js ) * (ie - is - 1) + (i - is);
				density[cell_index] = rho[grid->GetCellIndex(i, j, k)];
				x_vel[cell_index] = u[grid->GetCellIndex(i, j, k)];
				y_vel[cell_index] = v[grid->GetCellIndex(i, j, k)];
				z_vel[cell_index] = w[grid->GetCellIndex(i, j, k)];
				pressure[cell_index] = p[grid->GetCellIndex(i, j, k)];
			}
		}
	}

	int value_location[8] = { 1,1,1,0,0,0,0,0 };
	INTEGER4 file_format = 0;
	INTEGER4 Debug = 0;
	INTEGER4 VIsDouble = 1;
	INTEGER4 FileType = 0;

	int i = TECINI142((char*)zone_name.c_str(),
		(char*)var_name.c_str(),
		(char*)file_name.c_str(),
		(char*)".",
		&file_format,
		&FileType,
		&Debug,
		&VIsDouble);
	INTEGER4 ZoneType = 0;
	INTEGER4 IMax = ie - is;
	INTEGER4 JMax = je - js;
	INTEGER4 KMax = ke - ks;
	INTEGER4 ICellMax = ie - is - 1;
	INTEGER4 JCellMax = je - js - 1;
	INTEGER4 KCellMax = ke - ks - 1;
	INTEGER4 StrandID = 0;      /* StaticZone */
	INTEGER4 unused = 0;      // ParentZone is no longer used
	INTEGER4 IsBlock = 1;      /* Block */
	INTEGER4 NFConns = 0;
	INTEGER4 FNMode = 0;
	INTEGER4 TotalNumFaceNodes = 1;
	INTEGER4 TotalNumBndryFaces = 1;
	INTEGER4 TotalNumBndryConnections = 1;
	INTEGER4 ShrConn = 0;

	i = TECZNE142((char*)zone_name.c_str(),
		&ZoneType,
		&IMax,
		&JMax,
		&KMax,
		&ICellMax,
		&JCellMax,
		&KCellMax,
		&solution_time,
		&StrandID,
		&unused,
		&IsBlock,
		&NFConns,
		&FNMode,
		&TotalNumFaceNodes,
		&TotalNumBndryFaces,
		&TotalNumBndryConnections,
		NULL,
		value_location,
		NULL,
		&ShrConn);
	i = TECDAT142(&n_node, x.data(), &VIsDouble);
	i = TECDAT142(&n_node, y.data(), &VIsDouble);
	i = TECDAT142(&n_node, z.data(), &VIsDouble);
	i = TECDAT142(&n_cell, density.data(), &VIsDouble);
	i = TECDAT142(&n_cell, x_vel.data(), &VIsDouble);
	i = TECDAT142(&n_cell, y_vel.data(), &VIsDouble);
	i = TECDAT142(&n_cell, z_vel.data(), &VIsDouble);
	i = TECDAT142(&n_cell, pressure.data(), &VIsDouble);

	i = TECEND142();
}

