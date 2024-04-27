#include "Visual.h"
#include "GlobalData.h"
#include "Grid_Zaran_3D.h"
#include<TECIO.h>
#include <fstream>
#include <string>
#include"Log.h"
using namespace zaran;
void Visual::WriteTecplot(FieldSolver* solver)
{
	WriteTecplotBinary(solver);

	return;
	int currentIter = GlobalData::GetInt("currentIter");
	std::string filename = "result/" + std::to_string(currentIter) + ".dat";
	std::ofstream fout(filename);
	fout << "variables=x,y,z,rho,u,v,w,p,jacobi\n";
	Grid* grid = solver->GetGrid();
	int nInnerNum = grid->GetInnerNodeNum();
	int nBoundNum = grid->GetBoundNodeNum();
	int nTotalNum = grid->GetTotalNodeNum();
	NodeTopo* nodeTopo = grid->GetNodeTopo();

	auto& nodeCoord = nodeTopo->GetCoordinate();
	CellTopo* cellTopo = grid->GetCellTopo();
	auto& cell2node = cellTopo->GetNodeIndex();
	FieldData* data = solver->GetFieldData();
	double* density, * velocity_x, * velocity_y, * velocity_z, * pressure, * jacobi;
	data->GetData("density", density);
	data->GetData("velocity_x", velocity_x);
	data->GetData("velocity_y", velocity_y);
	data->GetData("velocity_w", velocity_z);
	data->GetData("pressure", pressure);
	data->GetData("coordTransJ", jacobi);
	fout << "ZONE T= grid_" << grid->GetName() << std::endl;
	fout << "N=" << grid->GetTotalNodeNum() << ", E= " << cell2node.size() << ", F=FEPOINT, ET=Brick" << std::endl;
	fout << "solutiontime= " << GlobalData::GetDouble("currentTime") << std::endl;
	for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
	{
		auto& currentCoord = nodeCoord[iNode];
		fout << currentCoord(0) << "  " << currentCoord(1) << "  " << currentCoord(2) << "  ";
		fout << density[iNode] << "  " << velocity_x[iNode] << "  " << velocity_y[iNode] << "  " << velocity_z[iNode] << "  " << pressure[iNode] << "  " << jacobi[iNode];
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
void zaran::Visual::WriteTecplotBinary(FieldSolver* solver)
{
	FieldData* data = solver->GetFieldData();
	double* density, * velocity_x, * velocity_y, * velocity_z, * pressure;
	data->GetData("density", density);
	data->GetData("velocity_x", velocity_x);
	data->GetData("velocity_y", velocity_y);
	data->GetData("velocity_w", velocity_z);
	data->GetData("pressure", pressure);
	Grid* grid = solver->GetGrid();
	CellTopo* cellTopo = grid->GetCellTopo();
	auto& cell2node = cellTopo->GetNodeIndex();
	NodeTopo* node_topo = grid->GetNodeTopo();
	auto& node_coord = node_topo->GetCoordinate();

	INTEGER4 node_num = grid->GetTotalNodeNum();
	INTEGER4 cell_num = cell2node.size();

	DArray x(node_num), y(node_num), z(node_num);

	for (int iNode = 0; iNode < node_num; ++iNode)
	{
		x[iNode] = node_coord[iNode].x();
		y[iNode] = node_coord[iNode].y();
		z[iNode] = node_coord[iNode].z();
	}
	INTEGER4 file_format = 0;
	INTEGER4 debug = 0;
	INTEGER4 vIsDouble = 1;
	INTEGER4 fileType = 0;
	string grid_name = "grid_" + grid->GetName();
	string var_name = "x y z rho u v w p";
	std::string file_name = "result/" + std::to_string(GlobalData::GetInt("currentIter")) + ".plt";
	int i = TECINI142(grid_name.c_str(),
		var_name.c_str(),
		file_name.c_str(),
		(char*)".",
		&file_format,
		&fileType,
		&debug,
		&vIsDouble);

	string zone_name = "grid_" + grid->GetName() + "_zone";
	INTEGER4 zone_type = 5;//Brick
	INTEGER4 face_num = 6;
	INTEGER4 iCellMax = 0;
	INTEGER4 jCellMax = 0;
	INTEGER4 kCellMax = 0;
	double solution_time = GlobalData::GetDouble("currentTime");
	INTEGER4 strandID = 0;
	INTEGER4 parentZn = 0;
	INTEGER4 isBlock = 1;
	INTEGER4 nFConns = 0;
	INTEGER4 FNMode = 0;
	int valueLocation[] = { 1,1,1,1,1,1,1,1 };
	int shrConn = 0;
	i = TECZNE142((char*)zone_name.c_str(),
		&zone_type,
		&node_num,
		&cell_num,
		&face_num,
		&iCellMax,
		&jCellMax,
		&kCellMax,
		&solution_time,
		&strandID,
		&parentZn,
		&isBlock,
		&nFConns,
		&FNMode,
		0,
		0,
		0,
		NULL,
		valueLocation,
		NULL,
		&shrConn);

	i = TECDAT142(&node_num, x.data(), &vIsDouble);
	i = TECDAT142(&node_num, y.data(), &vIsDouble);
	i = TECDAT142(&node_num, z.data(), &vIsDouble);
	i = TECDAT142(&node_num, density, &vIsDouble);
	i = TECDAT142(&node_num, velocity_x, &vIsDouble);
	i = TECDAT142(&node_num, velocity_y, &vIsDouble);
	i = TECDAT142(&node_num, velocity_z, &vIsDouble);
	i = TECDAT142(&node_num, pressure, &vIsDouble);
	INTEGER4 connectivityCount = cell_num * 8;
	Array<INTEGER4> cell_nodes(connectivityCount);
	for (int iCell = 0; iCell < cell_num; ++iCell)
	{
		for (int iNode = 0; iNode < cell2node[iCell].size(); ++iNode)
		{
			cell_nodes[iCell * 8 + iNode] = cell2node[iCell][iNode] + 1;
		}
	}

	i = TECNODE142(&connectivityCount, cell_nodes.data());

	/// bound face
	FaceTopo* face_topo = grid->GetFaceTopo();
	cell_num = face_topo->GetFaceNum();
	zone_name = "grid_" + grid->GetName() + "_bound";
	zone_type = 3;//Brick
	face_num = 6;
	iCellMax = 0;
	jCellMax = 0;
	kCellMax = 0;
	strandID = 0;
	parentZn = 0;
	isBlock = 1;
	nFConns = 0;
	FNMode = 0;
	i = TECZNE142((char*)zone_name.c_str(),
		&zone_type,
		&node_num,
		&cell_num,
		&face_num,
		&iCellMax,
		&jCellMax,
		&kCellMax,
		&solution_time,
		&strandID,
		&parentZn,
		&isBlock,
		&nFConns,
		&FNMode,
		0,
		0,
		0,
		NULL,
		valueLocation,
		NULL,
		&shrConn);

	i = TECDAT142(&node_num, x.data(), &vIsDouble);
	i = TECDAT142(&node_num, y.data(), &vIsDouble);
	i = TECDAT142(&node_num, z.data(), &vIsDouble);
	i = TECDAT142(&node_num, density, &vIsDouble);
	i = TECDAT142(&node_num, velocity_x, &vIsDouble);
	i = TECDAT142(&node_num, velocity_y, &vIsDouble);
	i = TECDAT142(&node_num, velocity_z, &vIsDouble);
	i = TECDAT142(&node_num, pressure, &vIsDouble);
	int node_num_per_cell = 4;
	connectivityCount = cell_num * node_num_per_cell;
	Array<INTEGER4> face_nodes(connectivityCount);
	for (int iFace = 0; iFace < cell_num; ++iFace)
	{
		int* face2node = face_topo->GetFace2Node(iFace);
		int n_node = face_topo->GetFaceNodeNum(iFace);
		for (int iNode = 0; iNode < n_node; ++iNode)
		{
			face_nodes[iFace * node_num_per_cell + iNode] = face2node[iNode] + 1;
		}
		if (n_node < 4)
		{
			for (int i = n_node; i < 4; ++i)
			{
				face_nodes[iFace * node_num_per_cell + i] = face2node[0] + 1;
			}
		}
	}

	i = TECNODE142(&connectivityCount, face_nodes.data());



	i = TECEND142();

}
void zaran::Visual::WriteTecplot2D(FieldSolver* solver)
{
	int currentIter = GlobalData::GetInt("currentIter");
	std::string filename = "result/" + std::to_string(currentIter) + ".dat";
	std::ofstream fout(filename);
	fout << "variables=x,y,rho,u,v,w,p,J,lim0,lim1,lim2,lim3,lim4,rhoGradX,rhoGradY,uGradX,uGradY,vGradX,vGradY,pGradX,pGradY\n";
	Grid* grid = solver->GetGrid();
	int nInnerNum = grid->GetInnerNodeNum();
	int nBoundNum = grid->GetBoundNodeNum();
	int nTotalNum = grid->GetTotalNodeNum();
	NodeTopo* nodeTopo = grid->GetNodeTopo();

	auto& nodeCoord = nodeTopo->GetCoordinate();
	CellTopo* cellTopo = grid->GetCellTopo();
	auto& cell2node = cellTopo->GetNodeIndex();
	auto& data = *solver->GetFieldData();
	double* rho, * u, * v, * w, * p, * jacobi, * limiterCoef0, * limiterCoef1, * limiterCoef2, * limiterCoef3, * limiterCoef4;
	double* rhoGradX, * rhoGradY, * uGradX, * uGradY, * vGradX, * vGradY, * pGradX, * pGradY;
	data.GetData("density", rho);
	data.GetData("velocity_x", u);
	data.GetData("velocity_y", v);
	data.GetData("velocity_w", w);
	data.GetData("pressure", p);
	data.GetData("coordTransJ", jacobi);
	data.GetData("limiterCoef0", limiterCoef0);
	data.GetData("limiterCoef1", limiterCoef1);
	data.GetData("limiterCoef2", limiterCoef2);
	data.GetData("limiterCoef3", limiterCoef3);
	data.GetData("limiterCoef4", limiterCoef4);
	data.GetData("rhoGradX", rhoGradX);
	data.GetData("rhoGradY", rhoGradY);
	data.GetData("uGradX", uGradX);
	data.GetData("uGradY", uGradY);
	data.GetData("vGradX", vGradX);
	data.GetData("vGradY", vGradY);
	data.GetData("pGradX", pGradX);
	data.GetData("pGradY", pGradY);
	fout << "ZONE T= grid_" << grid->GetName() << std::endl;
	fout << "N=" << grid->GetTotalNodeNum() << ", E= " << cell2node.size() << ", F=FEPOINT, ET=QUADRILATERAL" << std::endl;
	fout << "solutiontime= " << GlobalData::GetDouble("currentTime") << std::endl;
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
void zaran::Visual::WriteTecplotPoint(FieldSolver* solver)
{
	int currentIter = GlobalData::GetInt("currentIter");
	std::string filename = "result/" + std::to_string(currentIter) + ".dat";
	std::ofstream fout(filename);
	fout << "variables=x,y,z,rho,u,v,w,p\n";
	Grid* grid = solver->GetGrid();
	NodeTopo* nodeTopo = grid->GetNodeTopo();

	auto& nodeCoord = nodeTopo->GetCoordinate();
	auto& data = *solver->GetFieldData();
	double* rho, * u, * v, * w, * p;
	data.GetData("density", rho);
	data.GetData("velocity_x", u);
	data.GetData("velocity_y", v);
	data.GetData("velocity_w", w);
	data.GetData("pressure", p);
	//fout << "solutiontime= " << GlobalData::GetDouble("globalTime") << std::endl;
	for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
	{
		auto& currentCoord = nodeCoord[iNode];
		fout << currentCoord(0) << "  " << currentCoord(1) << "  " << currentCoord(2) << "  ";
		fout << rho[iNode] << "  " << u[iNode] << "  " << v[iNode] << "  " << w[iNode] << "  " << p[iNode];
		fout << std::endl;
	}
	fout.close();

	BoundaryMap* boundMap = grid->GetBoundaryMap();
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
void Visual::WriteVTK(FieldSolver* solver)
{
	// TO DO
}

void zaran::Visual::WriteTecplotZaran3D(FieldSolver* solver)
{
	int currentIter = GlobalData::GetInt("currentIter");
	std::string filename = "result/" + std::to_string(currentIter) + ".dat";
	std::ofstream fout(filename);
	fout << "variables=x,y,z,rho,u,v,w,p\n";
	auto& data = *solver->GetFieldData();
	double* rho, * u, * v, * w, * p, * jacobi;
	data.GetData("density", rho);
	data.GetData("velocity_x", u);
	data.GetData("velocity_y", v);
	data.GetData("velocity_w", w);
	data.GetData("pressure", p);
	data.GetData("coordTransJ", jacobi);
	Grid_Zaran_3D* grid = static_cast<Grid_Zaran_3D*>(solver->GetGrid());
	CellTopoZaran* cellTopo = grid->GetCellTopo();
	Array<CellType>& cell_type = cellTopo->GetType();
	auto& cell_center = cellTopo->GetCenterCoord();
	NodeTopo* node_topo = grid->GetNodeTopo();
	auto& node_coord = node_topo->GetCoordinate();
	int is, ie, js, je, ks, ke;
	grid->GetRange(is, ie, js, je, ks, ke);
	fout << "ZONE T= grid_" << grid->GetName() << std::endl;
	fout << "I=" << ie - is << ", J=" << je - js << ", K=" << ke - ks << ", DATAPACKING=BLOCK, VARLOCATION=([4-8]=CELLCENTERED)" << std::endl;
	fout << "solutiontime= " << GlobalData::GetDouble("currentTime") << std::endl;
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

void zaran::Visual::WriteTecplotZaran3DBinary(FieldSolver* solver)
{
	auto& data = *solver->GetFieldData();
	double* rho, * u, * v, * w, * p, * jacobi;
	data.GetData("density", rho);
	data.GetData("velocity_x", u);
	data.GetData("velocity_y", v);
	data.GetData("velocity_w", w);
	data.GetData("pressure", p);
	data.GetData("coordTransJ", jacobi);

	Grid_Zaran_3D* grid = static_cast<Grid_Zaran_3D*>(solver->GetGrid());
	CellTopoZaran* cellTopo = grid->GetCellTopo();
	Array<CellType>& cell_type = cellTopo->GetType();
	auto& cell_center = cellTopo->GetCenterCoord();
	NodeTopo* node_topo = grid->GetNodeTopo();
	auto& node_coord = node_topo->GetCoordinate();

	int is, ie, js, je, ks, ke;
	grid->GetRange(is, ie, js, je, ks, ke);
	int n_node = (ie - is) * (je - js) * (ke - ks);
	int n_cell = (ie - is - 1) * (je - js - 1) * (ke - ks - 1);
	DArray x(n_node), y(n_node), z(n_node), density(n_cell), x_vel(n_cell), y_vel(n_cell), z_vel(n_cell), pressure(n_cell);
	IArray cell_type_array(n_cell);
	int currentIter = GlobalData::GetInt("currentIter");
	std::string file_name = "result/" + std::to_string(currentIter) + ".plt";
	string zone_name = "grid_" + grid->GetName();
	string var_name = "x y z rho u v w p cell_type";
	double solution_time = GlobalData::GetDouble("currentTime");
	int node_index = 0;
	for (int k = ks; k < ke; ++k)
	{
		for (int j = js; j < je; ++j)
		{
			for (int i = is; i < ie; ++i)
			{
				node_index = (k - ks) * (je - js) * (ie - is) + (j - js) * (ie - is) + (i - is);
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
				cell_index = (k - ks) * (je - js - 1) * (ie - is - 1) + (j - js) * (ie - is - 1) + (i - is);
				density[cell_index] = rho[grid->GetCellIndex(i, j, k)];
				x_vel[cell_index] = u[grid->GetCellIndex(i, j, k)];
				y_vel[cell_index] = v[grid->GetCellIndex(i, j, k)];
				z_vel[cell_index] = w[grid->GetCellIndex(i, j, k)];
				pressure[cell_index] = p[grid->GetCellIndex(i, j, k)];
				cell_type_array[cell_index] = (int)cell_type[grid->GetCellIndex(i, j, k)];
			}
		}
	}

	int value_location[9] = { 1,1,1,0,0,0,0,0,0 };
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
	VIsDouble = 0;
	i = TECDAT142(&n_cell, cell_type_array.data(), &VIsDouble);
	i = TECEND142();
}

