#include "Create_Zaran_3D.h"
#include "Grid_Zaran_3D.h"
#include"ReadSTL.h"
#include"PolyData.h"
#include "Log.h"
#include"CellTopoInfoZaran.h"
#include<fstream>
void zaran::GridListFactoryZaran3D::Create(Ptr<GridList>& gridList)
{
	if (gridList.get() == nullptr)
	{
		gridList = std::make_shared<GridList>();
	}
	CreateStructPart(gridList);
	TagCell(gridList);
	CrateBoundPatch(gridList);
}

void zaran::GridListFactoryZaran3D::CreateStructPart(Ptr<GridList>& gridList)
{
	if (gridList.get() == nullptr)
	{
		gridList = std::make_shared<GridList>();
	}

	Ptr<Grid_Zaran_3D>grid = std::make_shared<Grid_Zaran_3D>();
	grid->SetDimension(Dimension::three);
	grid->SetIndex(0);
	grid->SetLevel(0);
	grid->SetName("noname");
	grid->SetType(GridType::Zaran_3D);
	auto& nodeTopo = grid->GetNodeTopo();
	auto& nodeCoord = nodeTopo->GetCoordinate();
	int xNodeNum = GlobalData::GetInt("nx");
	int yNodeNum = GlobalData::GetInt("ny");
	int zNodeNum = GlobalData::GetInt("nz");
	xNodeNum++;
	yNodeNum++;
	zNodeNum++;
	int nodeNum = 0;
	double xMin = GlobalData::GetDouble("xMin");
	double yMin = GlobalData::GetDouble("yMin");
	double zMin = GlobalData::GetDouble("zMin");
	double xMax = GlobalData::GetDouble("xMax");
	double yMax = GlobalData::GetDouble("yMax");
	double zMax = GlobalData::GetDouble("zMax");
	grid->SetBox(xMin, xMax, yMin, yMax, zMin, zMax);
	double dx = (xMax - xMin) / (xNodeNum - 1);
	double dy = (yMax - yMin) / (yNodeNum - 1);
	double dz = (zMax - zMin) / (zNodeNum - 1);
	int i, j, k;
	double x, y, z;
	grid->SetNodeNum((xNodeNum + 2), (yNodeNum + 2), (zNodeNum + 2));
	nodeCoord.resize((xNodeNum + 2) * (yNodeNum + 2) * (zNodeNum + 2));

	auto& nodeType = nodeTopo->GetType();
	nodeType.resize((xNodeNum + 2) * (yNodeNum + 2) * (zNodeNum + 2));
	auto& tempI = nodeTopo->GetTemplateI();
	auto& tempJ = nodeTopo->GetTemplateJ();
	auto& tempK = nodeTopo->GetTemplateK();
	tempI.resize((xNodeNum + 2) * (yNodeNum + 2) * (zNodeNum + 2));
	tempJ.resize((xNodeNum + 2) * (yNodeNum + 2) * (zNodeNum + 2));
	tempK.resize((xNodeNum + 2) * (yNodeNum + 2) * (zNodeNum + 2));
	int iNode = 0;
	for (k = 0; k < zNodeNum + 2; ++k)
	{
		z = k * dz + zMin - dz;
		for (j = 0; j < yNodeNum + 2; ++j)
		{
			y = j * dy + yMin - dy;
			for (i = 0; i < xNodeNum + 2; ++i)
			{
				x = i * dx + xMin - dx;
				iNode = grid->GetNodeIndex(i, j, k);
				nodeCoord[iNode] = { x,y,z };
				if (i == 0 || j == 0 || k == 0 || i == xNodeNum + 1 || j == yNodeNum + 1 || k == zNodeNum + 1)
				{
					nodeType[iNode] = NodeType::ghost;
					continue;
				}
				if (i == 1 || j == 1 || k == 1 || i == xNodeNum || j == yNodeNum || k == zNodeNum)
				{
					nodeType[iNode] = NodeType::inlet;
				}
				tempI[iNode] = { grid->GetNodeIndex(i - 1, j, k), grid->GetNodeIndex(i, j, k), grid->GetNodeIndex(i + 1, j, k) };
				tempJ[iNode] = { grid->GetNodeIndex(i, j - 1, k), grid->GetNodeIndex(i, j, k), grid->GetNodeIndex(i, j + 1, k) };
				tempK[iNode] = { grid->GetNodeIndex(i, j, k - 1), grid->GetNodeIndex(i, j, k), grid->GetNodeIndex(i, j, k + 1) };
			}
		}
	}
	auto& nodeNeighbor = nodeTopo->GetNeighborCloud();
	nodeNeighbor.resize((xNodeNum + 2) * (yNodeNum + 2) * (zNodeNum + 2));
	for (k = 0; k < zNodeNum + 2; ++k)
	{
		for (j = 0; j < yNodeNum + 2; ++j)
		{
			for (i = 0; i < xNodeNum + 2; ++i)
			{
				iNode = grid->GetNodeIndex(i, j, k);
				if (i == 0)
				{
					if (j == 0)
					{
						if (k == 0)
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i + 1, j, k),  grid->GetNodeIndex(i, j + 1, k),  grid->GetNodeIndex(i, j, k + 1) };
						}
						else if (k == zNodeNum + 1)
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i + 1, j, k),  grid->GetNodeIndex(i, j + 1, k), grid->GetNodeIndex(i, j, k - 1) };
						}
						else
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i + 1, j, k),  grid->GetNodeIndex(i, j + 1, k), grid->GetNodeIndex(i, j, k - 1), grid->GetNodeIndex(i, j, k + 1) };
						}
					}
					else if (j == yNodeNum + 1)
					{
						if (k == 0)
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i + 1, j, k), grid->GetNodeIndex(i, j - 1, k), grid->GetNodeIndex(i, j, k + 1) };

						}
						else if (k == zNodeNum + 1)
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i + 1, j, k), grid->GetNodeIndex(i, j - 1, k),  grid->GetNodeIndex(i, j, k - 1) };


						}
						else
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i + 1, j, k), grid->GetNodeIndex(i, j - 1, k),  grid->GetNodeIndex(i, j, k - 1), grid->GetNodeIndex(i, j, k + 1) };

						}
					}
					else
					{
						if (k == 0)
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i + 1, j, k), grid->GetNodeIndex(i, j - 1, k), grid->GetNodeIndex(i, j + 1, k),  grid->GetNodeIndex(i, j, k + 1) };

						}
						else if (k == zNodeNum + 1)
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i + 1, j, k), grid->GetNodeIndex(i, j - 1, k), grid->GetNodeIndex(i, j + 1, k), grid->GetNodeIndex(i, j, k - 1) };

						}
						else
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i + 1, j, k), grid->GetNodeIndex(i, j - 1, k), grid->GetNodeIndex(i, j + 1, k), grid->GetNodeIndex(i, j, k - 1), grid->GetNodeIndex(i, j, k + 1) };
						}
					}
				}
				else if (i == xNodeNum + 1)
				{
					if (j == 0)
					{
						if (k == 0)
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i - 1, j, k),  grid->GetNodeIndex(i, j + 1, k),  grid->GetNodeIndex(i, j, k + 1) };

						}
						else if (k == zNodeNum + 1)
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i - 1, j, k),  grid->GetNodeIndex(i, j + 1, k), grid->GetNodeIndex(i, j, k - 1) };
						}
						else
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i - 1, j, k),  grid->GetNodeIndex(i, j + 1, k), grid->GetNodeIndex(i, j, k - 1), grid->GetNodeIndex(i, j, k + 1) };
						}
					}
					else if (j == yNodeNum + 1)
					{
						if (k == 0)
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i - 1, j, k), grid->GetNodeIndex(i, j - 1, k), grid->GetNodeIndex(i, j, k + 1) };


						}
						else if (k == zNodeNum + 1)
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i - 1, j, k), grid->GetNodeIndex(i, j - 1, k),  grid->GetNodeIndex(i, j, k - 1) };

						}
						else
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i - 1, j, k), grid->GetNodeIndex(i, j - 1, k),  grid->GetNodeIndex(i, j, k - 1), grid->GetNodeIndex(i, j, k + 1) };

						}
					}
					else
					{
						if (k == 0)
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i - 1, j, k), grid->GetNodeIndex(i, j - 1, k), grid->GetNodeIndex(i, j + 1, k),  grid->GetNodeIndex(i, j, k + 1) };

						}
						else if (k == zNodeNum + 1)
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i - 1, j, k), grid->GetNodeIndex(i, j - 1, k), grid->GetNodeIndex(i, j + 1, k), grid->GetNodeIndex(i, j, k - 1) };

						}
						else
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i - 1, j, k), grid->GetNodeIndex(i, j - 1, k), grid->GetNodeIndex(i, j + 1, k), grid->GetNodeIndex(i, j, k - 1), grid->GetNodeIndex(i, j, k + 1) };

						}
					}
				}
				else
				{
					if (j == 0)
					{
						if (k == 0)
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i - 1, j, k), grid->GetNodeIndex(i + 1, j, k), grid->GetNodeIndex(i, j + 1, k),  grid->GetNodeIndex(i, j, k + 1) };

						}
						else if (k == zNodeNum + 1)
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i - 1, j, k), grid->GetNodeIndex(i + 1, j, k), grid->GetNodeIndex(i, j + 1, k), grid->GetNodeIndex(i, j, k - 1) };

						}
						else
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i - 1, j, k), grid->GetNodeIndex(i + 1, j, k),  grid->GetNodeIndex(i, j + 1, k), grid->GetNodeIndex(i, j, k - 1), grid->GetNodeIndex(i, j, k + 1) };

						}
					}
					else if (j == yNodeNum + 1)
					{
						if (k == 0)
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i - 1, j, k), grid->GetNodeIndex(i + 1, j, k), grid->GetNodeIndex(i, j - 1, k),   grid->GetNodeIndex(i, j, k + 1) };

						}
						else if (k == zNodeNum + 1)
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i - 1, j, k), grid->GetNodeIndex(i + 1, j, k), grid->GetNodeIndex(i, j - 1, k),  grid->GetNodeIndex(i, j, k - 1) };

						}
						else
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i - 1, j, k), grid->GetNodeIndex(i + 1, j, k), grid->GetNodeIndex(i, j - 1, k),  grid->GetNodeIndex(i, j, k - 1), grid->GetNodeIndex(i, j, k + 1) };

						}
					}
					else
					{
						if (k == 0)
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i - 1, j, k), grid->GetNodeIndex(i + 1, j, k), grid->GetNodeIndex(i, j - 1, k), grid->GetNodeIndex(i, j + 1, k),  grid->GetNodeIndex(i, j, k + 1) };
						}
						else if (k == zNodeNum + 1)
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i - 1, j, k), grid->GetNodeIndex(i + 1, j, k), grid->GetNodeIndex(i, j - 1, k), grid->GetNodeIndex(i, j + 1, k), grid->GetNodeIndex(i, j, k - 1) };

						}
						else
						{
							nodeNeighbor[iNode] = { grid->GetNodeIndex(i - 1, j, k), grid->GetNodeIndex(i + 1, j, k), grid->GetNodeIndex(i, j - 1, k), grid->GetNodeIndex(i, j + 1, k), grid->GetNodeIndex(i, j, k - 1), grid->GetNodeIndex(i, j, k + 1) };
						}
					}
				}
			}
		}
	}
	grid->SetInnerNodeNum((xNodeNum - 2) * (yNodeNum - 2) * (zNodeNum - 2));
	grid->SetTotalNodeNum((xNodeNum + 2) * (yNodeNum + 2) * (zNodeNum + 2));
	grid->SetBoundNodeNum(xNodeNum * yNodeNum * zNodeNum - grid->GetInnerNodeNum());
	grid->SetInterNodeInfo(std::make_shared<InterNodeInfo>());
	grid->SetBoundaryMap(std::make_shared<BoundaryMap>());

	auto& cellTopo = grid->GetCellTopo();
	auto& cell2node = cellTopo->GetNodeIndex();
	cell2node.resize((xNodeNum + 1) * (yNodeNum + 1) * (zNodeNum + 1));
	int iCell;
	for (k = 0; k < zNodeNum + 1; ++k)
	{
		for (j = 0; j < yNodeNum + 1; ++j)
		{

			for (i = 0; i < xNodeNum + 1; ++i)
			{
				iCell = grid->GetCellIndex(i, j, k);
				cell2node[iCell] = (IArray{ grid->GetNodeIndex(i, j, k), grid->GetNodeIndex(i + 1, j, k), grid->GetNodeIndex(i + 1, j + 1, k), grid->GetNodeIndex(i, j + 1, k), grid->GetNodeIndex(i, j, k + 1), grid->GetNodeIndex(i + 1, j, k + 1), grid->GetNodeIndex(i + 1, j + 1, k + 1), grid->GetNodeIndex(i, j + 1, k + 1) });
			}
		}
	}

	auto boundMap = grid->GetBoundaryMap();
	int nodeIndex, innerNodeIndex, ghostNodeIndex;
	DVector3D boundNorm;
	Boundary bound;
	//i方向两个面分别为入口和出口
	for (k = 1; k < zNodeNum + 1; ++k)
	{

		for (j = 1; j < yNodeNum + 1; ++j)
		{
			i = 1;
			nodeIndex = grid->GetNodeIndex(i, j, k);
			bound.SetNodeIndex(nodeIndex);
			innerNodeIndex = grid->GetNodeIndex(i + 1, j, k);
			bound.SetInnerNodeIndex(innerNodeIndex);
			ghostNodeIndex = grid->GetNodeIndex(i - 1, j, k);
			bound.SetGhostNodeIndex(ghostNodeIndex);
			boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
			bound.SetNorm(boundNorm);
			boundMap->AddBoundary("inlet", bound);

			i = xNodeNum;
			nodeIndex = grid->GetNodeIndex(i, j, k);
			bound.SetNodeIndex(nodeIndex);
			innerNodeIndex = grid->GetNodeIndex(i - 1, j, k);
			bound.SetInnerNodeIndex(innerNodeIndex);
			ghostNodeIndex = grid->GetNodeIndex(i + 1, j, k);
			bound.SetGhostNodeIndex(ghostNodeIndex);
			boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
			bound.SetNorm(boundNorm);
			boundMap->AddBoundary("outlet", bound);
		}
	}

	//其余四个面均为壁面
	for (i = 1; i < xNodeNum + 1; ++i)
	{
		//j方向
		for (k = 1; k < zNodeNum + 1; ++k)
		{
			j = 1;
			nodeIndex = grid->GetNodeIndex(i, j, k);
			bound.SetNodeIndex(nodeIndex);
			innerNodeIndex = grid->GetNodeIndex(i, j + 1, k);
			bound.SetInnerNodeIndex(innerNodeIndex);
			ghostNodeIndex = grid->GetNodeIndex(i, j - 1, k);
			bound.SetGhostNodeIndex(ghostNodeIndex);
			boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
			bound.SetNorm(boundNorm);
			boundMap->AddBoundary("wall", bound);
			j = yNodeNum;
			nodeIndex = grid->GetNodeIndex(i, j, k);
			bound.SetNodeIndex(nodeIndex);
			innerNodeIndex = grid->GetNodeIndex(i, j - 1, k);
			bound.SetInnerNodeIndex(innerNodeIndex);
			ghostNodeIndex = grid->GetNodeIndex(i, j + 1, k);
			bound.SetGhostNodeIndex(ghostNodeIndex);
			boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
			bound.SetNorm(boundNorm);
			boundMap->AddBoundary("wall", bound);
		}

		//j方向
		for (j = 1; j < yNodeNum + 1; ++j)
		{
			k = 1;
			nodeIndex = grid->GetNodeIndex(i, j, k);
			bound.SetNodeIndex(nodeIndex);
			innerNodeIndex = grid->GetNodeIndex(i, j, k + 1);
			bound.SetInnerNodeIndex(innerNodeIndex);
			ghostNodeIndex = grid->GetNodeIndex(i, j, k - 1);
			bound.SetGhostNodeIndex(ghostNodeIndex);
			boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
			bound.SetNorm(boundNorm);
			boundMap->AddBoundary("wall", bound);
			k = zNodeNum;
			nodeIndex = grid->GetNodeIndex(i, j, k);
			bound.SetNodeIndex(nodeIndex);
			innerNodeIndex = grid->GetNodeIndex(i, j, k - 1);
			bound.SetInnerNodeIndex(innerNodeIndex);
			ghostNodeIndex = grid->GetNodeIndex(i, j, k + 1);
			bound.SetGhostNodeIndex(ghostNodeIndex);
			boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
			bound.SetNorm(boundNorm);
			boundMap->AddBoundary("wall", bound);
		}
	}
	gridList->AddGrid(grid);
}

void zaran::GridListFactoryZaran3D::TagCell(Ptr<GridList>& gridList)
{
	string modelFileName = GlobalData::GetString("modelFileName");
	STLReader reader;
	reader.ReadSTLFile(modelFileName.c_str());
	PolyDataModel model;
	model.SetPolyData(reader.GetMesh(), 1e-6);
	if (model.IsClosed())
		ZaranLog::info("Import Model: {}, is closed!", modelFileName);
	else
		ZaranLog::info("Import Model: {}, is not closed!", modelFileName);
	Ptr<zaran::Grid_Zaran_3D> grid = std::static_pointer_cast<Grid_Zaran_3D>(gridList->GetGrid(0));
	auto cellTopo = std::static_pointer_cast<CellTopoInfoZaran>(grid->GetCellTopo());
	auto& cell_type = cellTopo->GetType();
	int ni, nj, nk;
	grid->GetNodeNum(ni, nj, nk);
	cell_type.resize((ni - 1) * (nj - 1) * (nk - 1));
	int i, j, k;
	double x_min, x_max, y_min, y_max, z_min, z_max;
	grid->GetBox(x_min, x_max, y_min, y_max, z_min, z_max);
	double dx = (x_max - x_min) / (ni - 1);
	double dy = (y_max - y_min) / (nj - 1);
	double dz = (z_max - z_min) / (nk - 1);
	DVector3D cell_center;
	for (int iCell = 0;iCell < cell_type.size();++iCell)
	{
		grid->GetCellIndex(iCell, i, j, k);
		if (cell_type[iCell] == CellType::Solid)
		{
			std::cout << "debug" << std::endl;
		}
		cell_center[0] = (i + 0.5) * dx + x_min;

		cell_center[1] = (j + 0.5) * dy + y_min;
		cell_center[2] = (k + 0.5) * dz + z_min;
		if (model.InModel(cell_center))
		{
			cell_type[iCell] = CellType::Solid;
		}
		else
		{
			cell_type[iCell] = CellType::Fluid;
		}
	}
	for (int iCell = 0;iCell < cell_type.size();++iCell)
	{
		grid->GetCellIndex(iCell, i, j, k);
		if (cell_type[iCell] == CellType::Fluid)
		{
			if (i > 0 && cell_type[grid->GetCellIndex(i - 1, j, k)] == CellType::Solid)
			{
				cell_type[iCell] = CellType::FluidSolid;
				continue;
			}
			if (i < ni && cell_type[grid->GetCellIndex(i + 1, j, k)] == CellType::Solid)
			{
				cell_type[iCell] = CellType::FluidSolid;
				continue;
			}
			if (j > 0 && cell_type[grid->GetCellIndex(i, j - 1, k)] == CellType::Solid)
			{
				cell_type[iCell] = CellType::FluidSolid;
				continue;
			}
			if (j < nj && cell_type[grid->GetCellIndex(i, j + 1, k)] == CellType::Solid)
			{
				cell_type[iCell] = CellType::FluidSolid;
				continue;
			}
			if (k > 0 && cell_type[grid->GetCellIndex(i, j, k - 1)] == CellType::Solid)
			{
				cell_type[iCell] = CellType::FluidSolid;
				continue;
			}
			if (k < nk && cell_type[grid->GetCellIndex(i, j, k + 1)] == CellType::Solid)
			{
				cell_type[iCell] = CellType::FluidSolid;
				continue;
			}
		}
	}
	std::ofstream fout("cell_type.dat");
	for (int iCell = 0;iCell < cell_type.size();++iCell)
	{
		grid->GetCellIndex(iCell, i, j, k);
		cell_center[0] = (i + 0.5) * dx + x_min;
		cell_center[1] = (j + 0.5) * dy + y_min;
		cell_center[2] = (k + 0.5) * dz + z_min;
		if (i == 25 || j == 25 || k == 25)
			fout << cell_center[0] << " " << cell_center[1] << " " << cell_center[2] << " " << int(cell_type[iCell]) << std::endl;
	}

}

void zaran::GridListFactoryZaran3D::CrateBoundPatch(Ptr<GridList>& gridList)
{
}
