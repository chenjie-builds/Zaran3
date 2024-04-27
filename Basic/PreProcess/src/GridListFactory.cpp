#include "GridListFactory.h"
#include"Grid_Struct_3D.h"
#include"Grid_Struct_2D.h"
#include "grid.h"
#include "log.h"
namespace zaran
{

	void GridListFactory::Create(Grid*& grid)
	{
		std::string createMethod = GlobalData::GetString("createGridMethod");
		if (createMethod == "Plot3D")
		{
			//CreateByTest(gridList);
			//CreateStructGrid3D(gridList);
			//CreateGridFNFDM2D(gridList);
			CreateStructGrid2D(grid);
		}
		else
		{
			Log::warn("Unsupported Grid Generate Method:{}, Please Check!", createMethod);
			system("pause");
		}
	}

	void GridListFactory::ReadPlot3D(Grid*& grid)
	{
		if (grid != nullptr)
		{
			delete[] grid;
		}
		grid = new Grid();
		grid->SetDimension(Dimension::three);
		grid->SetIndex(0);
		grid->SetLevel(0);
		grid->SetName("noname");
		grid->SetType(GridType::Unkown);
	}

	void GridListFactory::CreateByTest(Grid*& grid)
	{
		if (grid != nullptr)
		{
			delete[] grid;
		}
		grid = new Grid();
		grid->SetDimension(Dimension::three);
		grid->SetIndex(0);
		grid->SetLevel(0);
		grid->SetName("noname");
		grid->SetType(GridType::Unkown);
		NodeTopo* nodeTopo = grid->GetNodeTopo();

		auto& nodeCoord = nodeTopo->GetCoordinate();
		int xNodeNum = 31;
		int yNodeNum = 31;
		int zNodeNum = 31;
		int nodeNum = 0;
		double xMin = 0.0;
		double yMin = 0.0;
		double zMin = 0.0;
		double xMax = 1.0;
		double yMax = 1.0;
		double zMax = 1.0;
		double dx = (xMax - xMin) / (xNodeNum - 1);
		double dy = (yMax - yMin) / (yNodeNum - 1);
		double dz = (zMax - zMin) / (zNodeNum - 1);
		int i, j, k;
		double x, y, z;
		nodeCoord.resize((xNodeNum + 2) * (yNodeNum + 2) * (zNodeNum + 2));
		//生成结构网格ijk索引对应的节点编号
		Array<Array<IArray>> structNodeIndex;
		structNodeIndex.resize(xNodeNum + 2);
		for (i = 0; i < xNodeNum + 2; ++i)
		{
			structNodeIndex[i].resize(yNodeNum + 2);
			for (j = 0; j < yNodeNum + 2; ++j)
			{
				structNodeIndex[i][j].resize(zNodeNum + 2);
			}
		}
		//内部节点
		for (k = 2; k < zNodeNum; ++k)
		{
			for (j = 2; j < yNodeNum; ++j)
			{
				for (i = 2; i < xNodeNum; ++i)
				{
					structNodeIndex[i][j][k] = nodeNum++;
				}
			}
		}
		//边界节点
		for (k = 0; k < zNodeNum + 2; ++k)
		{
			for (j = 0; j < yNodeNum + 2; ++j)
			{

				for (i = 0; i < xNodeNum + 2; ++i)
				{
					if (k != 1 && k != zNodeNum && j != 1 && j != yNodeNum && i != 1 && i != xNodeNum)
					{
						continue;
					}
					if (k == 0 || k == zNodeNum + 1 || j == 0 || j == yNodeNum + 1 || i == 0 || i == xNodeNum + 1)
					{
						continue;
					}
					structNodeIndex[i][j][k] = nodeNum++;
				}
			}
		}
		//镜像节点
		for (k = 0; k < zNodeNum + 2; ++k)
		{
			for (j = 0; j < yNodeNum + 2; ++j)
			{

				for (i = 0; i < xNodeNum + 2; ++i)
				{
					if (k != 0 && k != zNodeNum + 1 && j != 0 && j != yNodeNum + 1 && i != 0 && i != xNodeNum + 1)
					{
						continue;
					}
					structNodeIndex[i][j][k] = nodeNum++;
				}
			}
		}
		auto& nodeType = nodeTopo->GetType();
		nodeType.resize((xNodeNum + 2) * (yNodeNum + 2) * (zNodeNum + 2));
		auto& tempI = nodeTopo->GetTemplateI();
		auto& tempJ = nodeTopo->GetTemplateJ();
		auto& tempK = nodeTopo->GetTemplateK();
		tempI.resize((xNodeNum + 2) * (yNodeNum + 2) * (zNodeNum + 2));
		tempJ.resize((xNodeNum + 2) * (yNodeNum + 2) * (zNodeNum + 2));
		tempK.resize((xNodeNum + 2) * (yNodeNum + 2) * (zNodeNum + 2));
		for (k = 0; k < zNodeNum + 2; ++k)
		{
			z = k * dz + zMin - dz;
			for (j = 0; j < yNodeNum + 2; ++j)
			{
				y = j * dy + yMin - dy;
				for (i = 0; i < xNodeNum + 2; ++i)
				{
					int iNode = structNodeIndex[i][j][k];
					x = i * dx + xMin - dx;
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
					tempI[iNode] = { structNodeIndex[i - 1][j][k], structNodeIndex[i][j][k], structNodeIndex[i + 1][j][k] };
					tempI[iNode] = { structNodeIndex[i][j - 1][k], structNodeIndex[i][j][k], structNodeIndex[i][j + 1][k] };
					tempI[iNode] = { structNodeIndex[i][j][k - 1], structNodeIndex[i][j][k], structNodeIndex[i][j][k + 1] };
				}
			}
		}
		auto& nodeNeighbor = nodeTopo->GetNeighborCloud();
		for (k = 0; k < zNodeNum + 2; ++k)
		{
			for (j = 0; j < yNodeNum + 2; ++j)
			{
				for (i = 0; i < xNodeNum + 2; ++i)
				{
					int iNode = structNodeIndex[i][j][k];
					if (i == 0)
					{
						if (j == 0)
						{
							if (k == 0)
							{
								nodeNeighbor[iNode] = { structNodeIndex[i + 1][j][k],  structNodeIndex[i][j + 1][k],  structNodeIndex[i][j][k + 1] };
							}
							else if (k == zNodeNum + 1)
							{
								nodeNeighbor[iNode] = { structNodeIndex[i + 1][j][k],  structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k - 1] };
							}
							else
							{
								nodeNeighbor[iNode] = { structNodeIndex[i + 1][j][k],  structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k - 1], structNodeIndex[i][j][k + 1] };
							}
						}
						else if (j == yNodeNum + 1)
						{
							if (k == 0)
							{
								nodeNeighbor[iNode] = { structNodeIndex[i + 1][j][k], structNodeIndex[i][j - 1][k], structNodeIndex[i][j][k + 1] };

							}
							else if (k == zNodeNum + 1)
							{
								nodeNeighbor[iNode] = { structNodeIndex[i + 1][j][k], structNodeIndex[i][j - 1][k],  structNodeIndex[i][j][k - 1] };

							}
							else
							{
								nodeNeighbor[iNode] = { structNodeIndex[i + 1][j][k], structNodeIndex[i][j - 1][k], structNodeIndex[i][j][k - 1], structNodeIndex[i][j][k + 1] };

							}
						}
						else
						{
							if (k == 0)
							{
								nodeNeighbor[iNode] = { structNodeIndex[i + 1][j][k], structNodeIndex[i][j - 1][k], structNodeIndex[i][j + 1][k],  structNodeIndex[i][j][k + 1] };

							}
							else if (k == zNodeNum + 1)
							{
								nodeNeighbor[iNode] = { structNodeIndex[i + 1][j][k], structNodeIndex[i][j - 1][k], structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k - 1] };

							}
							else
							{
								nodeNeighbor[iNode] = { structNodeIndex[i + 1][j][k], structNodeIndex[i][j - 1][k], structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k - 1], structNodeIndex[i][j][k + 1] };

							}
						}
					}
					else if (i == xNodeNum + 1)
					{
						if (j == 0)
						{
							if (k == 0)
							{
								nodeNeighbor[iNode] = { structNodeIndex[i - 1][j][k],  structNodeIndex[i][j + 1][k],  structNodeIndex[i][j][k + 1] };
							}
							else if (k == zNodeNum + 1)
							{
								nodeNeighbor[iNode] = { structNodeIndex[i - 1][j][k],  structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k - 1] };
							}
							else
							{
								nodeNeighbor[iNode] = { structNodeIndex[i - 1][j][k],  structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k - 1], structNodeIndex[i][j][k + 1] };
							}
						}
						else if (j == yNodeNum + 1)
						{
							if (k == 0)
							{
								nodeNeighbor[iNode] = { structNodeIndex[i - 1][j][k], structNodeIndex[i][j - 1][k],  structNodeIndex[i][j][k + 1] };

							}
							else if (k == zNodeNum + 1)
							{
								nodeNeighbor[iNode] = { structNodeIndex[i - 1][j][k], structNodeIndex[i][j - 1][k],  structNodeIndex[i][j][k - 1] };

							}
							else
							{
								nodeNeighbor[iNode] = { structNodeIndex[i - 1][j][k], structNodeIndex[i][j - 1][k],  structNodeIndex[i][j][k - 1], structNodeIndex[i][j][k + 1] };

							}
						}
						else
						{
							if (k == 0)
							{
								nodeNeighbor[iNode] = { structNodeIndex[i - 1][j][k], structNodeIndex[i][j - 1][k], structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k + 1] };

							}
							else if (k == zNodeNum + 1)
							{
								nodeNeighbor[iNode] = { structNodeIndex[i - 1][j][k], structNodeIndex[i][j - 1][k], structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k - 1] };

							}
							else
							{
								nodeNeighbor[iNode] = { structNodeIndex[i - 1][j][k], structNodeIndex[i][j - 1][k], structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k - 1], structNodeIndex[i][j][k + 1] };

							}
						}
					}
					else
					{
						if (j == 0)
						{
							if (k == 0)
							{
								nodeNeighbor[iNode] = { structNodeIndex[i - 1][j][k], structNodeIndex[i + 1][j][k], structNodeIndex[i][j + 1][k],structNodeIndex[i][j][k + 1] };

							}
							else if (k == zNodeNum + 1)
							{
								nodeNeighbor[iNode] = { structNodeIndex[i - 1][j][k], structNodeIndex[i + 1][j][k], structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k - 1] };

							}
							else
							{
								nodeNeighbor[iNode] = { structNodeIndex[i - 1][j][k], structNodeIndex[i + 1][j][k],  structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k - 1], structNodeIndex[i][j][k + 1] };

							}
						}
						else if (j == yNodeNum + 1)
						{
							if (k == 0)
							{
								nodeNeighbor[iNode] = { structNodeIndex[i - 1][j][k], structNodeIndex[i + 1][j][k], structNodeIndex[i][j - 1][k],   structNodeIndex[i][j][k + 1] };

							}
							else if (k == zNodeNum + 1)
							{
								nodeNeighbor[iNode] = { structNodeIndex[i - 1][j][k], structNodeIndex[i + 1][j][k], structNodeIndex[i][j - 1][k],  structNodeIndex[i][j][k - 1] };

							}
							else
							{
								nodeNeighbor[iNode] = { structNodeIndex[i - 1][j][k], structNodeIndex[i + 1][j][k], structNodeIndex[i][j - 1][k],  structNodeIndex[i][j][k - 1], structNodeIndex[i][j][k + 1] };

							}
						}
						else
						{
							if (k == 0)
							{
								nodeNeighbor[iNode] = { structNodeIndex[i - 1][j][k], structNodeIndex[i + 1][j][k], structNodeIndex[i][j - 1][k], structNodeIndex[i][j + 1][k],  structNodeIndex[i][j][k + 1] };

							}
							else if (k == zNodeNum + 1)
							{
								nodeNeighbor[iNode] = { structNodeIndex[i - 1][j][k], structNodeIndex[i + 1][j][k], structNodeIndex[i][j - 1][k], structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k - 1] };

							}
							else
							{
								nodeNeighbor[iNode] = { structNodeIndex[i - 1][j][k], structNodeIndex[i + 1][j][k], structNodeIndex[i][j - 1][k], structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k - 1], structNodeIndex[i][j][k + 1] };

							}
						}
					}
				}
			}
		}
		grid->SetInnerNodeNum((xNodeNum - 2) * (yNodeNum - 2) * (zNodeNum - 2));
		grid->SetTotalNodeNum((xNodeNum + 2) * (yNodeNum + 2) * (zNodeNum + 2));
		grid->SetBoundNodeNum(xNodeNum * yNodeNum * zNodeNum - grid->GetInnerNodeNum());
		CellTopo* cellTopo = grid->GetCellTopo();
		auto& cell2node = cellTopo->GetNodeIndex();
		cell2node.resize((xNodeNum - 1) * (yNodeNum - 1) * (zNodeNum - 1));
		auto& iterCell = cell2node.begin();
		for (k = 1; k < zNodeNum; ++k)
		{
			for (j = 1; j < yNodeNum; ++j)
			{

				for (i = 1; i < xNodeNum; ++i)
				{
					*iterCell = (IArray{ structNodeIndex[i][j][k], structNodeIndex[i + 1][j][k], structNodeIndex[i + 1][j + 1][k], structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k + 1], structNodeIndex[i + 1][j][k + 1], structNodeIndex[i + 1][j + 1][k + 1], structNodeIndex[i][j + 1][k + 1] });
					iterCell++;
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
				nodeIndex = structNodeIndex[i][j][k];
				bound.SetNodeIndex(nodeIndex);
				innerNodeIndex = structNodeIndex[i + 1][j][k];
				bound.SetInnerIndex(innerNodeIndex);
				ghostNodeIndex = structNodeIndex[i - 1][j][k];
				bound.SetGhostIndex(ghostNodeIndex);
				boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
				bound.SetNorm(boundNorm);
				boundMap->AddBoundary("inlet", bound);

				i = xNodeNum;
				nodeIndex = structNodeIndex[i][j][k];
				bound.SetNodeIndex(nodeIndex);
				innerNodeIndex = structNodeIndex[i - 1][j][k];
				bound.SetInnerIndex(innerNodeIndex);
				ghostNodeIndex = structNodeIndex[i + 1][j][k];
				bound.SetGhostIndex(ghostNodeIndex);
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
				nodeIndex = structNodeIndex[i][j][k];
				bound.SetNodeIndex(nodeIndex);
				innerNodeIndex = structNodeIndex[i][j + 1][k];
				bound.SetInnerIndex(innerNodeIndex);
				ghostNodeIndex = structNodeIndex[i][j - 1][k];
				bound.SetGhostIndex(ghostNodeIndex);
				boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
				bound.SetNorm(boundNorm);
				boundMap->AddBoundary("slipWall", bound);
				j = yNodeNum;
				nodeIndex = structNodeIndex[i][j][k];
				bound.SetNodeIndex(nodeIndex);
				innerNodeIndex = structNodeIndex[i][j - 1][k];
				bound.SetInnerIndex(innerNodeIndex);
				ghostNodeIndex = structNodeIndex[i][j + 1][k];
				bound.SetGhostIndex(ghostNodeIndex);
				boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
				bound.SetNorm(boundNorm);
				boundMap->AddBoundary("slipWall", bound);
			}

			//j方向
			for (j = 1; j < yNodeNum + 1; ++j)
			{
				k = 1;
				nodeIndex = structNodeIndex[i][j][k];
				bound.SetNodeIndex(nodeIndex);
				innerNodeIndex = structNodeIndex[i][j][k + 1];
				bound.SetInnerIndex(innerNodeIndex);
				ghostNodeIndex = structNodeIndex[i][j][k - 1];
				bound.SetGhostIndex(ghostNodeIndex);
				boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
				bound.SetNorm(boundNorm);
				boundMap->AddBoundary("slipWall", bound);
				k = zNodeNum;
				nodeIndex = structNodeIndex[i][j][k];
				bound.SetNodeIndex(nodeIndex);
				innerNodeIndex = structNodeIndex[i][j][k - 1];
				bound.SetInnerIndex(innerNodeIndex);
				ghostNodeIndex = structNodeIndex[i][j][k + 1];
				bound.SetGhostIndex(ghostNodeIndex);
				boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
				bound.SetNorm(boundNorm);
				boundMap->AddBoundary("slipWall", bound);
			}
		}
	}
	void GridListFactory::CreateStructGrid3D(Grid*& grid)
	{
		if (grid != nullptr)
		{
			delete grid;
		}
		Grid_Struct_3D* struct_grid = new Grid_Struct_3D();
		struct_grid->SetDimension(Dimension::three);
		struct_grid->SetIndex(0);
		struct_grid->SetLevel(0);
		struct_grid->SetName("noname");
		struct_grid->SetType(GridType::Structured);
		NodeTopo* nodeTopo = struct_grid->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		int xNodeNum = 31;
		int yNodeNum = 31;
		int zNodeNum = 31;
		int nodeNum = 0;
		double xMin = 0.0;
		double yMin = 0.0;
		double zMin = 0.0;
		double xMax = 1.0;
		double yMax = 1.0;
		double zMax = 1.0;
		double dx = (xMax - xMin) / (xNodeNum - 1);
		double dy = (yMax - yMin) / (yNodeNum - 1);
		double dz = (zMax - zMin) / (zNodeNum - 1);
		int i, j, k;
		double x, y, z;
		struct_grid->SetNodeNum((xNodeNum + 2), (yNodeNum + 2), (zNodeNum + 2));
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
					iNode = struct_grid->GetNodeIndex(i, j, k);
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
					tempI[iNode] = { struct_grid->GetNodeIndex(i - 1, j, k), struct_grid->GetNodeIndex(i, j, k), struct_grid->GetNodeIndex(i + 1, j, k) };
					tempJ[iNode] = { struct_grid->GetNodeIndex(i, j - 1, k), struct_grid->GetNodeIndex(i, j, k), struct_grid->GetNodeIndex(i, j + 1, k) };
					tempK[iNode] = { struct_grid->GetNodeIndex(i, j, k - 1), struct_grid->GetNodeIndex(i, j, k), struct_grid->GetNodeIndex(i, j, k + 1) };
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
					iNode = struct_grid->GetNodeIndex(i, j, k);
					if (i == 0)
					{
						if (j == 0)
						{
							if (k == 0)
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i + 1, j, k),  struct_grid->GetNodeIndex(i, j + 1, k),  struct_grid->GetNodeIndex(i, j, k + 1) };
							}
							else if (k == zNodeNum + 1)
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i + 1, j, k),  struct_grid->GetNodeIndex(i, j + 1, k), struct_grid->GetNodeIndex(i, j, k - 1) };
							}
							else
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i + 1, j, k),  struct_grid->GetNodeIndex(i, j + 1, k), struct_grid->GetNodeIndex(i, j, k - 1), struct_grid->GetNodeIndex(i, j, k + 1) };
							}
						}
						else if (j == yNodeNum + 1)
						{
							if (k == 0)
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i + 1, j, k), struct_grid->GetNodeIndex(i, j - 1, k), struct_grid->GetNodeIndex(i, j, k + 1) };

							}
							else if (k == zNodeNum + 1)
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i + 1, j, k), struct_grid->GetNodeIndex(i, j - 1, k),  struct_grid->GetNodeIndex(i, j, k - 1) };


							}
							else
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i + 1, j, k), struct_grid->GetNodeIndex(i, j - 1, k),  struct_grid->GetNodeIndex(i, j, k - 1), struct_grid->GetNodeIndex(i, j, k + 1) };

							}
						}
						else
						{
							if (k == 0)
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i + 1, j, k), struct_grid->GetNodeIndex(i, j - 1, k), struct_grid->GetNodeIndex(i, j + 1, k),  struct_grid->GetNodeIndex(i, j, k + 1) };

							}
							else if (k == zNodeNum + 1)
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i + 1, j, k), struct_grid->GetNodeIndex(i, j - 1, k), struct_grid->GetNodeIndex(i, j + 1, k), struct_grid->GetNodeIndex(i, j, k - 1) };

							}
							else
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i + 1, j, k), struct_grid->GetNodeIndex(i, j - 1, k), struct_grid->GetNodeIndex(i, j + 1, k), struct_grid->GetNodeIndex(i, j, k - 1), struct_grid->GetNodeIndex(i, j, k + 1) };
							}
						}
					}
					else if (i == xNodeNum + 1)
					{
						if (j == 0)
						{
							if (k == 0)
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i - 1, j, k),  struct_grid->GetNodeIndex(i, j + 1, k),  struct_grid->GetNodeIndex(i, j, k + 1) };

							}
							else if (k == zNodeNum + 1)
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i - 1, j, k),  struct_grid->GetNodeIndex(i, j + 1, k), struct_grid->GetNodeIndex(i, j, k - 1) };
							}
							else
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i - 1, j, k),  struct_grid->GetNodeIndex(i, j + 1, k), struct_grid->GetNodeIndex(i, j, k - 1), struct_grid->GetNodeIndex(i, j, k + 1) };
							}
						}
						else if (j == yNodeNum + 1)
						{
							if (k == 0)
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i - 1, j, k), struct_grid->GetNodeIndex(i, j - 1, k), struct_grid->GetNodeIndex(i, j, k + 1) };


							}
							else if (k == zNodeNum + 1)
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i - 1, j, k), struct_grid->GetNodeIndex(i, j - 1, k),  struct_grid->GetNodeIndex(i, j, k - 1) };

							}
							else
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i - 1, j, k), struct_grid->GetNodeIndex(i, j - 1, k),  struct_grid->GetNodeIndex(i, j, k - 1), struct_grid->GetNodeIndex(i, j, k + 1) };

							}
						}
						else
						{
							if (k == 0)
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i - 1, j, k), struct_grid->GetNodeIndex(i, j - 1, k), struct_grid->GetNodeIndex(i, j + 1, k),  struct_grid->GetNodeIndex(i, j, k + 1) };

							}
							else if (k == zNodeNum + 1)
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i - 1, j, k), struct_grid->GetNodeIndex(i, j - 1, k), struct_grid->GetNodeIndex(i, j + 1, k), struct_grid->GetNodeIndex(i, j, k - 1) };

							}
							else
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i - 1, j, k), struct_grid->GetNodeIndex(i, j - 1, k), struct_grid->GetNodeIndex(i, j + 1, k), struct_grid->GetNodeIndex(i, j, k - 1), struct_grid->GetNodeIndex(i, j, k + 1) };

							}
						}
					}
					else
					{
						if (j == 0)
						{
							if (k == 0)
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i - 1, j, k), struct_grid->GetNodeIndex(i + 1, j, k), struct_grid->GetNodeIndex(i, j + 1, k),  struct_grid->GetNodeIndex(i, j, k + 1) };

							}
							else if (k == zNodeNum + 1)
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i - 1, j, k), struct_grid->GetNodeIndex(i + 1, j, k), struct_grid->GetNodeIndex(i, j + 1, k), struct_grid->GetNodeIndex(i, j, k - 1) };

							}
							else
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i - 1, j, k), struct_grid->GetNodeIndex(i + 1, j, k),  struct_grid->GetNodeIndex(i, j + 1, k), struct_grid->GetNodeIndex(i, j, k - 1), struct_grid->GetNodeIndex(i, j, k + 1) };

							}
						}
						else if (j == yNodeNum + 1)
						{
							if (k == 0)
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i - 1, j, k), struct_grid->GetNodeIndex(i + 1, j, k), struct_grid->GetNodeIndex(i, j - 1, k),   struct_grid->GetNodeIndex(i, j, k + 1) };

							}
							else if (k == zNodeNum + 1)
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i - 1, j, k), struct_grid->GetNodeIndex(i + 1, j, k), struct_grid->GetNodeIndex(i, j - 1, k),  struct_grid->GetNodeIndex(i, j, k - 1) };

							}
							else
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i - 1, j, k), struct_grid->GetNodeIndex(i + 1, j, k), struct_grid->GetNodeIndex(i, j - 1, k),  struct_grid->GetNodeIndex(i, j, k - 1), struct_grid->GetNodeIndex(i, j, k + 1) };

							}
						}
						else
						{
							if (k == 0)
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i - 1, j, k), struct_grid->GetNodeIndex(i + 1, j, k), struct_grid->GetNodeIndex(i, j - 1, k), struct_grid->GetNodeIndex(i, j + 1, k),  struct_grid->GetNodeIndex(i, j, k + 1) };
							}
							else if (k == zNodeNum + 1)
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i - 1, j, k), struct_grid->GetNodeIndex(i + 1, j, k), struct_grid->GetNodeIndex(i, j - 1, k), struct_grid->GetNodeIndex(i, j + 1, k), struct_grid->GetNodeIndex(i, j, k - 1) };

							}
							else
							{
								nodeNeighbor[iNode] = { struct_grid->GetNodeIndex(i - 1, j, k), struct_grid->GetNodeIndex(i + 1, j, k), struct_grid->GetNodeIndex(i, j - 1, k), struct_grid->GetNodeIndex(i, j + 1, k), struct_grid->GetNodeIndex(i, j, k - 1), struct_grid->GetNodeIndex(i, j, k + 1) };
							}
						}
					}
				}
			}
		}
		struct_grid->SetInnerNodeNum((xNodeNum - 2) * (yNodeNum - 2) * (zNodeNum - 2));
		struct_grid->SetTotalNodeNum((xNodeNum + 2) * (yNodeNum + 2) * (zNodeNum + 2));
		struct_grid->SetBoundNodeNum(xNodeNum * yNodeNum * zNodeNum - struct_grid->GetInnerNodeNum());
		CellTopo* cellTopo = struct_grid->GetCellTopo();
		auto& cell2node = cellTopo->GetNodeIndex();
		cell2node.resize((xNodeNum + 1) * (yNodeNum + 1) * (zNodeNum + 1));
		int iCell;
		for (k = 0; k < zNodeNum + 1; ++k)
		{
			for (j = 0; j < yNodeNum + 1; ++j)
			{

				for (i = 0; i < xNodeNum + 1; ++i)
				{
					iCell = struct_grid->GetCellIndex(i, j, k);
					cell2node[iCell] = (IArray{ struct_grid->GetNodeIndex(i, j, k), struct_grid->GetNodeIndex(i + 1, j, k), struct_grid->GetNodeIndex(i + 1, j + 1, k), struct_grid->GetNodeIndex(i, j + 1, k), struct_grid->GetNodeIndex(i, j, k + 1), struct_grid->GetNodeIndex(i + 1, j, k + 1), struct_grid->GetNodeIndex(i + 1, j + 1, k + 1), struct_grid->GetNodeIndex(i, j + 1, k + 1) });
				}
			}
		}

		auto boundMap = struct_grid->GetBoundaryMap();
		int nodeIndex, innerNodeIndex, ghostNodeIndex;
		DVector3D boundNorm;
		Boundary bound;
		//i方向两个面分别为入口和出口
		for (k = 1; k < zNodeNum + 1; ++k)
		{

			for (j = 1; j < yNodeNum + 1; ++j)
			{
				i = 1;
				nodeIndex = struct_grid->GetNodeIndex(i, j, k);
				bound.SetNodeIndex(nodeIndex);
				innerNodeIndex = struct_grid->GetNodeIndex(i + 1, j, k);
				bound.SetInnerIndex(innerNodeIndex);
				ghostNodeIndex = struct_grid->GetNodeIndex(i - 1, j, k);
				bound.SetGhostIndex(ghostNodeIndex);
				boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
				bound.SetNorm(boundNorm);
				boundMap->AddBoundary("inlet", bound);

				i = xNodeNum;
				nodeIndex = struct_grid->GetNodeIndex(i, j, k);
				bound.SetNodeIndex(nodeIndex);
				innerNodeIndex = struct_grid->GetNodeIndex(i - 1, j, k);
				bound.SetInnerIndex(innerNodeIndex);
				ghostNodeIndex = struct_grid->GetNodeIndex(i + 1, j, k);
				bound.SetGhostIndex(ghostNodeIndex);
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
				nodeIndex = struct_grid->GetNodeIndex(i, j, k);
				bound.SetNodeIndex(nodeIndex);
				innerNodeIndex = struct_grid->GetNodeIndex(i, j + 1, k);
				bound.SetInnerIndex(innerNodeIndex);
				ghostNodeIndex = struct_grid->GetNodeIndex(i, j - 1, k);
				bound.SetGhostIndex(ghostNodeIndex);
				boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
				bound.SetNorm(boundNorm);
				boundMap->AddBoundary("wall", bound);
				j = yNodeNum;
				nodeIndex = struct_grid->GetNodeIndex(i, j, k);
				bound.SetNodeIndex(nodeIndex);
				innerNodeIndex = struct_grid->GetNodeIndex(i, j - 1, k);
				bound.SetInnerIndex(innerNodeIndex);
				ghostNodeIndex = struct_grid->GetNodeIndex(i, j + 1, k);
				bound.SetGhostIndex(ghostNodeIndex);
				boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
				bound.SetNorm(boundNorm);
				boundMap->AddBoundary("wall", bound);
			}

			//j方向
			for (j = 1; j < yNodeNum + 1; ++j)
			{
				k = 1;
				nodeIndex = struct_grid->GetNodeIndex(i, j, k);
				bound.SetNodeIndex(nodeIndex);
				innerNodeIndex = struct_grid->GetNodeIndex(i, j, k + 1);
				bound.SetInnerIndex(innerNodeIndex);
				ghostNodeIndex = struct_grid->GetNodeIndex(i, j, k - 1);
				bound.SetGhostIndex(ghostNodeIndex);
				boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
				bound.SetNorm(boundNorm);
				boundMap->AddBoundary("wall", bound);
				k = zNodeNum;
				nodeIndex = struct_grid->GetNodeIndex(i, j, k);
				bound.SetNodeIndex(nodeIndex);
				innerNodeIndex = struct_grid->GetNodeIndex(i, j, k - 1);
				bound.SetInnerIndex(innerNodeIndex);
				ghostNodeIndex = struct_grid->GetNodeIndex(i, j, k + 1);
				bound.SetGhostIndex(ghostNodeIndex);
				boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
				bound.SetNorm(boundNorm);
				boundMap->AddBoundary("wall", bound);
			}
		}
		grid = dynamic_cast<Grid*>(struct_grid);
	}
	void GridListFactory::CreateStructGrid2D(Grid*& grid)
	{
		if (grid != nullptr)
		{
			delete[] grid;
		}
		Grid_Struct_2D* struct_grid = new Grid_Struct_2D();
		struct_grid->SetDimension(Dimension::two);
		struct_grid->SetIndex(0);
		struct_grid->SetLevel(0);
		struct_grid->SetName("fnfdm-structred-grid-2d");
		struct_grid->SetType(GridType::Unkown);
		NodeTopo* nodeTopo = struct_grid->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		int xNodeNum = 101;
		int yNodeNum = 101;
		int nodeNum = 0;
		double xMin = 0.0;
		double yMin = 0.0;
		double xMax = 1.0;
		double yMax = 1.0;
		double dx = (xMax - xMin) / (xNodeNum - 1);
		double dy = (yMax - yMin) / (yNodeNum - 1);
		int i, j;
		double x, y;
		struct_grid->SetNodeNum((xNodeNum + 2), (yNodeNum + 2));
		nodeCoord.resize((xNodeNum + 2) * (yNodeNum + 2));
		auto& nodeType = nodeTopo->GetType();
		nodeType.resize((xNodeNum + 2) * (yNodeNum + 2));
		auto& tempI = nodeTopo->GetTemplateI();
		auto& tempJ = nodeTopo->GetTemplateJ();
		tempI.resize((xNodeNum + 2) * (yNodeNum + 2));
		tempJ.resize((xNodeNum + 2) * (yNodeNum + 2));
		//求节点下标的lamda表达式
		auto GetNodeIndex = [&](int i, int j)
			{
				return j * (xNodeNum + 2) + i;
			};
		int iNode = 0;
		for (j = 0; j < yNodeNum + 2; ++j)
		{
			y = j * dy + yMin - dy;
			for (i = 0; i < xNodeNum + 2; ++i)
			{
				x = i * dx + xMin - dx;
				iNode = struct_grid->GetNodeIndex(i, j);
				nodeCoord[iNode] = { x,y,0.0 };
				if (i == 0 || j == 0 || i == xNodeNum + 1 || j == yNodeNum + 1)
				{
					nodeType[iNode] = NodeType::ghost;
					continue;
				}
				if (i == 1 || j == 1 || i == xNodeNum || j == yNodeNum)
				{
					nodeType[iNode] = NodeType::inlet;
				}
				tempI[iNode] = { GetNodeIndex(i - 1, j), GetNodeIndex(i, j), GetNodeIndex(i + 1, j) };
				tempJ[iNode] = { GetNodeIndex(i, j - 1), GetNodeIndex(i, j), GetNodeIndex(i, j + 1) };
			}
		}
		auto& nodeNeighbor = nodeTopo->GetNeighborCloud();
		nodeNeighbor.resize((xNodeNum + 2) * (yNodeNum + 2));
		for (j = 0; j < yNodeNum + 2; ++j)
		{
			for (i = 0; i < xNodeNum + 2; ++i)
			{
				iNode = GetNodeIndex(i, j);
				if (i == 0)
				{
					if (j == 0)
					{
						nodeNeighbor[iNode] = { GetNodeIndex(i + 1, j),  GetNodeIndex(i, j + 1) };
					}
					else if (j == yNodeNum + 1)
					{
						nodeNeighbor[iNode] = { GetNodeIndex(i + 1, j),  GetNodeIndex(i, j - 1) };
					}
					else
					{
						nodeNeighbor[iNode] = { GetNodeIndex(i + 1, j),  GetNodeIndex(i, j - 1), GetNodeIndex(i, j + 1) };
					}
				}
				else if (i == xNodeNum + 1)
				{
					if (j == 0)
					{
						nodeNeighbor[iNode] = { GetNodeIndex(i - 1, j),  GetNodeIndex(i, j + 1) };

					}
					else if (j == yNodeNum + 1)
					{
						nodeNeighbor[iNode] = { GetNodeIndex(i - 1, j),  GetNodeIndex(i, j - 1) };
					}
				}
				else
				{
					if (j == 0)
					{
						nodeNeighbor[iNode] = { GetNodeIndex(i - 1, j), GetNodeIndex(i + 1, j), GetNodeIndex(i, j + 1) };

					}
					else if (j == yNodeNum + 1)
					{
						nodeNeighbor[iNode] = { GetNodeIndex(i - 1, j), GetNodeIndex(i + 1, j), GetNodeIndex(i, j - 1) };

					}
					else
					{
						nodeNeighbor[iNode] = { GetNodeIndex(i - 1, j), GetNodeIndex(i + 1, j), GetNodeIndex(i, j - 1), GetNodeIndex(i, j + 1) };
					}
				}
			}
		}
		struct_grid->SetInnerNodeNum((xNodeNum - 2) * (yNodeNum - 2));
		struct_grid->SetBoundNodeNum(xNodeNum * yNodeNum - struct_grid->GetInnerNodeNum());
		CellTopo* cellTopo = struct_grid->GetCellTopo();
		auto& cell2node = cellTopo->GetNodeIndex();
		cell2node.resize((xNodeNum + 1) * (yNodeNum + 1));
		int iCell;
		for (j = 0; j < yNodeNum + 1; ++j)
		{
			for (i = 0; i < xNodeNum + 1; ++i)
			{
				iCell = j * (xNodeNum + 1) + i;
				cell2node[iCell] = (IArray{ GetNodeIndex(i, j), GetNodeIndex(i + 1, j), GetNodeIndex(i + 1, j + 1), GetNodeIndex(i, j + 1) });
			}
		}

		BoundaryMap* boundMap = struct_grid->GetBoundaryMap();
		int nodeIndex, innerNodeIndex, ghostNodeIndex;
		DVector3D boundNorm;
		Boundary bound;
		//i方向两个面分别为入口和出口
		for (j = 1; j < yNodeNum + 1; ++j)
		{
			i = 1;
			nodeIndex = GetNodeIndex(i, j);
			bound.SetNodeIndex(nodeIndex);
			innerNodeIndex = GetNodeIndex(i + 1, j);
			bound.SetInnerIndex(innerNodeIndex);
			ghostNodeIndex = GetNodeIndex(i - 1, j);
			bound.SetGhostIndex(ghostNodeIndex);
			boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
			bound.SetNorm(boundNorm);
			boundMap->AddBoundary("inlet", bound);

			i = xNodeNum;
			nodeIndex = GetNodeIndex(i, j);
			bound.SetNodeIndex(nodeIndex);
			innerNodeIndex = GetNodeIndex(i - 1, j);
			bound.SetInnerIndex(innerNodeIndex);
			ghostNodeIndex = GetNodeIndex(i + 1, j);
			bound.SetGhostIndex(ghostNodeIndex);
			boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
			bound.SetNorm(boundNorm);
			boundMap->AddBoundary("outlet", bound);
		}
		//j方向两个面均为壁面
		for (i = 1; i < xNodeNum + 1; ++i)
		{
			j = 1;
			nodeIndex = GetNodeIndex(i, j);
			bound.SetNodeIndex(nodeIndex);
			innerNodeIndex = GetNodeIndex(i, j + 1);
			bound.SetInnerIndex(innerNodeIndex);
			ghostNodeIndex = GetNodeIndex(i, j - 1);
			bound.SetGhostIndex(ghostNodeIndex);
			boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
			bound.SetNorm(boundNorm);
			boundMap->AddBoundary("outlet", bound);
			j = yNodeNum;
			nodeIndex = GetNodeIndex(i, j);
			bound.SetNodeIndex(nodeIndex);
			innerNodeIndex = GetNodeIndex(i, j - 1);
			bound.SetInnerIndex(innerNodeIndex);
			ghostNodeIndex = GetNodeIndex(i, j + 1);
			bound.SetGhostIndex(ghostNodeIndex);
			boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
			bound.SetNorm(boundNorm);
			boundMap->AddBoundary("inlet", bound);
		}
		grid = dynamic_cast<Grid*>(struct_grid);
	}
	void GridListFactory::CreateGridFNFDM2D(Grid*& grid)
	{
		if (grid != nullptr)
		{
			delete[] grid;
		}
		grid = new Grid();
		grid->SetDimension(Dimension::two);
		grid->SetIndex(0);
		grid->SetLevel(0);
		grid->SetName("fnfdm-structred-grid");
		grid->SetType(GridType::Flexible);
		NodeTopo* nodeTopo = grid->GetNodeTopo();

		auto& nodeCoord = nodeTopo->GetCoordinate();
		int xNodeNum = 31;
		int yNodeNum = 31;
		int nodeNum = 0;
		double xMin = 0.0;
		double yMin = 0.0;
		double xMax = 1.0;
		double yMax = 1.0;
		double dx = (xMax - xMin) / (xNodeNum - 1);
		double dy = (yMax - yMin) / (yNodeNum - 1);
		int i, j;
		double x, y;
		grid->SetTotalNodeNum((xNodeNum + 2) * (yNodeNum + 2));
		nodeCoord.resize((xNodeNum + 2) * (yNodeNum + 2));
		auto& nodeType = nodeTopo->GetType();
		nodeType.resize((xNodeNum + 2) * (yNodeNum + 2));
		auto& tempI = nodeTopo->GetTemplateI();
		auto& tempJ = nodeTopo->GetTemplateJ();
		tempI.resize((xNodeNum + 2) * (yNodeNum + 2));
		tempJ.resize((xNodeNum + 2) * (yNodeNum + 2));
		//求节点下标的lamda表达式
		auto GetNodeIndex = [&](int i, int j)
			{
				return j * (xNodeNum + 2) + i;
			};
		int iNode = 0;
		for (j = 0; j < yNodeNum + 2; ++j)
		{
			y = j * dy + yMin - dy;
			for (i = 0; i < xNodeNum + 2; ++i)
			{
				x = i * dx + xMin - dx;
				iNode = GetNodeIndex(i, j);
				nodeCoord[iNode] = { x,y,0.0 };
				if (i == 0 || j == 0 || i == xNodeNum + 1 || j == yNodeNum + 1)
				{
					nodeType[iNode] = NodeType::ghost;
					continue;
				}
				if (i == 1 || j == 1 || i == xNodeNum || j == yNodeNum)
				{
					nodeType[iNode] = NodeType::inlet;
				}
				tempI[iNode] = { GetNodeIndex(i - 1, j), GetNodeIndex(i, j), GetNodeIndex(i + 1, j) };
				tempJ[iNode] = { GetNodeIndex(i, j - 1), GetNodeIndex(i, j), GetNodeIndex(i, j + 1) };
			}
		}
		auto& nodeNeighbor = nodeTopo->GetNeighborCloud();
		nodeNeighbor.resize((xNodeNum + 2) * (yNodeNum + 2));
		for (j = 0; j < yNodeNum + 2; ++j)
		{
			for (i = 0; i < xNodeNum + 2; ++i)
			{
				iNode = GetNodeIndex(i, j);
				if (i == 0)
				{
					if (j == 0)
					{
						nodeNeighbor[iNode] = { GetNodeIndex(i + 1, j),  GetNodeIndex(i, j + 1) };
					}
					else if (j == yNodeNum + 1)
					{
						nodeNeighbor[iNode] = { GetNodeIndex(i + 1, j),  GetNodeIndex(i, j - 1) };
					}
					else
					{
						nodeNeighbor[iNode] = { GetNodeIndex(i + 1, j),  GetNodeIndex(i, j - 1), GetNodeIndex(i, j + 1) };
					}
				}
				else if (i == xNodeNum + 1)
				{
					if (j == 0)
					{
						nodeNeighbor[iNode] = { GetNodeIndex(i - 1, j),  GetNodeIndex(i, j + 1) };

					}
					else if (j == yNodeNum + 1)
					{
						nodeNeighbor[iNode] = { GetNodeIndex(i - 1, j),  GetNodeIndex(i, j - 1) };
					}
				}
				else
				{
					if (j == 0)
					{
						nodeNeighbor[iNode] = { GetNodeIndex(i - 1, j), GetNodeIndex(i + 1, j), GetNodeIndex(i, j + 1) };

					}
					else if (j == yNodeNum + 1)
					{
						nodeNeighbor[iNode] = { GetNodeIndex(i - 1, j), GetNodeIndex(i + 1, j), GetNodeIndex(i, j - 1) };

					}
					else
					{
						nodeNeighbor[iNode] = { GetNodeIndex(i - 1, j), GetNodeIndex(i + 1, j), GetNodeIndex(i, j - 1), GetNodeIndex(i, j + 1) };
					}
				}
			}
		}
		grid->SetInnerNodeNum((xNodeNum - 2) * (yNodeNum - 2));
		grid->SetBoundNodeNum(xNodeNum * yNodeNum - grid->GetInnerNodeNum());
		CellTopo* cellTopo = grid->GetCellTopo();
		auto& cell2node = cellTopo->GetNodeIndex();
		cell2node.resize((xNodeNum + 1) * (yNodeNum + 1));
		int iCell;
		for (j = 0; j < yNodeNum + 1; ++j)
		{
			for (i = 0; i < xNodeNum + 1; ++i)
			{
				iCell = j * (xNodeNum + 1) + i;
				cell2node[iCell] = (IArray{ GetNodeIndex(i, j), GetNodeIndex(i + 1, j), GetNodeIndex(i + 1, j + 1), GetNodeIndex(i, j + 1) });
			}
		}

		BoundaryMap* boundMap = grid->GetBoundaryMap();
		int nodeIndex, innerNodeIndex, ghostNodeIndex;
		DVector3D boundNorm;
		Boundary bound;
		//i方向两个面分别为入口和出口
		for (j = 1; j < yNodeNum + 1; ++j)
		{
			i = 1;
			nodeIndex = GetNodeIndex(i, j);
			bound.SetNodeIndex(nodeIndex);
			innerNodeIndex = GetNodeIndex(i + 1, j);
			bound.SetInnerIndex(innerNodeIndex);
			ghostNodeIndex = GetNodeIndex(i - 1, j);
			bound.SetGhostIndex(ghostNodeIndex);
			boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
			bound.SetNorm(boundNorm);
			boundMap->AddBoundary("inlet", bound);

			i = xNodeNum;
			nodeIndex = GetNodeIndex(i, j);
			bound.SetNodeIndex(nodeIndex);
			innerNodeIndex = GetNodeIndex(i - 1, j);
			bound.SetInnerIndex(innerNodeIndex);
			ghostNodeIndex = GetNodeIndex(i + 1, j);
			bound.SetGhostIndex(ghostNodeIndex);
			boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
			bound.SetNorm(boundNorm);
			boundMap->AddBoundary("outlet", bound);
		}
		//j方向两个面均为壁面
		for (i = 1; i < xNodeNum + 1; ++i)
		{
			j = 1;
			nodeIndex = GetNodeIndex(i, j);
			bound.SetNodeIndex(nodeIndex);
			innerNodeIndex = GetNodeIndex(i, j + 1);
			bound.SetInnerIndex(innerNodeIndex);
			ghostNodeIndex = GetNodeIndex(i, j - 1);
			bound.SetGhostIndex(ghostNodeIndex);
			boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
			bound.SetNorm(boundNorm);
			boundMap->AddBoundary("slipWall", bound);
			j = yNodeNum;
			nodeIndex = GetNodeIndex(i, j);
			bound.SetNodeIndex(nodeIndex);
			innerNodeIndex = GetNodeIndex(i, j - 1);
			bound.SetInnerIndex(innerNodeIndex);
			ghostNodeIndex = GetNodeIndex(i, j + 1);
			bound.SetGhostIndex(ghostNodeIndex);
			boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
			bound.SetNorm(boundNorm);
			boundMap->AddBoundary("slipWall", bound);
		}
	}
}