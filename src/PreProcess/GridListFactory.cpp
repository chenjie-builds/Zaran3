#include "GridListFactory.h"
#include"StructGrid.h"
#include "grid.h"
#include "log.h"
namespace zaran
{

	void GridListFactory::Create(Ptr<GridList>& gridList)
	{
		std::string createMethod = GlobalData::GetString("createGridMethod");
		if (createMethod == "Plot3D")
		{
			//CreateByTest(gridList);
			CreateStructGrid3D(gridList);
		}
		else
		{
			ZaranLog::warn("Unsupported Grid Generate Method:{}, Please Check!", createMethod);
			system("pause");
		}
	}

	void GridListFactory::ReadPlot3D(Ptr<GridList>& gridList)
	{
		if (gridList.get() == nullptr)
		{
			gridList = std::make_shared<GridList>();
		}
		Ptr<Grid>testGrid = std::make_shared<Grid>();
		testGrid->SetDimension(Dimension::three);
		testGrid->SetIndex(0);
		testGrid->SetLevel(0);
		testGrid->SetName("noname");
		testGrid->SetType(GridType::NOTSET);
		gridList->AddGrid(testGrid);
	}

	void GridListFactory::CreateByTest(Ptr<GridList>& gridList)
	{
		if (gridList.get() == nullptr)
		{
			gridList = std::make_shared<GridList>();
		}
		Ptr<Grid>testGrid = std::make_shared<Grid>();
		testGrid->SetDimension(Dimension::three);
		testGrid->SetIndex(0);
		testGrid->SetLevel(0);
		testGrid->SetName("noname");
		testGrid->SetType(GridType::NOTSET);
		auto& nodeTopo = testGrid->GetNodeTopo();
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
		testGrid->SetInnerNodeNum((xNodeNum - 2) * (yNodeNum - 2) * (zNodeNum - 2));
		testGrid->SetTotalNodeNum((xNodeNum + 2) * (yNodeNum + 2) * (zNodeNum + 2));
		testGrid->SetBoundNodeNum(xNodeNum * yNodeNum * zNodeNum - testGrid->GetInnerNodeNum());
		testGrid->SetInterNodeInfo(std::make_shared<InterNodeInfo>());
		testGrid->SetBoundaryMap(std::make_shared<BoundaryMap>());

		auto& cellTopo = testGrid->GetCellTopo();
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

		auto boundMap = testGrid->GetBoundaryMap();
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
				bound.SetInnerNodeIndex(innerNodeIndex);
				ghostNodeIndex = structNodeIndex[i - 1][j][k];
				bound.SetGhostNodeIndex(ghostNodeIndex);
				boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
				bound.SetNorm(boundNorm);
				boundMap->AddBoundary("inlet", bound);

				i = xNodeNum;
				nodeIndex = structNodeIndex[i][j][k];
				bound.SetNodeIndex(nodeIndex);
				innerNodeIndex = structNodeIndex[i - 1][j][k];
				bound.SetInnerNodeIndex(innerNodeIndex);
				ghostNodeIndex = structNodeIndex[i + 1][j][k];
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
				nodeIndex = structNodeIndex[i][j][k];
				bound.SetNodeIndex(nodeIndex);
				innerNodeIndex = structNodeIndex[i][j + 1][k];
				bound.SetInnerNodeIndex(innerNodeIndex);
				ghostNodeIndex = structNodeIndex[i][j - 1][k];
				bound.SetGhostNodeIndex(ghostNodeIndex);
				boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
				bound.SetNorm(boundNorm);
				boundMap->AddBoundary("wall", bound);
				j = yNodeNum;
				nodeIndex = structNodeIndex[i][j][k];
				bound.SetNodeIndex(nodeIndex);
				innerNodeIndex = structNodeIndex[i][j - 1][k];
				bound.SetInnerNodeIndex(innerNodeIndex);
				ghostNodeIndex = structNodeIndex[i][j + 1][k];
				bound.SetGhostNodeIndex(ghostNodeIndex);
				boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
				bound.SetNorm(boundNorm);
				boundMap->AddBoundary("wall", bound);
			}

			//j方向
			for (j = 1; j < yNodeNum + 1; ++j)
			{
				k = 1;
				nodeIndex = structNodeIndex[i][j][k];
				bound.SetNodeIndex(nodeIndex);
				innerNodeIndex = structNodeIndex[i][j][k + 1];
				bound.SetInnerNodeIndex(innerNodeIndex);
				ghostNodeIndex = structNodeIndex[i][j][k - 1];
				bound.SetGhostNodeIndex(ghostNodeIndex);
				boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
				bound.SetNorm(boundNorm);
				boundMap->AddBoundary("wall", bound);
				k = zNodeNum;
				nodeIndex = structNodeIndex[i][j][k];
				bound.SetNodeIndex(nodeIndex);
				innerNodeIndex = structNodeIndex[i][j][k - 1];
				bound.SetInnerNodeIndex(innerNodeIndex);
				ghostNodeIndex = structNodeIndex[i][j][k + 1];
				bound.SetGhostNodeIndex(ghostNodeIndex);
				boundNorm = nodeCoord[ghostNodeIndex] - nodeCoord[nodeIndex];
				bound.SetNorm(boundNorm);
				boundMap->AddBoundary("wall", bound);
			}
		}
		gridList->AddGrid(testGrid);

	}
	void GridListFactory::CreateStructGrid3D(Ptr<GridList>& gridList)
	{
		if (gridList.get() == nullptr)
		{
			gridList = std::make_shared<GridList>();
		}

		Ptr<StructGrid>grid = std::make_shared<StructGrid>();
		grid->SetDimension(Dimension::three);
		grid->SetIndex(0);
		grid->SetLevel(0);
		grid->SetName("noname");
		grid->SetType(GridType::NOTSET);
		auto& nodeTopo = grid->GetNodeTopo();
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
}