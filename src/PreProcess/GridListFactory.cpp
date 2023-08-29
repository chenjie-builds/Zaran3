#include "GridListFactory.h"
#include "grid.h"
#include "log.h"
using namespace zaran;
void GridListFactory::Create(Ptr<GridList>& gridList)
{
	std::string createMethod = GlobalData::GetString("createGridMethod");
	if (createMethod == "Plot3D")
	{
		CreateByTest(gridList);
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
	auto& node = testGrid->GetNodeTopoInfo();
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
	node.resize((xNodeNum + 2) * (yNodeNum + 2) * (zNodeNum + 2));
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
	IArray nodeTempI(3), nodeTempJ(3), nodeTempK(3);
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
				auto& currentNode = node[iNode];
				currentNode.SetCoordinate(x, y, z);
				currentNode.SetTag(1);
				if (i == 0 || j == 0 || k == 0 || i == xNodeNum + 1 || j == yNodeNum + 1 || k == zNodeNum + 1)
				{
					currentNode.SetType(NodeType::ghost);
					continue;
				}
				if (i == 1 || j == 1 || k == 1 || i == xNodeNum || j == yNodeNum || k == zNodeNum)
				{
					currentNode.SetType(NodeType::inlet);
				}
				nodeTempI = { structNodeIndex[i - 1][j][k], structNodeIndex[i][j][k], structNodeIndex[i + 1][j][k] };
				nodeTempJ = { structNodeIndex[i][j - 1][k], structNodeIndex[i][j][k], structNodeIndex[i][j + 1][k] };
				nodeTempK = { structNodeIndex[i][j][k - 1], structNodeIndex[i][j][k], structNodeIndex[i][j][k + 1] };
				currentNode.SetNeighborTemplateI(nodeTempI);
				currentNode.SetNeighborTemplateJ(nodeTempJ);
				currentNode.SetNeighborTemplateK(nodeTempK);
			}
		}
	}
	IArray neiborCloud(6);
	for (k = 0; k < zNodeNum + 2; ++k)
	{
		for (j = 0; j < yNodeNum + 2; ++j)
		{
			for (i = 0; i < xNodeNum + 2; ++i)
			{
				int iNode = structNodeIndex[i][j][k];
				auto& currentNode = node[iNode];
				if (i == 0)
				{
					if (j == 0)
					{
						if (k == 0)
						{
							neiborCloud = { structNodeIndex[i + 1][j][k],  structNodeIndex[i][j + 1][k],  structNodeIndex[i][j][k + 1] };
						}
						else if (k == zNodeNum + 1)
						{
							neiborCloud = { structNodeIndex[i + 1][j][k],  structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k - 1] };
						}
						else
						{
							neiborCloud = { structNodeIndex[i + 1][j][k],  structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k - 1], structNodeIndex[i][j][k + 1] };
						}
					}
					else if (j == yNodeNum + 1)
					{
						if (k == 0)
						{
							neiborCloud = { structNodeIndex[i + 1][j][k], structNodeIndex[i][j - 1][k], structNodeIndex[i][j][k + 1] };

						}
						else if (k == zNodeNum + 1)
						{
							neiborCloud = { structNodeIndex[i + 1][j][k], structNodeIndex[i][j - 1][k],  structNodeIndex[i][j][k - 1] };

						}
						else
						{
							neiborCloud = { structNodeIndex[i + 1][j][k], structNodeIndex[i][j - 1][k], structNodeIndex[i][j][k - 1], structNodeIndex[i][j][k + 1] };

						}
					}
					else
					{
						if (k == 0)
						{
							neiborCloud = { structNodeIndex[i + 1][j][k], structNodeIndex[i][j - 1][k], structNodeIndex[i][j + 1][k],  structNodeIndex[i][j][k + 1] };

						}
						else if (k == zNodeNum + 1)
						{
							neiborCloud = { structNodeIndex[i + 1][j][k], structNodeIndex[i][j - 1][k], structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k - 1] };

						}
						else
						{
							neiborCloud = { structNodeIndex[i + 1][j][k], structNodeIndex[i][j - 1][k], structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k - 1], structNodeIndex[i][j][k + 1] };

						}
					}
				}
				else if (i == xNodeNum + 1)
				{
					if (j == 0)
					{
						if (k == 0)
						{
							neiborCloud = { structNodeIndex[i - 1][j][k],  structNodeIndex[i][j + 1][k],  structNodeIndex[i][j][k + 1] };
						}
						else if (k == zNodeNum + 1)
						{
							neiborCloud = { structNodeIndex[i - 1][j][k],  structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k - 1] };
						}
						else
						{
							neiborCloud = { structNodeIndex[i - 1][j][k],  structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k - 1], structNodeIndex[i][j][k + 1] };
						}
					}
					else if (j == yNodeNum + 1)
					{
						if (k == 0)
						{
							neiborCloud = { structNodeIndex[i - 1][j][k], structNodeIndex[i][j - 1][k],  structNodeIndex[i][j][k + 1] };

						}
						else if (k == zNodeNum + 1)
						{
							neiborCloud = { structNodeIndex[i - 1][j][k], structNodeIndex[i][j - 1][k],  structNodeIndex[i][j][k - 1] };

						}
						else
						{
							neiborCloud = { structNodeIndex[i - 1][j][k], structNodeIndex[i][j - 1][k],  structNodeIndex[i][j][k - 1], structNodeIndex[i][j][k + 1] };

						}
					}
					else
					{
						if (k == 0)
						{
							neiborCloud = { structNodeIndex[i - 1][j][k], structNodeIndex[i][j - 1][k], structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k + 1] };

						}
						else if (k == zNodeNum + 1)
						{
							neiborCloud = { structNodeIndex[i - 1][j][k], structNodeIndex[i][j - 1][k], structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k - 1] };

						}
						else
						{
							neiborCloud = { structNodeIndex[i - 1][j][k], structNodeIndex[i][j - 1][k], structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k - 1], structNodeIndex[i][j][k + 1] };

						}
					}
				}
				else
				{
					if (j == 0)
					{
						if (k == 0)
						{
							neiborCloud = { structNodeIndex[i - 1][j][k], structNodeIndex[i + 1][j][k], structNodeIndex[i][j + 1][k],structNodeIndex[i][j][k + 1] };

						}
						else if (k == zNodeNum + 1)
						{
							neiborCloud = { structNodeIndex[i - 1][j][k], structNodeIndex[i + 1][j][k], structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k - 1] };

						}
						else
						{
							neiborCloud = { structNodeIndex[i - 1][j][k], structNodeIndex[i + 1][j][k],  structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k - 1], structNodeIndex[i][j][k + 1] };

						}
					}
					else if (j == yNodeNum + 1)
					{
						if (k == 0)
						{
							neiborCloud = { structNodeIndex[i - 1][j][k], structNodeIndex[i + 1][j][k], structNodeIndex[i][j - 1][k],   structNodeIndex[i][j][k + 1] };

						}
						else if (k == zNodeNum + 1)
						{
							neiborCloud = { structNodeIndex[i - 1][j][k], structNodeIndex[i + 1][j][k], structNodeIndex[i][j - 1][k],  structNodeIndex[i][j][k - 1] };

						}
						else
						{
							neiborCloud = { structNodeIndex[i - 1][j][k], structNodeIndex[i + 1][j][k], structNodeIndex[i][j - 1][k],  structNodeIndex[i][j][k - 1], structNodeIndex[i][j][k + 1] };

						}
					}
					else
					{
						if (k == 0)
						{
							neiborCloud = { structNodeIndex[i - 1][j][k], structNodeIndex[i + 1][j][k], structNodeIndex[i][j - 1][k], structNodeIndex[i][j + 1][k],  structNodeIndex[i][j][k + 1] };

						}
						else if (k == zNodeNum + 1)
						{
							neiborCloud = { structNodeIndex[i - 1][j][k], structNodeIndex[i + 1][j][k], structNodeIndex[i][j - 1][k], structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k - 1] };

						}
						else
						{
							neiborCloud = { structNodeIndex[i - 1][j][k], structNodeIndex[i + 1][j][k], structNodeIndex[i][j - 1][k], structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k - 1], structNodeIndex[i][j][k + 1] };

						}
					}
				}
				currentNode.SetNeighborCloud(neiborCloud);
			}
		}
	}
	testGrid->SetInnerNodeNumber((xNodeNum - 2) * (yNodeNum - 2) * (zNodeNum - 2));
	testGrid->SetTotalNodeNumber((xNodeNum + 2) * (yNodeNum + 2) * (zNodeNum + 2));
	testGrid->SetBoundNodeNumber(xNodeNum * yNodeNum * zNodeNum - testGrid->GetInnerNodeNum());
	testGrid->SetInterNodeInfo(std::make_shared<InterNodeInfo>());
	testGrid->SetBoundaryMap(std::make_shared<BoundaryMap>());

	auto& cell = testGrid->GetCellTopoInfo();
	cell.resize((xNodeNum - 1) * (yNodeNum - 1) * (zNodeNum - 1));
	auto& iterCell = cell.begin();
	for (k = 1; k < zNodeNum; ++k)
	{
		for (j = 1; j < yNodeNum; ++j)
		{

			for (i = 1; i < xNodeNum; ++i)
			{
				iterCell->SetNode(IArray{ structNodeIndex[i][j][k], structNodeIndex[i + 1][j][k], structNodeIndex[i + 1][j + 1][k], structNodeIndex[i][j + 1][k], structNodeIndex[i][j][k + 1], structNodeIndex[i + 1][j][k + 1], structNodeIndex[i + 1][j + 1][k + 1], structNodeIndex[i][j + 1][k + 1] });
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
			boundNorm = node[ghostNodeIndex].GetCoordinate() - node[nodeIndex].GetCoordinate();
			bound.SetNorm(boundNorm);
			boundMap->AddBoundary("inlet", bound);

			i = xNodeNum;
			nodeIndex = structNodeIndex[i][j][k];
			bound.SetNodeIndex(nodeIndex);
			innerNodeIndex = structNodeIndex[i - 1][j][k];
			bound.SetInnerNodeIndex(innerNodeIndex);
			ghostNodeIndex = structNodeIndex[i + 1][j][k];
			bound.SetGhostNodeIndex(ghostNodeIndex);
			boundNorm = node[ghostNodeIndex].GetCoordinate() - node[nodeIndex].GetCoordinate();
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
			boundNorm = node[ghostNodeIndex].GetCoordinate() - node[nodeIndex].GetCoordinate();
			bound.SetNorm(boundNorm);
			boundMap->AddBoundary("wall", bound);
			j = yNodeNum;
			nodeIndex = structNodeIndex[i][j][k];
			bound.SetNodeIndex(nodeIndex);
			innerNodeIndex = structNodeIndex[i][j - 1][k];
			bound.SetInnerNodeIndex(innerNodeIndex);
			ghostNodeIndex = structNodeIndex[i][j + 1][k];
			bound.SetGhostNodeIndex(ghostNodeIndex);
			boundNorm = node[ghostNodeIndex].GetCoordinate() - node[nodeIndex].GetCoordinate();
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
			boundNorm = node[ghostNodeIndex].GetCoordinate() - node[nodeIndex].GetCoordinate();
			bound.SetNorm(boundNorm);
			boundMap->AddBoundary("wall", bound);
			k = zNodeNum;
			nodeIndex = structNodeIndex[i][j][k];
			bound.SetNodeIndex(nodeIndex);
			innerNodeIndex = structNodeIndex[i][j][k - 1];
			bound.SetInnerNodeIndex(innerNodeIndex);
			ghostNodeIndex = structNodeIndex[i][j][k + 1];
			bound.SetGhostNodeIndex(ghostNodeIndex);
			boundNorm = node[ghostNodeIndex].GetCoordinate() - node[nodeIndex].GetCoordinate();
			bound.SetNorm(boundNorm);
			boundMap->AddBoundary("wall", bound);
		}
	}
	gridList->AddGrid(testGrid);

}
