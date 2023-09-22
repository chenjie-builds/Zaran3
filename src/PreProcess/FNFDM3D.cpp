#include"FNFDM3D.h"
#include"Log.h"
#include<set>
#include<fstream>
using namespace zaran;
zaran::GridListFactoryFNFDM3D::GridListFactoryFNFDM3D()
{
	m_fileName = "sysu.dat";
}
void zaran::GridListFactoryFNFDM3D::Create(Ptr<GridList>& gridList)
{
	if (!gridList)
		gridList = std::make_shared<GridList>();
	ReadFile(gridList);
}

void zaran::GridListFactoryFNFDM3D::ReadFile(Ptr<GridList>& gridList)
{
	Ptr < Grid > grid = std::make_shared<Grid>();
	grid->SetDimension(Dimension::three);
	gridList->AddGrid(grid);
	auto& nodeTopo = grid->GetNodeTopo();
	std::ifstream fin(m_fileName);
	//读取所有节点坐标
	fin >> m_NodeNum;
	grid->SetTotalNodeNum(m_NodeNum);
	auto& nodeCoord = nodeTopo->GetCoordinate();
	nodeCoord.resize(m_NodeNum);
	for (size_t i = 0; i < m_NodeNum; i++)
	{
		auto& currentCoord = nodeCoord[i];
		fin >> currentCoord[0] >> currentCoord[1] >> currentCoord[2];
	}
	//读取所有内部节点邻居节点
	int innerNodeNum = 0;
	fin >> innerNodeNum;
	int innerNodeIndex;
	IArray neiborNodeIndex(6);
	auto& nodeType = nodeTopo->GetType();
	nodeType.resize(m_NodeNum);
	auto& temp_i = nodeTopo->GetTemplateI();
	auto& temp_j = nodeTopo->GetTemplateJ();
	auto& temp_k = nodeTopo->GetTemplateK();
	temp_i.resize(m_NodeNum);
	temp_j.resize(m_NodeNum);
	temp_k.resize(m_NodeNum);
	auto& nodeNeibor = nodeTopo->GetNeighborCloud();
	nodeNeibor.resize(m_NodeNum);
	for (size_t i = 0; i < innerNodeNum; i++)
	{
		fin >> innerNodeIndex;
		fin >> neiborNodeIndex[0] >> neiborNodeIndex[1] >> neiborNodeIndex[2] >> neiborNodeIndex[3] >> neiborNodeIndex[4] >> neiborNodeIndex[5];
		innerNodeIndex -= 1;
		neiborNodeIndex[0] -= 1;
		neiborNodeIndex[1] -= 1;
		neiborNodeIndex[2] -= 1;
		neiborNodeIndex[3] -= 1;
		neiborNodeIndex[4] -= 1;
		neiborNodeIndex[5] -= 1;
		nodeType[innerNodeIndex] = NodeType::inner;
		nodeNeibor[innerNodeIndex] = neiborNodeIndex;
		temp_i[innerNodeIndex] = IArray{ neiborNodeIndex[0],innerNodeIndex,neiborNodeIndex[1] };
		temp_j[innerNodeIndex] = IArray{ neiborNodeIndex[2],innerNodeIndex,neiborNodeIndex[3] };
		temp_k[innerNodeIndex] = IArray{ neiborNodeIndex[4],innerNodeIndex,neiborNodeIndex[5] };
	}
	//读取所有边界节点邻居节点
	auto& boundMap = grid->GetBoundaryMap();
	m_BoundNodeNum = 0;
	int nBound;
	fin >> nBound;
	m_BoundNodeNum += nBound;
	int tempIndex1, tempIndex2;
	int boundNodeIndex, connectNodeIndex;
	Boundary tempBound;
	for (size_t i = 0; i < nBound; i++)
	{
		fin >> boundNodeIndex >> tempIndex1;
		boundNodeIndex -= 1;
		boundMap->AddBoundary("inlet", Boundary{ boundNodeIndex,0,0,DVector3D{} });
		nodeType[boundNodeIndex] = NodeType::inlet;
	}
	fin >> nBound;
	m_BoundNodeNum += nBound;
	for (size_t i = 0; i < nBound; i++)
	{
		fin >> boundNodeIndex >> connectNodeIndex >> tempIndex1;
		boundNodeIndex -= 1;
		connectNodeIndex -= 1;
		boundMap->AddBoundary("outlet", Boundary{ boundNodeIndex,connectNodeIndex,0,DVector3D{} });
		nodeType[boundNodeIndex] = NodeType::outlet;
	}
	fin >> nBound;
	m_BoundNodeNum += nBound;
	for (size_t i = 0; i < nBound; i++)
	{
		fin >> boundNodeIndex >> connectNodeIndex >> tempIndex1;
		boundNodeIndex -= 1;
		connectNodeIndex -= 1;
		boundMap->AddBoundary("outlet", Boundary{ boundNodeIndex,connectNodeIndex,0,DVector3D{} });
		nodeType[boundNodeIndex] = NodeType::outlet;
	}
	fin >> nBound;
	m_BoundNodeNum += nBound;
	for (size_t i = 0; i < nBound; i++)
	{
		fin >> boundNodeIndex >> connectNodeIndex >> tempIndex1;
		boundNodeIndex -= 1;
		connectNodeIndex -= 1;
		boundMap->AddBoundary("outlet", Boundary{ boundNodeIndex,connectNodeIndex,0,DVector3D{} });
		nodeType[boundNodeIndex] = NodeType::outlet;
	}
	fin >> nBound;
	m_BoundNodeNum += nBound;
	for (size_t i = 0; i < nBound; i++)
	{
		fin >> boundNodeIndex >> connectNodeIndex >> tempIndex1;
		boundNodeIndex -= 1;
		connectNodeIndex -= 1;
		boundMap->AddBoundary("outlet", Boundary{ boundNodeIndex,connectNodeIndex,0,DVector3D{} });
		nodeType[boundNodeIndex] = NodeType::outlet;
	}
	fin >> nBound;
	m_BoundNodeNum += nBound;
	for (size_t i = 0; i < nBound; i++)
	{
		if (i == nBound - 1)
			i = i;
		fin >> boundNodeIndex >> connectNodeIndex >> tempIndex1;
		boundNodeIndex -= 1;
		connectNodeIndex -= 1;
		boundMap->AddBoundary("outlet", Boundary{ boundNodeIndex,connectNodeIndex,0,DVector3D{} });
		nodeType[boundNodeIndex] = NodeType::outlet;
	}
	fin >> nBound;
	m_BoundNodeNum += nBound;
	for (size_t i = 0; i < nBound; i++)
	{
		fin >> boundNodeIndex >> connectNodeIndex >> tempIndex1;
		boundNodeIndex -= 1;
		connectNodeIndex -= 1;
		DVector3D wallNorm = nodeCoord[connectNodeIndex] - nodeCoord[boundNodeIndex];
		boundMap->AddBoundary("slipWall", Boundary{ boundNodeIndex,connectNodeIndex,0,wallNorm });
		nodeType[boundNodeIndex] = NodeType::slipWall;
	}
	Array<std::set<int>> nodeNeiborSet(m_NodeNum);
	for (int iNode = 0; iNode < m_NodeNum; iNode++)
	{
		if (nodeType[iNode] != NodeType::inner)
			continue;
		auto& currentNeibor = nodeNeibor[iNode];
		auto& neiborSet = nodeNeiborSet[iNode];
		for (auto& iNeibor : currentNeibor)
		{
			neiborSet.insert(iNeibor);
		}
		for (auto& iNeibor : currentNeibor)
		{
			auto& neiborNeibor = nodeNeibor[iNeibor];
			for (auto& iNeiborNeibor : neiborNeibor)
			{
				neiborSet.insert(iNeiborNeibor);
			}
		}
		neiborSet.erase(iNode);
	}
	for (int iNode = 0; iNode < m_NodeNum; iNode++)
	{
		if (nodeType[iNode] != NodeType::inner)
			continue;
		auto& currentNeibor = nodeNeibor[iNode];
		auto& neiborSet = nodeNeiborSet[iNode];
		currentNeibor.resize(neiborSet.size());
		int i = 0;
		for (auto& iNeibor : neiborSet)
		{
			currentNeibor[i] = iNeibor;
			i++;
		}
	}
	fin.close();
	fin.open("cell.dat");
	auto& cellTopo = grid->GetCellTopo();
	int cellNum;
	fin >> cellNum;
	auto& cell_node = cellTopo->GetNodeIndex();
	cell_node.resize(cellNum);
	IArray cellNeiborNodeIndex(8);
	for (int iCell = 0; iCell < cellNum; iCell++)
	{
		fin >> cellNeiborNodeIndex[0] >> cellNeiborNodeIndex[1] >> cellNeiborNodeIndex[2] >> cellNeiborNodeIndex[3]
			>> cellNeiborNodeIndex[4] >> cellNeiborNodeIndex[5] >> cellNeiborNodeIndex[6] >> cellNeiborNodeIndex[7];
		cellNeiborNodeIndex[0] -= 1;
		cellNeiborNodeIndex[1] -= 1;
		cellNeiborNodeIndex[2] -= 1;
		cellNeiborNodeIndex[3] -= 1;
		cellNeiborNodeIndex[4] -= 1;
		cellNeiborNodeIndex[5] -= 1;
		cellNeiborNodeIndex[6] -= 1;
		cellNeiborNodeIndex[7] -= 1;
		cell_node[iCell] = cellNeiborNodeIndex;
	}
	fin.close();


}
