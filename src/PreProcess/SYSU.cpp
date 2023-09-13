#include"SYSU.h"
#include"Log.h"
#include<fstream>
using namespace zaran;
zaran::GridListFactorySYSU::GridListFactorySYSU()
{
	m_fileName = "sysu.dat";
}
void zaran::GridListFactorySYSU::Create(Ptr<GridList>& gridList)
{
	if (!gridList)
		gridList = std::make_shared<GridList>();
	ReadFile(gridList);
}

void zaran::GridListFactorySYSU::ReadFile(Ptr<GridList>& gridList)
{
	Ptr < Grid > grid = std::make_shared<Grid>();
	gridList->AddGrid(grid);
	auto& nodeTopoVec = grid->GetNodeTopoInfo();
	std::ifstream fin(m_fileName);
	//读取所有节点坐标
	fin >> m_NodeNum;
	grid->SetTotalNodeNumber(m_NodeNum);
	nodeTopoVec.resize(m_NodeNum);
	for (size_t i = 0; i < m_NodeNum; i++)
	{
		auto& currentNode = nodeTopoVec[i];
		auto& currentCoord = currentNode.GetCoordinate();
		fin >> currentCoord[0] >> currentCoord[1] >> currentCoord[2];
	}
	//读取所有内部节点邻居节点
	int innerNodeNum = 0;
	fin >> innerNodeNum;
	int innerNodeIndex;
	IArray neiborNodeIndex(6);
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
		auto& currentNode = nodeTopoVec[innerNodeIndex];
		currentNode.SetType(NodeType::inner);
		currentNode.SetNeighborCloud(neiborNodeIndex);
		currentNode.SetNeighborTemplateI(IArray{ neiborNodeIndex[0],innerNodeIndex,neiborNodeIndex[1] });
		currentNode.SetNeighborTemplateJ(IArray{ neiborNodeIndex[2],innerNodeIndex,neiborNodeIndex[3] });
		currentNode.SetNeighborTemplateK(IArray{ neiborNodeIndex[4],innerNodeIndex,neiborNodeIndex[5] });
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
		nodeTopoVec[boundNodeIndex].SetType(NodeType::inlet);
	}
	fin >> nBound;
	m_BoundNodeNum += nBound;
	for (size_t i = 0; i < nBound; i++)
	{
		fin >> boundNodeIndex >> connectNodeIndex >> tempIndex1;
		boundNodeIndex -= 1;
		connectNodeIndex -= 1;
		boundMap->AddBoundary("outlet", Boundary{ boundNodeIndex,connectNodeIndex,0,DVector3D{} });
		nodeTopoVec[boundNodeIndex].SetType(NodeType::outlet);
	}
	fin >> nBound;
	m_BoundNodeNum += nBound;
	for (size_t i = 0; i < nBound; i++)
	{
		fin >> boundNodeIndex >> connectNodeIndex >> tempIndex1;
		boundNodeIndex -= 1;
		connectNodeIndex -= 1;
		boundMap->AddBoundary("outlet", Boundary{ boundNodeIndex,connectNodeIndex,0,DVector3D{} });
		nodeTopoVec[boundNodeIndex].SetType(NodeType::outlet);
	}
	fin >> nBound;
	m_BoundNodeNum += nBound;
	for (size_t i = 0; i < nBound; i++)
	{
		fin >> boundNodeIndex >> connectNodeIndex >> tempIndex1;
		boundNodeIndex -= 1;
		connectNodeIndex -= 1;
		boundMap->AddBoundary("outlet", Boundary{ boundNodeIndex,connectNodeIndex,0,DVector3D{} });
		nodeTopoVec[boundNodeIndex].SetType(NodeType::outlet);
	}
	fin >> nBound;
	m_BoundNodeNum += nBound;
	for (size_t i = 0; i < nBound; i++)
	{
		fin >> boundNodeIndex >> connectNodeIndex >> tempIndex1;
		boundNodeIndex -= 1;
		connectNodeIndex -= 1;
		boundMap->AddBoundary("outlet", Boundary{ boundNodeIndex,connectNodeIndex,0,DVector3D{} });
		nodeTopoVec[boundNodeIndex].SetType(NodeType::outlet);
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
		nodeTopoVec[boundNodeIndex].SetType(NodeType::outlet);
	}
	fin >> nBound;
	m_BoundNodeNum += nBound;
	for (size_t i = 0; i < nBound; i++)
	{
		fin >> boundNodeIndex >> connectNodeIndex >> tempIndex1;
		boundNodeIndex -= 1;
		connectNodeIndex -= 1;
		DVector3D wallNorm = nodeTopoVec[connectNodeIndex].GetCoordinate() - nodeTopoVec[boundNodeIndex].GetCoordinate();
		boundMap->AddBoundary("slipWall", Boundary{ boundNodeIndex,connectNodeIndex,0,wallNorm });
		nodeTopoVec[boundNodeIndex].SetType(NodeType::slipWall);
	}

	//for (int iNode = 0; iNode < grid->GetTotalNodeNum(); iNode++)
	//{
	//	auto& currentNode = grid->GetNodeTopoInfo()[iNode];
	//	auto& coord = currentNode.GetCoordinate();
	//	if (abs(coord[0] + 0.61) < SMALL_NUMBER)
	//	{
	//		if (currentNode.GetType() == NodeType::inlet)
	//			continue;
	//		currentNode.SetType(NodeType::inlet);
	//		boundMap->AddBoundary("inlet", Boundary{ iNode ,0,0,DVector3D{} });
	//	}
	//}

	fin.close();
	fin.open("cell.dat");
	auto& cellTopoVec = grid->GetCellTopoInfo();
	int cellNum;
	fin >> cellNum;
	cellTopoVec.resize(cellNum);
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
		cellTopoVec[iCell].SetNode(cellNeiborNodeIndex);
	}

	fin.close();


}
