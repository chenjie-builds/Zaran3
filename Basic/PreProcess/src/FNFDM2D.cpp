#include "FNFDM2D.h"
#include"Log.h"
#include<fstream>
#include<set>
namespace zaran
{
	void GridListFactoryFNFDM2D::Create(Ptr<GridList>& gridList)
	{
		if (!gridList)
		{
			gridList = std::make_shared<GridList>();
		}
		Ptr<Grid> gridPtr = std::make_shared<Grid>();
		gridPtr->SetDimension(Dimension::two);
		gridPtr->SetName("FNFDM2D");
		gridList->AddGrid(gridPtr);
		std::string fileName = "grid.dat";
		std::ifstream fin;
		fin.open(fileName, std::ios_base::in);
		if (!fin.is_open())
		{
			ZaranLog::error("Grid File:{} is NOT exist! Please Check!", fileName);
			system("pause");
		}
		int nNode;
		fin >> nNode;
		gridPtr->SetTotalNodeNum(nNode);
		auto& nodeTopo = gridPtr->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& nodeType = nodeTopo->GetType();
		nodeType.resize(nNode);
		nodeCoord.resize(nNode);
		for (int iNode = 0; iNode < nNode; iNode++)
		{
			fin >> nodeCoord[iNode][0] >> nodeCoord[iNode][1];
			nodeType[iNode] = NodeType::undefined;

		}
		int nInnerNode;
		fin >> nInnerNode;
		IArray neiborNodeIndex(4);
		auto& temp_i = nodeTopo->GetTemplateI();
		auto& temp_j = nodeTopo->GetTemplateJ();
		temp_i.resize(nNode);
		temp_j.resize(nNode);
		auto& nodeNeibor = nodeTopo->GetNeighborCloud();
		nodeNeibor.resize(nNode);
		int nodeIndex;
		for (int iNode = 0; iNode < nInnerNode; iNode++)
		{
			fin >> nodeIndex;
			if (nodeIndex == 3389)
			{
				int a = 0;
			}
			nodeNeibor[nodeIndex].resize(4);
			fin >> neiborNodeIndex[0] >> neiborNodeIndex[1] >> neiborNodeIndex[2] >> neiborNodeIndex[3];
			nodeType[nodeIndex] = NodeType::inner;
			temp_i[nodeIndex] = IArray{ neiborNodeIndex[0],nodeIndex,neiborNodeIndex[1] };
			temp_j[nodeIndex] = IArray{ neiborNodeIndex[2],nodeIndex,neiborNodeIndex[3] };
			nodeNeibor[nodeIndex] = neiborNodeIndex;
		}
		auto& boundMap = gridPtr->GetBoundaryMap();
		int nBoundNode;
		fin >> nBoundNode;
		int boundNodeIndex, boundType, innerNodeIndex;
		DVector3D boundNorm;
		for (int iNode = 0; iNode < nBoundNode; iNode++)
		{
			fin >> boundNodeIndex >> boundType >> innerNodeIndex >> boundNorm[0] >> boundNorm[1] >> boundNorm[2];
			if (boundType == 0)
			{
				nodeType[boundNodeIndex] = NodeType::inlet;
				boundMap->AddBoundary("inlet", Boundary{ boundNodeIndex,innerNodeIndex,0,boundNorm });
			}
			else if (boundType == 1)
			{
				nodeType[boundNodeIndex] = NodeType::outlet;
				boundMap->AddBoundary("outlet", Boundary{ boundNodeIndex,innerNodeIndex,0,boundNorm });
			}
			else if (boundType == 2)
			{
				nodeType[boundNodeIndex] = NodeType::slipWall;
				boundMap->AddBoundary("slipWall", Boundary{ boundNodeIndex,innerNodeIndex,0,boundNorm });
			}
			else if (boundType == 3)
			{
				nodeType[boundNodeIndex] = NodeType::hole;
				boundMap->AddBoundary("hole", Boundary{ boundNodeIndex,innerNodeIndex,0,boundNorm });
			}
			else if (boundType == 4)
			{
				nodeType[boundNodeIndex] = NodeType::userDefined;
				boundMap->AddBoundary("userDefined", Boundary{ boundNodeIndex,innerNodeIndex,0,boundNorm });
			}
			else
			{
				ZaranLog::info("Boundary Type:{} is NOT exist! Please Check!", boundType);
			}
		}
		int nCell;
		fin >> nCell;
		auto& cellTopo = gridPtr->GetCellTopo();
		auto& cell_node = cellTopo->GetNodeIndex();
		cell_node.resize(nCell);
		IArray cellNeiborNodeIndex(4);
		for (int iCell = 0; iCell < nCell; iCell++)
		{
			fin >> cellNeiborNodeIndex[0] >> cellNeiborNodeIndex[1] >> cellNeiborNodeIndex[2] >> cellNeiborNodeIndex[3];
			cell_node[iCell] = cellNeiborNodeIndex;
		}
		fin.close();

		Array<std::set<int>> nodeNeiborSet(nNode);
		for (int iNode = 0; iNode < nNode; iNode++)
		{
			// if (nodeType[iNode] != NodeType::inner)
			// 	continue;
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
		for (int iNode = 0; iNode < nNode; iNode++)
		{
			// if (nodeType[iNode] != NodeType::inner)
			// 	continue;
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




	}

}