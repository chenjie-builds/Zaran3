#include"FNFDM3D.h"
#include"Log.h"
#include"MathBasic.h"
#include<set>
#include<fstream>
#include<vtkKdTreePointLocator.h>
#include<vtkNew.h>
#include<vtkPoints.h>
#include<vtkPolyData.h>
#include<vtkVertexGlyphFilter.h>
namespace zaran
{

	GridListFactoryFNFDM3D::GridListFactoryFNFDM3D()
	{
		m_node_file_name = "node.dat";
		m_ele_file_name = "cell.dat";
		m_bnd_file_name = "bnd.dat";
	}
	void GridListFactoryFNFDM3D::Create(Ptr<GridList>& gridList)
	{
		if (!gridList)
			gridList = std::make_shared<GridList>();
		ReadFile(gridList);
		//	SortNeiborNode(gridList);
	}

	void GridListFactoryFNFDM3D::ReadFile(Ptr<GridList>& gridList)
	{
		Ptr < Grid > grid = std::make_shared<Grid>();
		grid->SetDimension(Dimension::three);
		gridList->AddGrid(grid);
		auto& nodeTopo = grid->GetNodeTopo();
		std::ifstream fin(m_node_file_name);
		//读取所有节点坐标
		fin >> m_NodeNum;
		ZaranLog::info("Total node num:{}", m_NodeNum);
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
		ZaranLog::info("Inner node num:{}", innerNodeNum);
		int innerNodeIndex;
		IArray neibor_index(6);
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
		double delta=1e-5;
		for (size_t i = 0; i < innerNodeNum; i++)
		{
			fin >> innerNodeIndex;
			fin >> neibor_index[0] >> neibor_index[1] >> neibor_index[2] >> neibor_index[3] >> neibor_index[4] >> neibor_index[5];
			innerNodeIndex -= 1;
			neibor_index[0] -= 1;
			neibor_index[1] -= 1;
			neibor_index[2] -= 1;
			neibor_index[3] -= 1;
			neibor_index[4] -= 1;
			neibor_index[5] -= 1;
			if(innerNodeIndex==709761)
				ZaranLog::info("inner node index:{}", innerNodeIndex);
			nodeType[innerNodeIndex] = NodeType::inner;
			nodeNeibor[innerNodeIndex] = neibor_index;
			Array<DVector3D> vec(3);
			vec[0] = nodeCoord[neibor_index[1]] - nodeCoord[neibor_index[0]];
			vec[1] = nodeCoord[neibor_index[3]] - nodeCoord[neibor_index[2]];
			vec[2] = nodeCoord[neibor_index[5]] - nodeCoord[neibor_index[4]];
			double angle = AngleOfTwoArray3D(vec[0].data(), vec[1].data());
			//i,j方向平行
			if (abs(angle) < delta)
			{
				std::swap(neibor_index[1], neibor_index[3]);
			}
			else if (abs(angle - PI) < delta)
			{
				std::swap(neibor_index[1], neibor_index[2]);
			}

			angle = AngleOfTwoArray3D(vec[0].data(), vec[2].data());
			//i,k方向平行
			if (abs(angle) < delta)
			{
				std::swap(neibor_index[1], neibor_index[5]);
			}
			else if (abs(angle - PI) < delta)
			{
				std::swap(neibor_index[1], neibor_index[4]);
			}
			angle = AngleOfTwoArray3D(vec[1].data(), vec[2].data());
			//j,k方向平行
			if (abs(angle) < delta)
			{
				std::swap(neibor_index[3], neibor_index[5]);
			}
			else if (abs(angle - PI) < delta)
			{
				std::swap(neibor_index[3], neibor_index[4]);
			}
			//检查是否是右手坐标系
			vec[0] = nodeCoord[neibor_index[1]] - nodeCoord[neibor_index[0]];
			vec[1] = nodeCoord[neibor_index[3]] - nodeCoord[neibor_index[2]];
			vec[2] = nodeCoord[neibor_index[5]] - nodeCoord[neibor_index[4]];
			if (vec[0].cross(vec[1]).dot(vec[2]) < 0)
			{
				std::swap(neibor_index[4], neibor_index[5]);
			}
			temp_i[innerNodeIndex] = IArray{ neibor_index[0],innerNodeIndex,neibor_index[1] };
			temp_j[innerNodeIndex] = IArray{ neibor_index[2],innerNodeIndex,neibor_index[3] };
			temp_k[innerNodeIndex] = IArray{ neibor_index[4],innerNodeIndex,neibor_index[5] };


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
		ZaranLog::info("x- boundary node num:{}", nBound);
		for (size_t i = 0; i < nBound; i++)
		{
			fin >> boundNodeIndex >> connectNodeIndex >> tempIndex1;
			boundNodeIndex -= 1;
			connectNodeIndex -= 1;
			boundMap->AddBoundary("inlet", Boundary{ boundNodeIndex,connectNodeIndex,0,DVector3D{} });
			nodeType[boundNodeIndex] = NodeType::inlet;
		}
		fin >> nBound;
		m_BoundNodeNum += nBound;
		ZaranLog::info("x+ boundary node num:{}", nBound);
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
		ZaranLog::info("y- boundary node num:{}", nBound);
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
		ZaranLog::info("y+ boundary node num:{}", nBound);
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
		ZaranLog::info("z- boundary node num:{}", nBound);
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
		ZaranLog::info("z+ boundary node num:{}", nBound);
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
		ZaranLog::info("wall boundary node num:{}", nBound);
		for (size_t i = 0; i < nBound; i++)
		{
			fin >> boundNodeIndex >> connectNodeIndex >> tempIndex1;
			boundNodeIndex -= 1;
			connectNodeIndex -= 1;
			DVector3D wallNorm = nodeCoord[connectNodeIndex] - nodeCoord[boundNodeIndex];
			boundMap->AddBoundary("slipWall", Boundary{ boundNodeIndex,connectNodeIndex,0,wallNorm });
			nodeType[boundNodeIndex] = NodeType::slipWall;
		}
		ZaranLog::info("extend inner node neibor node");
		ExtendNeighborNode(gridList);

		//扩展内部节点邻居节点，用于计算梯度

		//查找邻居节点，看自身是否是其邻居，如不是，则加进去
		ZaranLog::info("Add self to neibor node's neibor node");
		bool find_current;
		for (int iNode = 0;iNode < m_NodeNum;iNode++)
		{
			if (nodeType[iNode] != NodeType::inner)
				continue;
			find_current = false;
			auto& currentNeibor = nodeNeibor[iNode];
			for (int iNeibor = 0;iNeibor < currentNeibor.size();iNeibor++)
			{
				auto& neighbor_neighbor = nodeNeibor[currentNeibor[iNeibor]];
				for (int jNeibor = 0;jNeibor < neighbor_neighbor.size();jNeibor++)
				{
					if (neighbor_neighbor[jNeibor] == iNode)
						find_current = true;
				}
				if (find_current == false)
				{
					neighbor_neighbor.push_back(iNode);
				}
			}
		}
		fin.close();
		ReadCellFile(gridList);
	}

	void GridListFactoryFNFDM3D::SortNeiborNode(Ptr<GridList>& gridList)
	{
		auto& grid = gridList->GetGrid(0);
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& nodeNeibor = nodeTopo->GetNeighborCloud();
		auto& nodeType = nodeTopo->GetType();
		auto& temp_i = nodeTopo->GetTemplateI();
		auto& temp_j = nodeTopo->GetTemplateJ();
		auto& temp_k = nodeTopo->GetTemplateK();
		struct node_pair
		{
			int node1, node2;
		};
		std::map<double, node_pair> node_pair_map;
		for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
		{
			node_pair_map.clear();
			if (nodeType[iNode] != NodeType::inner)
				continue;
			auto& currentNodeCoord = nodeCoord[iNode];
			auto& currentNeibor = nodeNeibor[iNode];
			//ZaranLog::info("node:{}\n", iNode);
			if (iNode == 15896)
			{
				ZaranLog::info("node before:{}", iNode);
				ZaranLog::info("{} {} {}", currentNodeCoord.x(), currentNodeCoord.y(), currentNodeCoord.z());
				for (auto& i : currentNeibor)
				{
					ZaranLog::info("{} {} {}", nodeCoord[i].x(), nodeCoord[i].y(), nodeCoord[i].z());
				}
			}
			//删除邻居节点中与当地节点距离小于小量的节点
			for (int i = 0; i < currentNeibor.size(); ++i)
			{
				if ((nodeCoord[currentNeibor[i]] - currentNodeCoord).norm() < 1e-6)
				{
					currentNeibor.erase(currentNeibor.begin() + i);
					--i;
				}

			}
			// 生成点对map
			// 以点对与iNode连线的夹角为key
			// 以点对为value
			// 点对为所有与iNode相邻的点对
			for (int i = 0; i < currentNeibor.size(); ++i)
			{
				for (int j = i + 1; j < currentNeibor.size(); ++j)
				{
					node_pair temp;
					temp.node1 = currentNeibor[i];
					temp.node2 = currentNeibor[j];
					DVector3D vec1 = nodeCoord[temp.node1] - nodeCoord[iNode];
					DVector3D vec2 = nodeCoord[temp.node2] - nodeCoord[iNode];
					double angle = AngleOfTwoArray3D(vec1.data(), vec2.data());
					node_pair_map[angle] = temp;
				}
			}
			//取出map中最后一个点对，即夹角最大的点对
			node_pair main_pair = node_pair_map.rbegin()->second;
			//从邻居节点中删除这个点对
			currentNeibor.erase(std::find(currentNeibor.begin(), currentNeibor.end(), main_pair.node1));
			currentNeibor.erase(std::find(currentNeibor.begin(), currentNeibor.end(), main_pair.node2));
			//主方向向量
			DVector3D main_vec = nodeCoord[main_pair.node1] - nodeCoord[main_pair.node2];

			//求出所有邻居节点在以主方向向量为法向量，经过当地节点的平面上的投影
			map<int, DVector3D> node_proj_map;
			for (int i = 0; i < currentNeibor.size(); ++i)
			{
				DVector3D vec = nodeCoord[currentNeibor[i]] - nodeCoord[iNode];
				vec = vec - vec.dot(main_vec) * main_vec / (main_vec.norm() * main_vec.norm());
				//如果投影向量的模长小于1e-6，删除该邻居节点
				if (vec.norm() < 1e-6)
				{
					currentNeibor.erase(currentNeibor.begin() + i);
					--i;
					continue;
				}
				node_proj_map[currentNeibor[i]] = vec;
			}
			//以第一个邻居节点投影向量为基准向量，求出基准向量以法向量为旋转轴旋转到其他投影向量的角度, 0~2pi
			map<double, int> node_angle_map;
			for (int i = 0; i < currentNeibor.size(); ++i)
			{
				if (i == 0)
				{
					node_angle_map[0] = currentNeibor[i];
					continue;
				}
				DVector3D vec = node_proj_map[currentNeibor[i]];
				double angle = AngleOfTwoArray3D(node_proj_map[currentNeibor[0]].data(), vec.data());
				if ((node_proj_map[currentNeibor[0]].cross(vec).dot(main_vec) < 0))
					angle = 2 * PI - angle;
				//如果map中已经有这个角度，比较两个角度对应的邻居节点的距离，删除距离大的邻居节点
				if (node_angle_map.find(angle) != node_angle_map.end())
				{
					if (node_proj_map[currentNeibor[i]].norm() > node_proj_map[node_angle_map[angle]].norm())
					{
						node_angle_map[angle] = currentNeibor[i];
						currentNeibor.erase(std::find(currentNeibor.begin(), currentNeibor.end(), node_angle_map[angle]));
					}
					else
					{
						currentNeibor.erase(std::find(currentNeibor.begin(), currentNeibor.end(), currentNeibor[i]));
					}
					--i;
					continue;
				}
				node_angle_map[angle] = currentNeibor[i];
			}
			//根据角度排序后的邻居节点
			currentNeibor.clear();
			for (auto& i : node_angle_map)
			{
				currentNeibor.push_back(i.second);
			}
			if (iNode == 15896)
			{
				ZaranLog::info("node after:{}", iNode);
				ZaranLog::info("{} {} {}", currentNodeCoord.x(), currentNodeCoord.y(), currentNodeCoord.z());
				for (auto& i : currentNeibor)
				{
					ZaranLog::info("{} {} {}", nodeCoord[i].x(), nodeCoord[i].y(), nodeCoord[i].z());
				}
			}
			//求出邻居节点与当地节点之间的距离
			map<int, double> node_dis_map;
			for (int i = 0; i < currentNeibor.size(); ++i)
			{
				node_dis_map[currentNeibor[i]] = node_proj_map[currentNeibor[i]].norm();
			}
			node_pair_map.clear();
			// 获取下一个点的lamda表达式
			auto get_next_node = [&](int iNode, IArray neiborNode)-> int
				{
					if (iNode == neiborNode.size() - 1)
						return 0;
					else
						return iNode + 1;
				};
			// 获取上一个点的lamda表达式
			auto get_last_node = [&](int iNode, IArray neiborNode)-> int
				{
					if (iNode == 0)
						return neiborNode.size() - 1;
					else
						return iNode - 1;
				};

			for (int i = 0; i < currentNeibor.size(); ++i)
			{
				node_pair temp;
				temp.node1 = currentNeibor[i];
				temp.node2 = currentNeibor[get_next_node(i, currentNeibor)];
				DVector3D vec1, vec2;
				vec1 = node_proj_map[temp.node1];
				vec2 = node_proj_map[temp.node2];
				double angle = AngleOfTwoArray3D(vec1.data(), vec2.data()) + GetRand(0.0, 1.0) * EPSILON_NUMBER;
				if (vec1.cross(vec2).dot(main_vec) < 0)
					angle = 2 * PI - angle;
				node_pair_map[angle] = temp;
			}


			if (node_pair_map.size() != currentNeibor.size())
				ZaranLog::info("node_pair_map.size()!=currentNeibor.size()");
			//取出map中第一个点对，即夹角最小的点对
			//删除这个点对中距离最大的点
			while (currentNeibor.size() > 4)
			{

				auto& temp_pair = node_pair_map.begin()->second;
				int remove_node, remove_index;
				node_pair temp;
				if (node_dis_map[temp_pair.node1] > node_dis_map[temp_pair.node2])
				{
					remove_node = temp_pair.node1;
					remove_index = std::find(currentNeibor.begin(), currentNeibor.end(), remove_node) - currentNeibor.begin();
					for (auto i : node_pair_map)
					{
						if (i.second.node1 == currentNeibor[get_last_node(remove_index, currentNeibor)] && i.second.node2 == currentNeibor[remove_index])
						{
							node_pair_map.erase(i.first);
							break;
						}
					}
				}
				else
				{
					remove_node = temp_pair.node2;
					remove_index = std::find(currentNeibor.begin(), currentNeibor.end(), remove_node) - currentNeibor.begin();
					for (auto i : node_pair_map)
					{
						if (i.second.node1 == currentNeibor[remove_index] && i.second.node2 == currentNeibor[get_next_node(remove_index, currentNeibor)])
						{
							node_pair_map.erase(i.first);
							break;
						}
					}
				}
				node_pair_map.erase(node_pair_map.begin());
				temp.node1 = currentNeibor[get_last_node(remove_index, currentNeibor)];
				temp.node2 = currentNeibor[get_next_node(remove_index, currentNeibor)];
				currentNeibor.erase(std::find(currentNeibor.begin(), currentNeibor.end(), remove_node));
				node_dis_map.erase(remove_node);
				DVector3D vec1, vec2;
				vec1 = node_proj_map[temp.node1];
				vec2 = node_proj_map[temp.node2];
				double angle = AngleOfTwoArray3D(vec1.data(), vec2.data()) + GetRand(0.0, 1.0) * EPSILON_NUMBER;
				if (vec1.cross(vec2).dot(main_vec) < 0)
					angle = 2 * PI - angle;
				node_pair_map[angle] = temp;
			}
			temp_k[iNode][0] = main_pair.node2;
			temp_k[iNode][1] = iNode;
			temp_k[iNode][2] = main_pair.node1;
			temp_i[iNode][0] = currentNeibor[2];
			temp_i[iNode][1] = iNode;
			temp_i[iNode][2] = currentNeibor[0];
			temp_j[iNode][0] = currentNeibor[3];
			temp_j[iNode][1] = iNode;
			temp_j[iNode][2] = currentNeibor[1];
			currentNeibor.push_back(main_pair.node1);
			currentNeibor.push_back(main_pair.node2);
			if (iNode == 15896)
			{
				ZaranLog::info("node after:{}", iNode);
				ZaranLog::info("{} {} {}", currentNodeCoord.x(), currentNodeCoord.y(), currentNodeCoord.z());
				for (auto& i : currentNeibor)
				{
					ZaranLog::info("{} {} {}", nodeCoord[i].x(), nodeCoord[i].y(), nodeCoord[i].z());
				}
			}
		}


	}

	void GridListFactoryFNFDM3D::ExtendNeighborNode(Ptr<GridList>& gridList)
	{
		// 构建节点KD树
		auto& grid = gridList->GetGrid(0);
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& nodeNeibor = nodeTopo->GetNeighborCloud();
		auto& nodeType = nodeTopo->GetType();
		vtkNew<vtkPoints> points;
		for (int i = 0; i < nodeCoord.size(); ++i)
		{
			points->InsertNextPoint(nodeCoord[i].data());
		}
		vtkNew<vtkPolyData> polydata;
		polydata->SetPoints(points);
		vtkNew<vtkVertexGlyphFilter> glyphFilter;
		glyphFilter->SetInputData(polydata);
		glyphFilter->Update();
		vtkNew<vtkKdTreePointLocator> kdTree;
		kdTree->SetDataSet(glyphFilter->GetOutput());
		kdTree->BuildLocator();
		int point_num = kdTree->GetDataSet()->GetNumberOfPoints();
		ZaranLog::info("point num:{}", point_num);
		int neibor_num_before = 0;
		int neibor_num_after = 0;
		// 扩展内部节点邻居节点，用于计算梯度
		for (int iNode = 0; iNode < m_NodeNum; iNode++)
		{
			if (nodeType[iNode] != NodeType::inner)
				continue;
			auto& currentNeibor = nodeNeibor[iNode];
			neibor_num_before = currentNeibor.size();

			double max_distance = 0;
			for (auto& iNeibor : currentNeibor)
			{
				max_distance = Max(max_distance, (nodeCoord[iNeibor] - nodeCoord[iNode]).norm());
			}
			//以当前节点为中心，以最大距离为半径，找到范围内的节点
			vtkNew<vtkIdList> result;
			kdTree->FindPointsWithinRadius(max_distance * 1.1, nodeCoord[iNode].data(), result);
			currentNeibor.clear();
			for (int i = 0; i < result->GetNumberOfIds(); ++i)
			{
				currentNeibor.push_back(result->GetId(i));
			}
			// 删除邻居节点中与当地节点距离小于小量的节点
			for (int i = 0; i < currentNeibor.size(); ++i)
			{
				if ((nodeCoord[currentNeibor[i]] - nodeCoord[iNode]).norm() < 1e-6)
				{
					currentNeibor.erase(currentNeibor.begin() + i);
					--i;
				}
			}
			neibor_num_after = currentNeibor.size();
			if (neibor_num_after < neibor_num_before)
				ZaranLog::info("node:{} neibor num before:{} neibor num after:{}", iNode, neibor_num_before, neibor_num_after);
		}

	}

	void GridListFactoryFNFDM3D::ReadCellFile(Ptr<GridList>& gridList)
	{
		auto& grid = gridList->GetGrid(0);
		std::ifstream fin;
		fin.open("cell.dat");
		auto& cellTopo = grid->GetCellTopo();
		int cellNum;
		fin >> cellNum;
		ZaranLog::info("Total cell num:{}", cellNum);
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

	void GridListFactoryFNFDM3D::ReadBoundFile(Ptr<GridList>& gridList)
	{
		auto& grid = gridList->GetGrid(0);
		//TODO:读取边界文件
	}

}