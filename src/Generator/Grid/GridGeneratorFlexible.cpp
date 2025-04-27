#include "GridGeneratorFlexible.h"
#include "Log.h"
#include "MathBasic.h"
#include <fstream>
namespace zaran
{

	GridGeneratorFlexible::GridGeneratorFlexible(const string& node_file_name, const string& ele_file_name, const string& bnd_file_name) :m_node_file_name(node_file_name), m_ele_file_name(ele_file_name), m_bnd_file_name(bnd_file_name)
	{
	}

	void GridGeneratorFlexible::CreateGrid(dynamic_array<shared_ptr<GridBase>>& grid_list)
	{
		ReadNodeFile();
		ReadCellFile();
		ReadBoundFile();
		//SortNeiborNode();
		// ExtendNeighborNode();
		CheckNode();
		CheckUnkownNode();
		SetBoundNeighbor();
		AddSelfToNeighbor();
		CheckNeighborNum();
		grid_list.resize(1);
		grid_list[0] = make_shared<GridFN>("FNFDM", 0, 3);
		auto grid = std::static_pointer_cast<GridFN>(grid_list[0]);
		ConvertToGrid(grid);
	}

	void GridGeneratorFlexible::ReadNodeFile()
	{
		std::ifstream fin(m_node_file_name);
		int node_num = 0;
		// read all node coord
		fin >> node_num;
		Log::info("Total node num:{}", node_num);
		m_node_coord.resize(node_num);
		m_node_type.resize(node_num);
		m_node_neibor.resize(node_num);
		for (size_t i = 0; i < node_num; i++)
		{
			m_node_coord[i].resize(3);
			fin >> m_node_coord[i][0] >> m_node_coord[i][1] >> m_node_coord[i][2];
		}
		// read inner node's neighbor node
		int innerNodeNum = 0;
		fin >> innerNodeNum;
		Log::info("Inner node num:{}", innerNodeNum);
		int inner_node_idx;
		for (size_t i = 0; i < node_num; i++)
		{
			m_node_type[i] = NodeType::undefined;
		}
		double delta = 1e-5;
		for (size_t i = 0; i < innerNodeNum; i++)
		{
			fin >> inner_node_idx;
			inner_node_idx -= 1;
			int neighbor_num = 6;
			auto& neighbor_index = m_node_neibor[inner_node_idx];
			neighbor_index.resize(neighbor_num);
			for (int j = 0; j < neighbor_num; j++)
			{
				fin >> neighbor_index[j];
				neighbor_index[j] -= 1;
			}
			m_node_type[inner_node_idx] = NodeType::inner;
		}
		// read bound node info
		map<int, string> bound_map = {
			{0, "none"},
			{1, "inlet"},
			{2, "outlet"},
			{3, "farfield"},
			{4, "wall"} };
		index_type bound_num;
		fin >> bound_num;
		m_bound_node.resize(bound_num);
		int bound_type;
		int ref_node_idx;
		Log::info("Boundary node num:{}", bound_num);
		for (size_t i = 0; i < bound_num; i++)
		{
			auto& bound_node = m_bound_node[i];
			fin >> bound_node.bound_index >> bound_type >> bound_node.ref_index;
			bound_node.bound_index -= 1;
			bound_node.ref_index -= 1;
			bound_node.type = bound_map[bound_type];
			fin >> bound_node.normal[0] >> bound_node.normal[1] >> bound_node.normal[2];
			if (bound_node.type == "inlet")
				m_node_type[bound_node.bound_index] = NodeType::inlet;
			else if (bound_node.type == "outlet")
				m_node_type[bound_node.bound_index] = NodeType::outlet;
			else if (bound_node.type == "farfield")
				m_node_type[bound_node.bound_index] = NodeType::farfield;
			else if (bound_node.type == "wall")
				m_node_type[bound_node.bound_index] = NodeType::wall;
		}
		fin.close();
	}


	void GridGeneratorFlexible::SortNeiborNode()
	{
		struct node_pair
		{
			int node1, node2;
		};
		int node_num = m_node_coord.size();
		std::map<double, node_pair> node_pair_map;
		for (int iNode = 0; iNode < node_num; ++iNode)
		{
			if(iNode==1380)
				int a = 0;
			node_pair_map.clear();
			if (m_node_type[iNode] != NodeType::inner)
				continue;
			auto& node_coord = m_node_coord[iNode];
			auto& neighbor = m_node_neibor[iNode];
			// 删除邻居节点中与当地节点距离小于小量的节�?
			for (int i = 0; i < neighbor.size(); ++i)
			{
				if (DistanceOfTwoPoints(node_coord.data(), m_node_coord[neighbor[i]].data()) < EPSILON_NUMBER)
				{
					neighbor.erase(neighbor.begin() + i);
					--i;
				}
			}
			// 生成点对map
			// 以点对与iNode连线的夹角为key
			// 以点对为value
			// 点对为所有与iNode相邻的点�?
			dynamic_array<double> vec1(3), vec2(3);
			for (int i = 0; i < neighbor.size(); ++i)
			{
				for (int j = i + 1; j < neighbor.size(); ++j)
				{
					node_pair temp;
					temp.node1 = neighbor[i];
					temp.node2 = neighbor[j];
					for (int k = 0; k < 3; ++k)
					{
						vec1[k] = m_node_coord[temp.node1][k] - node_coord[k];
						vec2[k] = m_node_coord[temp.node2][k] - node_coord[k];
					}
					double angle = AngleOfTwoArray3D(vec1.data(), vec2.data());
					node_pair_map[angle] = temp;
				}
			}
			// 取出map中最后一个点对，即夹角最大的点对
			node_pair main_pair = node_pair_map.rbegin()->second;
			// 从邻居节点中删除这个点对
			neighbor.erase(std::find(neighbor.begin(), neighbor.end(), main_pair.node1));
			neighbor.erase(std::find(neighbor.begin(), neighbor.end(), main_pair.node2));
			// 主方向向�?
			Eigen::Vector3d main_vec;
			for (int k = 0; k < 3; ++k)
			{
				main_vec[k] = m_node_coord[main_pair.node1][k] - m_node_coord[main_pair.node2][k];
			}
			main_vec.normalize();
			// 求出所有邻居节点在以主方向向量为法向量，经过当地节点的平面上的投影
			map<int, Eigen::Vector3d> node_proj_map;
			Eigen::Vector3d vec;
			for (int i = 0; i < neighbor.size(); ++i)
			{
				for (int k = 0; k < 3; ++k)
				{
					vec[k] = m_node_coord[neighbor[i]][k] - node_coord[k];
				}
				vec -= vec.dot(main_vec) * main_vec;
				double vec_norm = vec.norm();
				// 如果投影向量的模长小�?1e-6，删除该邻居节点
				if (vec_norm < EPSILON_NUMBER)
				{
					neighbor.erase(neighbor.begin() + i);
					--i;
					continue;
				}
				node_proj_map[neighbor[i]] = vec;
			}
			// 以第一个邻居节点投影向量为基准向量，求出基准向量以法向量为旋转轴旋转到其他投影向量的角�?, 0~2pi
			map<double, int> node_angle_map;
			for (int i = 0; i < neighbor.size(); ++i)
			{
				if (i == 0)
				{
					node_angle_map[0] = neighbor[i];
					continue;
				}
				Eigen::Vector3d vec = node_proj_map[neighbor[i]];
				double angle = AngleOfTwoArray3D(node_proj_map[neighbor[0]].data(), vec.data());
				if ((node_proj_map[neighbor[0]].cross(vec).dot(main_vec) < 0))
					angle = 2 * PI - angle;
				// 如果map中已经有这个角度，比较两个角度对应的邻居节点的距离，删除距离大的邻居节点
				if (node_angle_map.find(angle) != node_angle_map.end())
				{
					if (node_proj_map[neighbor[i]].norm() > node_proj_map[node_angle_map[angle]].norm())
					{
						node_angle_map[angle] = neighbor[i];
						neighbor.erase(std::find(neighbor.begin(), neighbor.end(), node_angle_map[angle]));
					}
					else
					{
						neighbor.erase(std::find(neighbor.begin(), neighbor.end(), neighbor[i]));
					}
					--i;
					continue;
				}
				node_angle_map[angle] = neighbor[i];
			}
			// 根据角度排序后的邻居节点
			neighbor.clear();
			for (auto& i : node_angle_map)
			{
				neighbor.push_back(i.second);
			}
			// 求出邻居节点与当地节点之间的距离
			map<int, double> node_dis_map;
			for (int i = 0; i < neighbor.size(); ++i)
			{
				node_dis_map[neighbor[i]] = node_proj_map[neighbor[i]].norm();
			}
			node_pair_map.clear();
			// 获取下一个点的lamda表达�?
			auto get_next_node = [&](index_type iNode, dynamic_array<index_type> neiborNode) -> int
				{
					if (iNode == neiborNode.size() - 1)
						return 0;
					else
						return iNode + 1;
				};
			// 获取上一个点的lamda表达�?
			auto get_last_node = [&](index_type iNode, dynamic_array<index_type> neiborNode) -> int
				{
					if (iNode == 0)
						return neiborNode.size() - 1;
					else
						return iNode - 1;
				};

			for (size_t i = 0; i < neighbor.size(); ++i)
			{
				node_pair temp;
				temp.node1 = neighbor[i];
				temp.node2 = neighbor[get_next_node(i, neighbor)];
				Eigen::Vector3d vec1, vec2;
				vec1 = node_proj_map[temp.node1];
				vec2 = node_proj_map[temp.node2];
				double angle = AngleOfTwoArray3D(vec1.data(), vec2.data()) + GetRand(0.0, 1.0) * EPSILON_NUMBER;
				if (vec1.cross(vec2).dot(main_vec) < 0)
					angle = 2 * PI - angle;
				node_pair_map[angle] = temp;
			}

			if (node_pair_map.size() != neighbor.size())
				Log::info("node_pair_map.size()!=currentNeibor.size()");
			// 取出map中第一个点对，即夹角最小的点对
			// 删除这个点对中距离最大的�?
			while (neighbor.size() > 4)
			{

				auto& temp_pair = node_pair_map.begin()->second;
				int remove_node, remove_index;
				node_pair temp;
				if (node_dis_map[temp_pair.node1] > node_dis_map[temp_pair.node2])
				{
					remove_node = temp_pair.node1;
					remove_index = std::find(neighbor.begin(), neighbor.end(), remove_node) - neighbor.begin();
					for (auto& i : node_pair_map)
					{
						if (i.second.node1 == neighbor[get_last_node(remove_index, neighbor)] &&
							i.second.node2 == neighbor[remove_index])
						{
							node_pair_map.erase(i.first);
							break;
						}
					}
				}
				else
				{
					remove_node = temp_pair.node2;
					remove_index = std::find(neighbor.begin(), neighbor.end(), remove_node) - neighbor.begin();
					for (auto& i : node_pair_map)
					{
						if (i.second.node1 == neighbor[remove_index] &&
							i.second.node2 == neighbor[get_next_node(remove_index, neighbor)])
						{
							node_pair_map.erase(i.first);
							break;
						}
					}
				}
				node_pair_map.erase(node_pair_map.begin());
				temp.node1 = neighbor[get_last_node(remove_index, neighbor)];
				temp.node2 = neighbor[get_next_node(remove_index, neighbor)];
				neighbor.erase(std::find(neighbor.begin(), neighbor.end(), remove_node));
				node_dis_map.erase(remove_node);
				Eigen::Vector3d vec1, vec2;
				vec1 = node_proj_map[temp.node1];
				vec2 = node_proj_map[temp.node2];
				double angle = AngleOfTwoArray3D(vec1.data(), vec2.data()) + GetRand(0.0, 1.0) * EPSILON_NUMBER;
				if (vec1.cross(vec2).dot(main_vec) < 0)
					angle = 2 * PI - angle;
				node_pair_map[angle] = temp;
			}
			//! 前四个点没有排序
			neighbor.push_back(main_pair.node1);
			neighbor.push_back(main_pair.node2);
		}
	}

	void GridGeneratorFlexible::ReadCellFile()
	{
		std::ifstream fin(m_ele_file_name);
		index_type cellNum;
		fin >> cellNum;
		m_cell_node.resize(cellNum);
		for (int iCell = 0; iCell < cellNum; iCell++)
		{
			int cell_node_num = 8;
			m_cell_node[iCell].resize(cell_node_num);
			for (int i = 0; i < cell_node_num; i++)
			{
				fin >> m_cell_node[iCell][i];
				m_cell_node[iCell][i] -= 1;
			}
		}
		fin.close();

	}

	void GridGeneratorFlexible::ReadBoundFile()
	{

	}

	void GridGeneratorFlexible::CheckNode()
	{
		double delta = 1e-5;
		double min_angle = LARGE_NUMBER;
		int min_angle_index;
		int min_angle_neibor_index1, min_angle_neibor_index2;
		int node_num = m_node_coord.size();
		// 各个方向的最小夹角, i,j,k方向
		double min_angle_i = LARGE_NUMBER, min_angle_j = LARGE_NUMBER, min_angle_k = LARGE_NUMBER;
		// 各个方向最小夹角的节点索引
		int min_angle_index_i, min_angle_index_j, min_angle_index_k;
		for (size_t iNode = 0; iNode < node_num; iNode++)
		{
			if (m_node_type[iNode] != NodeType::inner)
				continue;
			int neighbor_num = m_node_neibor[iNode].size();
			auto& neighbor_index = m_node_neibor[iNode];
			dynamic_array<Eigen::Vector3d> vec(3);
			for (int i = 0; i < 3; i++)
			{
				vec[0][i] = m_node_coord[neighbor_index[1]][i] - m_node_coord[neighbor_index[0]][i];
				vec[1][i] = m_node_coord[neighbor_index[3]][i] - m_node_coord[neighbor_index[2]][i];
				vec[2][i] = m_node_coord[neighbor_index[5]][i] - m_node_coord[neighbor_index[4]][i];
			}
			// 检查是否是右手坐标系
			if (vec[0].cross(vec[1]).dot(vec[2]) < 0)
			{
				//Log::info("Invalid jacobi matrix, node index:{}, neighbor index:{}, {}, {}, {}, {}, {}", iNode, neighbor_index[0], neighbor_index[1], neighbor_index[2], neighbor_index[3], neighbor_index[4], neighbor_index[5]);
				std::swap(neighbor_index[4], neighbor_index[5]);
			}
			double angle = AngleOfTwoArray3D(vec[0].data(), vec[1].data());
			// i,j方向平行
			if (abs(angle) < delta)
			{
				std::swap(neighbor_index[1], neighbor_index[3]);
			}
			else if (abs(angle - PI) < delta)
			{
				std::swap(neighbor_index[1], neighbor_index[2]);
			}

			angle = AngleOfTwoArray3D(vec[0].data(), vec[2].data());
			// i,k方向平行
			if (abs(angle) < delta)
			{
				std::swap(neighbor_index[1], neighbor_index[5]);
			}
			else if (abs(angle - PI) < delta)
			{
				std::swap(neighbor_index[1], neighbor_index[4]);
			}
			angle = AngleOfTwoArray3D(vec[1].data(), vec[2].data());
			// j,k方向平行
			if (abs(angle) < delta)
			{
				std::swap(neighbor_index[3], neighbor_index[5]);
			}
			else if (abs(angle - PI) < delta)
			{
				std::swap(neighbor_index[3], neighbor_index[4]);
			}
			// 检查是否是右手坐标系
			for (int i = 0; i < 3; i++)
			{
				vec[0][i] = m_node_coord[neighbor_index[1]][i] - m_node_coord[neighbor_index[0]][i];
				vec[1][i] = m_node_coord[neighbor_index[3]][i] - m_node_coord[neighbor_index[2]][i];
				vec[2][i] = m_node_coord[neighbor_index[5]][i] - m_node_coord[neighbor_index[4]][i];
			}
			if (vec[0].cross(vec[1]).dot(vec[2]) < 0)
			{
				std::swap(neighbor_index[4], neighbor_index[5]);
			}
			// 计算坐标轴之间的最小夹角
			double temp_min_angle = LARGE_NUMBER;
			int temp_index1, temp_index2;
			for (int iNeigh = 0; iNeigh < 6; ++iNeigh)
			{
				for (int k = 0; k < 3; ++k)
				{
					vec[0][k] = m_node_coord[neighbor_index[iNeigh]][k] - m_node_coord[iNode][k];
				}
				for (int jNeigh = iNeigh + 1; jNeigh < 6; ++jNeigh)
				{
					for (int k = 0; k < 3; ++k)
					{
						vec[1][k] = m_node_coord[neighbor_index[jNeigh]][k] - m_node_coord[iNode][k];
					}
					double angle = AngleOfTwoArray3D(vec[0].data(), vec[1].data());
					if (angle < temp_min_angle)
					{
						temp_min_angle = angle;
						temp_index1 = neighbor_index[iNeigh];
						temp_index2 = neighbor_index[jNeigh];
					}
					if (iNeigh == 0 && jNeigh == 1)
					{
						if (angle < min_angle_i)
						{
							min_angle_i = angle;
							min_angle_index_i = iNode;
						}
					}
					if (iNeigh == 2 && jNeigh == 3)
					{
						if (angle < min_angle_j)
						{
							min_angle_j = angle;
							min_angle_index_j = iNode;
						}
					}
					if (iNeigh == 4 && jNeigh == 5)
					{
						if (angle < min_angle_k)
						{
							min_angle_k = angle;
							min_angle_index_k = iNode;
						}
					}
				}
			}
			if (temp_min_angle < min_angle)
			{
				min_angle = temp_min_angle;
				min_angle_index = iNode;
				min_angle_neibor_index1 = temp_index1;
				min_angle_neibor_index2 = temp_index2;
			}
		}
		Log::info("Min axis Angle: {}, Node index: {}", min_angle, min_angle_index);
		Log::info("-------- Coord: {}, {}, {}", m_node_coord[min_angle_index][0], m_node_coord[min_angle_index][1],
			m_node_coord[min_angle_index][2]);
		Log::info("-------- Neighbor: {}, {}, {}, {}, {}, {}", m_node_neibor[min_angle_index][0],
			m_node_neibor[min_angle_index][1], m_node_neibor[min_angle_index][2], m_node_neibor[min_angle_index][3],
			m_node_neibor[min_angle_index][4], m_node_neibor[min_angle_index][5]);
		Log::info("-------- Axis Node 1: {}, Node 2: {}", min_angle_neibor_index1, min_angle_neibor_index2);
		Log::info("Min angle i: {}, Node index: {}", min_angle_i, min_angle_index_i);
		Log::info("-------- Coord: {}, {}, {}", m_node_coord[min_angle_index_i][0], m_node_coord[min_angle_index_i][1],
			m_node_coord[min_angle_index_i][2]);
		Log::info("-------- Neighbor: {}, {}, {}, {}, {}, {}", m_node_neibor[min_angle_index_i][0], m_node_neibor[min_angle_index_i][1], m_node_neibor[min_angle_index_i][2], m_node_neibor[min_angle_index_i][3], m_node_neibor[min_angle_index_i][4], m_node_neibor[min_angle_index_i][5]);
		Log::info("Min angle j: {}, Node index: {}", min_angle_j, min_angle_index_j);
		Log::info("-------- Coord: {}, {}, {}", m_node_coord[min_angle_index_j][0], m_node_coord[min_angle_index_j][1],
			m_node_coord[min_angle_index_j][2]);
		Log::info("-------- Neighbor: {}, {}, {}, {}, {}, {}", m_node_neibor[min_angle_index_j][0], m_node_neibor[min_angle_index_j][1], m_node_neibor[min_angle_index_j][2], m_node_neibor[min_angle_index_j][3], m_node_neibor[min_angle_index_j][4], m_node_neibor[min_angle_index_j][5]);
		Log::info("Min angle k: {}, Node index: {}", min_angle_k, min_angle_index_k);
		Log::info("-------- Coord: {}, {}, {}", m_node_coord[min_angle_index_k][0], m_node_coord[min_angle_index_k][1],
			m_node_coord[min_angle_index_k][2]);
		Log::info("-------- Neighbor: {}, {}, {}, {}, {}, {}", m_node_neibor[min_angle_index_k][0], m_node_neibor[min_angle_index_k][1], m_node_neibor[min_angle_index_k][2], m_node_neibor[min_angle_index_k][3], m_node_neibor[min_angle_index_k][4], m_node_neibor[min_angle_index_k][5]);

	}

	void GridGeneratorFlexible::CheckUnkownNode()
	{
		int node_num = m_node_coord.size();
		for (size_t i = 0; i < node_num; i++)
		{
			if (m_node_type[i] == NodeType::undefined)
			{
				Log::info("node:{} type is undefined", i);
			}
		}
		Log::info("Check undefined node done");

	}

	void GridGeneratorFlexible::CheckNeighborNum()
	{
		int node_num = m_node_coord.size();
		int min_neibor_num = 1E5;
		int max_neibor_num = 0;
		int min_neibor_index = 0;
		int max_neibor_index = 0;
		for (int iNode = 0; iNode < node_num; iNode++)
		{
			auto& nodeNeibor = m_node_neibor[iNode];
			if (min_neibor_num > nodeNeibor.size())
			{
				min_neibor_num = nodeNeibor.size();
				min_neibor_index = iNode;
			}
			if (max_neibor_num < nodeNeibor.size())
			{
				max_neibor_num = nodeNeibor.size();
				max_neibor_index = iNode;
			}
		}
		Log::info("min neibor num:{} max neibor num:{}", min_neibor_num, max_neibor_num);
		Log::info("min neibor index:{} max neibor index:{}", min_neibor_index, max_neibor_index);
	}

	void GridGeneratorFlexible::AddSelfToNeighbor()
	{
		int node_num = m_node_coord.size();
		for (int iNode = 0; iNode < node_num; iNode++)
		{
			bool find_current = false;
			auto& current_neighbor = m_node_neibor[iNode];
			for (int iNeibor = 0; iNeibor < current_neighbor.size(); iNeibor++)
			{
				auto& neighbors_neighbor = m_node_neibor[current_neighbor[iNeibor]];
				if (std::find(neighbors_neighbor.begin(), neighbors_neighbor.end(), iNode) == neighbors_neighbor.end())
				{
					neighbors_neighbor.push_back(iNode);
				}
			}
		}
		Log::info("Add self to neibor node's neibor node done");
	}

	void GridGeneratorFlexible::SetBoundNeighbor()
	{

	}

	void GridGeneratorFlexible::ConvertToGrid(shared_ptr<GridFN>grid)
	{
		index_type node_num = m_node_coord.size();
		dynamic_array<index_type> neighbor_node_num(node_num);
		dynamic_array<index_type> neighbor_face_num(node_num);
		dynamic_array<index_type> neighbor_cell_num(node_num);
		for (index_type iNode = 0; iNode < node_num; iNode++)
		{
			neighbor_node_num[iNode] = m_node_neibor[iNode].size();
			neighbor_face_num[iNode] = 0;
			neighbor_cell_num[iNode] = 0;
		}
		dynamic_array<dynamic_array<index_type>> neighbor_face(node_num, dynamic_array<index_type>(0));
		dynamic_array<dynamic_array<index_type>> neighbor_cell(node_num, dynamic_array<index_type>(0));
		auto& node = grid->GetNode();
		node.SetNodeNum(node_num);
		node.SetNeighborNode(m_node_neibor);
		node.SetCoord(m_node_coord);
		node.SetType(m_node_type);
		node.SetNeighborFace(neighbor_face);
		node.SetNeighborCell(neighbor_cell);

		auto& cell = grid->GetCell();
		cell.SetCellNum(m_cell_node.size());
		cell.SetNode(m_cell_node);
		dynamic_array<dynamic_array<index_type>> cell_node_face(m_cell_node.size(), dynamic_array<index_type>(0));
		cell.SetFace(cell_node_face);
		dynamic_array<double> center_coord(3);
		for (index_type iCell = 0; iCell < m_cell_node.size(); iCell++)
		{
			center_coord[0] = center_coord[1] = center_coord[2] = 0;
			for (index_type iNode = 0; iNode < m_cell_node[iCell].size(); iNode++)
			{
				for (index_type i = 0; i < 3; i++)
				{
					center_coord[i] += m_node_coord[m_cell_node[iCell][iNode]][i];
				}
			}
			for (index_type i = 0; i < 3; i++)
			{
				center_coord[i] /= m_cell_node[iCell].size();
			}
			cell.SetCenterCoord(iCell, center_coord.data());
		}

		FaceFN& face = grid->GetFace();

		dynamic_array<index_type> face_node_num(m_bound_face.size());
		for (index_type iFace = 0; iFace < m_bound_face.size(); iFace++)
		{
			face_node_num[iFace] = m_bound_face[iFace].face_node.size();
		}
		face.Allocate(m_bound_face.size(), face_node_num.data());
		for (index_type iFace = 0; iFace < m_bound_face.size(); iFace++)
		{
			face.SetFace2Node(iFace, m_bound_face[iFace].face_node.data(), m_bound_face[iFace].face_node.size());
			face.SetNormal(iFace, m_bound_face[iFace].normal.data());
			face.SetArea(iFace, m_bound_face[iFace].area);
		}

		BoundManagerFN& boundary_map = grid->GetBoundaryMap();
		for (int iFace = 0; iFace < m_bound_node.size(); iFace++)
		{
			boundary_map.AddBoundary(m_bound_node[iFace].type, BoundFN(m_bound_node[iFace].bound_index, m_bound_node[iFace].ref_index, 0, m_bound_node[iFace].normal));
		}
	}

} // namespace zaran