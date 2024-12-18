#include "GridGeneratorFlexibleSYSU.h"
#include "BasicType.h"
#include "Log.h"
#include "MathBasic.h"
#include <fstream>
#include <set>
#include <vtkKdTreePointLocator.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkVertexGlyphFilter.h>
#include <algorithm>
namespace zaran
{

	GridBuilderSYSU_FN::GridBuilderSYSU_FN(const string& node_file_name, const string& ele_file_name, const string& bnd_file_name) : m_node_file_name(node_file_name), m_ele_file_name(ele_file_name), m_bnd_file_name(bnd_file_name)
	{
	}

	void GridBuilderSYSU_FN::CreateGrid(dynamic_array<shared_ptr<GridBase>>& grid_list)
	{
		ReadNodeFile();
		ReadCellFile();
		ReadBoundFile();
		// SortNeiborNode();
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

	void GridBuilderSYSU_FN::ReadNodeFile()
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
		int innerNodeIndex;
		for (size_t i = 0; i < node_num; i++)
		{
			m_node_type[i] = NodeType::undefined;
		}
		double delta = 1e-5;
		for (size_t i = 0; i < innerNodeNum; i++)
		{
			fin >> innerNodeIndex;
			innerNodeIndex -= 1;
			int neighbor_num = 6;
			auto& neighbor_index = m_node_neibor[innerNodeIndex];
			neighbor_index.resize(neighbor_num);
			for (int j = 0; j < neighbor_num; j++)
			{
				fin >> neighbor_index[j];
				neighbor_index[j] -= 1;
			}
			m_node_type[innerNodeIndex] = NodeType::inner;
		}
		// read bound node info
		int total_bound_node_num = 0;
		int nBound;
		fin >> nBound;
		total_bound_node_num += nBound;
		m_bound_node.resize(total_bound_node_num);
		int tempIndex1, tempIndex2;
		int boundNodeIndex, connectNodeIndex;
		Log::info("x- boundary node num:{}", nBound);
		for (size_t i = 0; i < nBound; i++)
		{
			auto& bound_node = m_bound_node[i];
			fin >> bound_node.bound_index >> bound_node.ref_index >> tempIndex1;
			bound_node.bound_index -= 1;
			bound_node.ref_index -= 1;
			bound_node.type = "inlet";
			double mod_norm = 0;
			for (int j = 0; j < 3; j++)
			{
				bound_node.normal[j] = m_node_coord[bound_node.bound_index][j] - m_node_coord[bound_node.ref_index][j];
				mod_norm += bound_node.normal[j] * bound_node.normal[j];
			}
			mod_norm = sqrt(mod_norm);
			for (int j = 0; j < 3; j++)
			{
				bound_node.normal[j] /= mod_norm;
			}
			m_node_type[bound_node.bound_index] = NodeType::inlet;
		}
		fin >> nBound;
		total_bound_node_num += nBound;
		m_bound_node.resize(total_bound_node_num);
		Log::info("x+ boundary node num:{}", nBound);
		for (size_t i = 0; i < nBound; i++)
		{
			auto& bound_node = m_bound_node[i + total_bound_node_num - nBound];
			fin >> bound_node.bound_index >> bound_node.ref_index >> tempIndex1;
			bound_node.bound_index -= 1;
			bound_node.ref_index -= 1;
			bound_node.type = "outlet";
			double mod_norm = 0;
			for (int j = 0; j < 3; j++)
			{
				bound_node.normal[j] = m_node_coord[bound_node.bound_index][j] - m_node_coord[bound_node.ref_index][j];
				mod_norm += bound_node.normal[j] * bound_node.normal[j];
			}
			mod_norm = sqrt(mod_norm);
			for (int j = 0; j < 3; j++)
			{
				bound_node.normal[j] /= mod_norm;
			}
			m_node_type[bound_node.bound_index] = NodeType::outlet;
		}
		fin >> nBound;
		total_bound_node_num += nBound;
		m_bound_node.resize(total_bound_node_num);
		Log::info("y- boundary node num:{}", nBound);
		for (size_t i = 0; i < nBound; i++)
		{
			auto& bound_node = m_bound_node[i + total_bound_node_num - nBound];
			fin >> bound_node.bound_index >> bound_node.ref_index >> tempIndex1;
			bound_node.bound_index -= 1;
			bound_node.ref_index -= 1;
			bound_node.type = "outlet";
			double mod_norm = 0;
			for (int j = 0; j < 3; j++)
			{
				bound_node.normal[j] = m_node_coord[bound_node.bound_index][j] - m_node_coord[bound_node.ref_index][j];
				mod_norm += bound_node.normal[j] * bound_node.normal[j];
			}
			mod_norm = sqrt(mod_norm);
			for (int j = 0; j < 3; j++)
			{
				bound_node.normal[j] /= mod_norm;
			}
			m_node_type[bound_node.bound_index] = NodeType::outlet;
		}
		fin >> nBound;
		total_bound_node_num += nBound;
		m_bound_node.resize(total_bound_node_num);
		Log::info("y+ boundary node num:{}", nBound);
		for (size_t i = 0; i < nBound; i++)
		{
			auto& bound_node = m_bound_node[i + total_bound_node_num - nBound];
			fin >> bound_node.bound_index >> bound_node.ref_index >> tempIndex1;
			bound_node.bound_index -= 1;
			bound_node.ref_index -= 1;
			bound_node.type = "outlet";
			double mod_norm = 0;
			for (int j = 0; j < 3; j++)
			{
				bound_node.normal[j] = m_node_coord[bound_node.bound_index][j] - m_node_coord[bound_node.ref_index][j];
				mod_norm += bound_node.normal[j] * bound_node.normal[j];
			}
			mod_norm = sqrt(mod_norm);
			for (int j = 0; j < 3; j++)
			{
				bound_node.normal[j] /= mod_norm;
			}
			m_node_type[bound_node.bound_index] = NodeType::outlet;
		}
		fin >> nBound;
		total_bound_node_num += nBound;
		m_bound_node.resize(total_bound_node_num);
		Log::info("z- boundary node num:{}", nBound);
		for (size_t i = 0; i < nBound; i++)
		{
			auto& bound_node = m_bound_node[i + total_bound_node_num - nBound];
			fin >> bound_node.bound_index >> bound_node.ref_index >> tempIndex1;
			bound_node.bound_index -= 1;
			bound_node.ref_index -= 1;
			bound_node.type = "outlet";
			double mod_norm = 0;
			for (int j = 0; j < 3; j++)
			{
				bound_node.normal[j] = m_node_coord[bound_node.bound_index][j] - m_node_coord[bound_node.ref_index][j];
				mod_norm += bound_node.normal[j] * bound_node.normal[j];
			}
			mod_norm = sqrt(mod_norm);
			for (int j = 0; j < 3; j++)
			{
				bound_node.normal[j] /= mod_norm;
			}
			m_node_type[bound_node.bound_index] = NodeType::outlet;
		}
		fin >> nBound;
		total_bound_node_num += nBound;
		m_bound_node.resize(total_bound_node_num);
		Log::info("z+ boundary node num:{}", nBound);
		for (size_t i = 0; i < nBound; i++)
		{
			auto& bound_node = m_bound_node[i + total_bound_node_num - nBound];
			fin >> bound_node.bound_index >> bound_node.ref_index >> tempIndex1;
			bound_node.bound_index -= 1;
			bound_node.ref_index -= 1;
			bound_node.type = "outlet";
			double mod_norm = 0;
			for (int j = 0; j < 3; j++)
			{
				bound_node.normal[j] = m_node_coord[bound_node.bound_index][j] - m_node_coord[bound_node.ref_index][j];
				mod_norm += bound_node.normal[j] * bound_node.normal[j];
			}
			mod_norm = sqrt(mod_norm);
			for (int j = 0; j < 3; j++)
			{
				bound_node.normal[j] /= mod_norm;
			}
			m_node_type[bound_node.bound_index] = NodeType::outlet;
		}
		fin >> nBound;
		total_bound_node_num += nBound;
		m_bound_node.resize(total_bound_node_num);
		Log::info("wall boundary node num:{}", nBound);
		for (size_t i = 0; i < nBound; i++)
		{
			auto& bound_node = m_bound_node[i + total_bound_node_num - nBound];
			fin >> bound_node.bound_index >> bound_node.ref_index >> tempIndex1;
			bound_node.bound_index -= 1;
			bound_node.ref_index -= 1;
			bound_node.type = "wall";
			double mod_norm = 0;
			for (int j = 0; j < 3; j++)
			{
				bound_node.normal[j] = m_node_coord[bound_node.bound_index][j] - m_node_coord[bound_node.ref_index][j];
				mod_norm += bound_node.normal[j] * bound_node.normal[j];
			}
			mod_norm = sqrt(mod_norm);
			for (int j = 0; j < 3; j++)
			{
				bound_node.normal[j] /= mod_norm;
			}
			m_node_type[bound_node.bound_index] = NodeType::wall;
		}
		fin.close();
	}

	void GridBuilderSYSU_FN::SortNeiborNode()
	{
		struct node_pair
		{
			int node1, node2;
		};
		int node_num = m_node_coord.size();
		std::map<double, node_pair> node_pair_map;
		for (int iNode = 0; iNode < node_num; ++iNode)
		{
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

	void GridBuilderSYSU_FN::ExtendNeighborNode()
	{
		// 构建节点KD?
		// 初始化vtk?
		vtkNew<vtkPoints> points;
		for (int i = 0; i < m_node_coord.size(); ++i)
		{
			points->InsertNextPoint(m_node_coord[i].data());
		}
		// 把vtk点转换为vtkploydata
		vtkNew<vtkPolyData> polydata;
		polydata->SetPoints(points);
		// 将vtkploydata转换为vtk点定位器
		vtkNew<vtkVertexGlyphFilter> glyphFilter;
		glyphFilter->SetInputData(polydata);
		glyphFilter->Update();
		// 构建KD�?
		vtkNew<vtkKdTreePointLocator> kdTree;
		kdTree->SetDataSet(glyphFilter->GetOutput());
		kdTree->BuildLocator();
		int point_num = kdTree->GetDataSet()->GetNumberOfPoints();
		int neibor_num_before = 0;
		int neibor_num_after = 0;

		int node_num = m_node_coord.size();
		// 扩展内部节点邻居节点，用于计算梯度
		std::set<int> neibor_set;
		for (int iNode = 0; iNode < node_num; iNode++)
		{
			if (m_node_type[iNode] != NodeType::inner)
				continue;
			neibor_set.clear();
			auto& nodeNeibor = m_node_neibor[iNode];
			neibor_num_before = nodeNeibor.size();
			double max_distance = 0;
			for (auto& iNeibor : nodeNeibor)
			{
				max_distance =
					Max(max_distance, DistanceOfTwoPoints(m_node_coord[iNode].data(), m_node_coord[iNeibor].data()));
				neibor_set.insert(iNeibor);
			}
			// 以当前节点为中心，以最大距离为半径，找到范围内的节点
			double search_radius = max_distance * 1.00001;
			vtkNew<vtkIdList> result;
			while (result->GetNumberOfIds() < 6)
			{
				kdTree->FindPointsWithinRadius(search_radius, m_node_coord[iNode].data(), result);
				if (result->GetNumberOfIds() > 50)
				{
					search_radius *= 0.9;
					result->Reset();
				}
				else
					search_radius *= 1.1;
			}
			for (int i = 0; i < result->GetNumberOfIds(); ++i)
			{
				neibor_set.insert(result->GetId(i));
			}
			neibor_set.erase(iNode);
			nodeNeibor.resize(6); // 差分模板不改变
			for (auto& i : neibor_set)
			{
				// 如果nodeNeibor中没有该节点，添加该节点
				if (std::find(nodeNeibor.begin(), nodeNeibor.end(), i) == nodeNeibor.end())
					nodeNeibor.push_back(i);
			}
			neibor_num_after = nodeNeibor.size();
			if (neibor_num_after < neibor_num_before)
				Log::info("node:{} neibor num before:{} neibor num after:{}", iNode, neibor_num_before, neibor_num_after);
		}
	}

	void GridBuilderSYSU_FN::ReadCellFile()
	{
		std::ifstream fin;
		fin.open("cell.dat");
		int nodeNum;
		int cellNum;
		std::string separator = " \r\n\t#$;\"";
		std::string line;
		while (std::getline(fin, line))
		{
			if (line.empty())
				continue;
			std::transform(line.begin(), line.end(), line.begin(), [](unsigned char c) { return std::toupper(c); });
			line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
			size_t start_id_node_num = line.find("N=");
			if (start_id_node_num == std::string::npos)
				continue;
			start_id_node_num += 2;
			size_t end_id_node_num = line.find_first_of(',');
			std::string node_num_str = line.substr(start_id_node_num, end_id_node_num - start_id_node_num);
			line.erase(0, end_id_node_num + 1);
			size_t start_id_cell_num = line.find("E=");
			if (start_id_cell_num == std::string::npos)
				continue;
			start_id_cell_num += 2;
			size_t end_id_cell_num = line.find_first_of(',');
			std::string cell_num_str = line.substr(start_id_cell_num, end_id_cell_num - start_id_cell_num);
			nodeNum = std::stoi(node_num_str);
			cellNum = std::stoi(cell_num_str);
			break;
		}
		// 跳过节点坐标
		for (int i = 0; i < nodeNum; i++)
		{
			fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}

		Log::info("Total cell num:{}", cellNum);
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

	void GridBuilderSYSU_FN::CheckNode()
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
	void GridBuilderSYSU_FN::CheckUnkownNode()
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
	void GridBuilderSYSU_FN::CheckNeighborNum()
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
	void GridBuilderSYSU_FN::AddSelfToNeighbor()
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
	void GridBuilderSYSU_FN::SetBoundNeighbor()
	{
		return;
		// 构建节点KD?
		// 初始化vtk?
		vtkNew<vtkPoints> points;
		for (int i = 0; i < m_node_coord.size(); ++i)
		{
			points->InsertNextPoint(m_node_coord[i].data());
		}
		// 把vtk点转换为vtkploydata
		vtkNew<vtkPolyData> polydata;
		polydata->SetPoints(points);
		// 将vtkploydata转换为vtk点定位器
		vtkNew<vtkVertexGlyphFilter> glyphFilter;
		glyphFilter->SetInputData(polydata);
		glyphFilter->Update();
		// 构建KD�?
		vtkNew<vtkKdTreePointLocator> kdTree;
		kdTree->SetDataSet(glyphFilter->GetOutput());
		kdTree->BuildLocator();
		int point_num = kdTree->GetDataSet()->GetNumberOfPoints();
		int neibor_num_before = 0;
		int neibor_num_after = 0;
		int min_neibor_num = 1E5;
		int max_neibor_num = 0;
		int node_num = m_node_coord.size();
		// 扩展内部节点邻居节点，用于计算梯度
		std::set<int> neibor_set;
		// 消除边界点对应的邻居点也是边界点的情况
		for (int iBound = 0; iBound < m_bound_node.size(); iBound++)
		{
			auto& bound_node = m_bound_node[iBound];
			auto bound_index = bound_node.bound_index;
			auto ref_index = bound_node.ref_index;
			auto& bound_neibor = m_node_neibor[bound_index];
			auto& ref_neibor = m_node_neibor[ref_index];
			if (bound_neibor.size() == 0)
			{
				bound_neibor.push_back(ref_index);
			}
			if (m_node_type[ref_index] == NodeType::inner)
			{
				if (std::find(bound_neibor.begin(), bound_neibor.end(), ref_index) == bound_neibor.end())
				{
					bound_neibor.push_back(ref_index);
				}
			}
			if (bound_neibor.size() > 6)
			{
				continue;
			}
			else
			{
				neibor_set.clear();
				auto& nodeNeibor = m_node_neibor[bound_index];
				neibor_num_before = nodeNeibor.size();
				double max_distance = 0;
				for (auto& iNeibor : nodeNeibor)
				{
					max_distance =
						Max(max_distance, DistanceOfTwoPoints(m_node_coord[bound_index].data(), m_node_coord[iNeibor].data()));
					neibor_set.insert(iNeibor);
				}
				// 以当前节点为中心，以最大距离为半径，找到范围内的节点
				double search_radius = max_distance * 1.00001;
				vtkNew<vtkIdList> result;
				while (true)
				{
					kdTree->FindPointsWithinRadius(search_radius, m_node_coord[bound_index].data(), result);
					int inner_num = 0;
					for (int i = 0; i < result->GetNumberOfIds(); ++i)
					{
						if (m_node_type[result->GetId(i)] == NodeType::inner)
						{
							inner_num++;
						}
					}
					// Log::info("iBound={}, neighbor num: {}, inner num: {}, search radius: {:6E}", iBound, result->GetNumberOfIds(), inner_num, search_radius);
					if (inner_num > 6 || result->GetNumberOfIds() > 100)
						search_radius *= 0.95;
					else if ((inner_num != 0 && result->GetNumberOfIds() > 10) || (inner_num >= 3 && result->GetNumberOfIds() >= 7) || result->GetNumberOfIds() > 30)
						break;
					else
						search_radius *= 1.1;
				}
				for (int i = 0; i < result->GetNumberOfIds(); ++i)
				{
					neibor_set.insert(result->GetId(i));
					// if (m_node_type[result->GetId(i)] == NodeType::inner)
					// {
					//     auto &inner_neibor = m_node_neibor[result->GetId(i)];
					//     for (auto &iNeibor : inner_neibor)
					//     {
					//         neibor_set.insert(iNeibor);
					//     }
					// }
				}
				neibor_set.erase(bound_index);
				for (auto& i : neibor_set)
				{
					// 如果nodeNeibor中没有该节点，添加该节点
					if (std::find(nodeNeibor.begin(), nodeNeibor.end(), i) == nodeNeibor.end())
						nodeNeibor.push_back(i);
				}
				Log::info("iBound :{}/{}, neighbor num: {}", iBound, m_bound_node.size(), nodeNeibor.size());

				// Log::info("bound_index:{} neibor num: {}", bound_index, nodeNeibor.size());
				// std::cout << "bound index:" << bound_index << ": ";
				// double innner_num = 0;
				// for (int iNeigh = 0; iNeigh < nodeNeibor.size(); iNeigh++)
				// {
				//     std::cout << nodeNeibor[iNeigh] << " ";
				//     if (m_node_type[nodeNeibor[iNeigh]] == NodeType::inner)
				//     {
				//         innner_num++;
				//     }
				// }
				// std::cout << "inner num:" << innner_num << std::endl;
				neibor_num_after = nodeNeibor.size();
				min_neibor_num = Min(min_neibor_num, neibor_num_after);
				max_neibor_num = Max(max_neibor_num, neibor_num_after);
				if (neibor_num_after < neibor_num_before)
					Log::info("node:{} neibor num before:{} neibor num after:{}", bound_index, neibor_num_before, neibor_num_after);
			}
		}
		Log::info("Add inner node's neibor node to bound done");
	}
	void GridBuilderSYSU_FN::ConvertToGrid(shared_ptr<GridFN> grid)
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
		NodeFN* node = new NodeFN(node_num, neighbor_node_num.data(), neighbor_face_num.data(), neighbor_cell_num.data());
		for (size_t iNode = 0; iNode < m_node_coord.size(); iNode++)
		{
			node->SetCoord(iNode, m_node_coord[iNode].data());
			node->SetType(iNode, m_node_type[iNode]);
			node->SetNeighborNode(iNode, neighbor_node_num[iNode], m_node_neibor[iNode].data());
			node->SetNeighborFace(iNode, neighbor_face_num[iNode], nullptr);
			node->SetNeighborCell(iNode, neighbor_cell_num[iNode], nullptr);
		}
		grid->SetNode(node);
		CellFN* cell = new CellFN(m_cell_node.size());
		cell->SetNode(m_cell_node);
		dynamic_array<dynamic_array<index_type>> cell_node_face(m_cell_node.size());
		cell->SetFace(cell_node_face);
		double center[3];
		for (index_type iCell = 0; iCell < m_cell_node.size(); iCell++)
		{
			center[0] = center[1] = center[2] = 0;
			for (index_type iNode = 0; iNode < m_cell_node[iCell].size(); iNode++)
			{
				for (index_type i = 0; i < 3; i++)
				{
					center[i] += m_node_coord[m_cell_node[iCell][iNode]][i];
				}
			}
			for (index_type i = 0; i < 3; i++)
			{
				center[i] /= m_cell_node[iCell].size();
			}
			cell->SetCenter(iCell, center);
		}
		grid->SetCell(cell);

		FaceFN* face = new FaceFN();

		dynamic_array<index_type> face_node_num(m_bound_face.size());
		for (index_type iFace = 0; iFace < m_bound_face.size(); iFace++)
		{
			face_node_num[iFace] = m_bound_face[iFace].face_node.size();
		}
		face->Allocate(m_bound_face.size(), face_node_num.data());
		for (index_type iFace = 0; iFace < m_bound_face.size(); iFace++)
		{
			face->SetFace2Node(iFace, m_bound_face[iFace].face_node.data(), m_bound_face[iFace].face_node.size());
			face->SetNormal(iFace, m_bound_face[iFace].normal.data());
			face->SetArea(iFace, m_bound_face[iFace].area);
		}
		grid->SetFace(face);

		BoundManagerFN* boundary_map = new BoundManagerFN();
		for (int iFace = 0; iFace < m_bound_node.size(); iFace++)
		{
			boundary_map->AddBoundary(m_bound_node[iFace].type, BoundFN(m_bound_node[iFace].bound_index, m_bound_node[iFace].ref_index, 0, m_bound_node[iFace].normal));
		}
		grid->SetBoundaryMap(boundary_map);
	}
	void GridBuilderSYSU_FN::ReadBoundFile()
	{
		std::ifstream fin;
		fin.open("bound.dat");
		int nodeNum;
		int bound_face_num;
		std::string separator = " \r\n\t#$;\"";
		std::string line;
		while (std::getline(fin, line))
		{
			if (line.empty())
				continue;
			std::transform(line.begin(), line.end(), line.begin(), [](unsigned char c) { return std::toupper(c); });
			line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
			size_t start_id_node_num = line.find("N=");
			if (start_id_node_num == std::string::npos)
				continue;
			start_id_node_num += 2;
			size_t end_id_node_num = line.find_first_of(',');
			std::string node_num_str = line.substr(start_id_node_num, end_id_node_num - start_id_node_num);
			line.erase(0, end_id_node_num + 1);
			size_t start_id_cell_num = line.find("E=");
			if (start_id_cell_num == std::string::npos)
				continue;
			start_id_cell_num += 2;
			size_t end_id_cell_num = line.find_first_of(',');
			std::string cell_num_str = line.substr(start_id_cell_num, end_id_cell_num - start_id_cell_num);
			nodeNum = std::stoi(node_num_str);
			bound_face_num = std::stoi(cell_num_str);
			break;
		}
		// 跳过节点坐标
		for (int i = 0; i < nodeNum; i++)
		{
			fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}
		m_bound_face.resize(bound_face_num);
		for (int iBound = 0; iBound < bound_face_num; iBound++)
		{
			auto& bound_face = m_bound_face[iBound];
			int boundNodeNum = 4;
			bound_face.face_node.resize(boundNodeNum);
			for (int i = 0; i < boundNodeNum; i++)
			{
				fin >> bound_face.face_node[i];
				bound_face.face_node[i] -= 1;
			}
			if (bound_face.face_node[2] == bound_face.face_node[3])
			{
				bound_face.face_node.pop_back();
			}
		}
		fin.close();
		double v1[3], v2[3];
		for (int iBound = 0; iBound < bound_face_num; iBound++)
		{
			auto& bound_face = m_bound_face[iBound];
			bound_face.normal.resize(3);
			if (bound_face.face_node.size() == 3)
			{
				bound_face.area =
					TriangleArea(m_node_coord[bound_face.face_node[0]].data(), m_node_coord[bound_face.face_node[1]].data(),
						m_node_coord[bound_face.face_node[2]].data());
				for (int i = 0; i < 3; i++)
				{
					v1[i] = m_node_coord[bound_face.face_node[1]][i] - m_node_coord[bound_face.face_node[0]][i];
					v2[i] = m_node_coord[bound_face.face_node[2]][i] - m_node_coord[bound_face.face_node[0]][i];
				}
			}
			else
			{
				bound_face.area = QuadrangleArea(
					m_node_coord[bound_face.face_node[0]].data(), m_node_coord[bound_face.face_node[1]].data(),
					m_node_coord[bound_face.face_node[2]].data(), m_node_coord[bound_face.face_node[3]].data());
				for (int i = 0; i < 3; i++)
				{
					v1[i] = m_node_coord[bound_face.face_node[2]][i] - m_node_coord[bound_face.face_node[0]][i];
					v2[i] = m_node_coord[bound_face.face_node[3]][i] - m_node_coord[bound_face.face_node[1]][i];
				}
			}
			CrossProduct(v1, v2, bound_face.normal.data());
			double norm = sqrt(bound_face.normal[0] * bound_face.normal[0] + bound_face.normal[1] * bound_face.normal[1] +
				bound_face.normal[2] * bound_face.normal[2]);
			for (int i = 0; i < 3; i++)
			{
				bound_face.normal[i] /= norm;
			}
		}
	}

} // namespace zaran