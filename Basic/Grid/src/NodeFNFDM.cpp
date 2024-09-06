#include "NodeFNFDM.h"
namespace zaran
{
	NodeFN::NodeFN(int node_num)		
	{
		SetCount(node_num);
		m_coordinate = new double[3 * node_num];
		m_type = new NodeType[node_num];
		m_neighor_node_num = new int[node_num];
		m_neighor_node_index = new int[node_num];
		m_neighor_face_num = new int[node_num];
		m_neighor_face_index = new int[node_num];
		m_neighor_cell_num = new int[node_num];
		m_neighor_cell_index = new int[node_num];
		m_neighor_node = nullptr;
		m_neighor_face = nullptr;
		m_neighor_cell = nullptr;
	}

	NodeFN::NodeFN(int node_num,const int* neighbor_node_num,const int* neighbor_face_num, const int* neighbor_cell_num)
	{
		SetCount(node_num);
		m_coordinate = new double[3 * node_num];
		m_type = new NodeType[node_num];
		m_neighor_node_num = new int[node_num];
		for (int i = 0; i < node_num; i++)
		{
			m_neighor_node_num[i] = neighbor_node_num[i];
		}
		m_neighor_node_index = new int[node_num];
		m_neighor_face_num = new int[node_num];
		m_neighor_face_index = new int[node_num];
		m_neighor_cell_num = new int[node_num];
		m_neighor_cell_index = new int[node_num];
		int sum_node, sum_face, sum_cell;
		sum_node = sum_face = sum_cell = 0;
		for (int i = 0; i < node_num; i++)
		{
			m_neighor_node_index[i] = sum_node;
			m_neighor_face_index[i] = sum_face;
			m_neighor_cell_index[i] = sum_cell;
			sum_node += neighbor_node_num[i];
			sum_face += neighbor_face_num[i];
			sum_cell += neighbor_cell_num[i];
		}
		m_neighor_node = new int[sum_node];
		m_neighor_face = new int[sum_face];
		m_neighor_cell = new int[sum_cell];
	}
	NodeFN::~NodeFN()
	{
		delete[] m_coordinate;
		delete[] m_type;
		delete[] m_neighor_node_num;
		delete[] m_neighor_node_index;
		delete[] m_neighor_face_num;
		delete[] m_neighor_face_index;
		delete[] m_neighor_cell_num;
		delete[] m_neighor_cell_index;
		delete[] m_neighor_node;
		delete[] m_neighor_face;
		delete[] m_neighor_cell;
	}
	void NodeFN::SetCoord(int idx, const double* coord)
	{
		for (int i = 0; i < 3; i++)
		{
			m_coordinate[idx * 3 + i] = coord[i];
		}
	}
	void NodeFN::SetType(int idx, NodeType& type)
	{
		m_type[idx] = type;
	}
	void NodeFN::SetNeighborNode(int idx, int neighbor_num, const int* neighbor_cloud)
	{
		if (neighbor_num != m_neighor_node_num[idx])
		{
			int delta = neighbor_num - m_neighor_node_num[idx];
			m_neighor_node_num[idx] = neighbor_num;
			for (int i = idx + 1; i < GetCount(); i++)
			{
				m_neighor_node_index[i] += delta;
			}
			int sum = 0;
			for (int i = 0; i < GetCount(); i++)
			{
				sum += m_neighor_node_num[i];
			}
			//重新分配内存
			int* temp = new int[sum];
			//拷贝数据 0-->idx
			for (int i = 0; i < idx; i++)
			{
				for (int j = 0; j < m_neighor_node_num[i]; j++)
				{
					temp[m_neighor_node_index[i] + j] = m_neighor_node[m_neighor_node_index[i] + j];
				}
			}
			//拷贝数据 idx
			for (int i = 0; i < neighbor_num; i++)
			{
				temp[m_neighor_node_index[idx] + i] = neighbor_cloud[i];
			}
			//拷贝数据 idx+1-->node_num
			for (int i = idx + 1; i < GetCount(); i++)
			{
				for (int j = 0; j < m_neighor_node_num[i]; j++)
				{
					temp[m_neighor_node_index[i] + j] = m_neighor_node[m_neighor_node_index[i] - delta + j];
				}
			}
			delete[] m_neighor_node;
			m_neighor_node = temp;
		}
		else
		{
			for (int i = 0; i < neighbor_num; i++)
			{
				m_neighor_node[m_neighor_node_index[idx] + i] = neighbor_cloud[i];
			}
		}
	}
	void NodeFN::SetNeighborFace(int idx, int neighbor_num, const int* neighbor_face)
	{
		for (int i = 0; i < neighbor_num; i++)
		{
			m_neighor_face[m_neighor_face_index[idx] + i] = neighbor_face[i];
		}
	}
	void NodeFN::SetNeighborCell(int idx, int neighbor_num, const int* neighbor_cell)
	{
		for (int i = 0; i < neighbor_num; i++)
		{
			m_neighor_cell[m_neighor_cell_index[idx] + i] = neighbor_cell[i];
		}
	}
	const NodeType& NodeFN::GetType(const int& index)const
	{
		return m_type[index];
	}
	const double* NodeFN::GetCoord(const int& index)const
	{
		return m_coordinate + index * 3;
	}
	const int* NodeFN::GetNeighborNode(const int& idx)const
	{
		return m_neighor_node + m_neighor_node_index[idx];
	}
	int NodeFN::GetNeighborNodeNum(const int& idx)const
	{
		return m_neighor_node_num[idx];
	}
	const int* NodeFN::GetNeighborFace(const int& idx)const
	{
		return m_neighor_face + m_neighor_face_index[idx];
	}
	int NodeFN::GetNeighborFaceNum(const int& idx)const
	{
		return m_neighor_face_num[idx];
	}
	const int* NodeFN::GetNeighborCell(const int& idx)const
	{
		return m_neighor_cell + m_neighor_cell_index[idx];
	}
	int NodeFN::GetNeighborCellNum(const int& idx)const
	{
		return m_neighor_cell_num[idx];
	}

}