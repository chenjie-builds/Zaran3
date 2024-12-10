#include "NodeFNFDM.h"
#include "Log.h"
namespace zaran
{
	NodeFN::NodeFN(index_type node_num)
	{
		SetCount(node_num);
		m_coordinate = new double[3 * node_num];
		m_type = new NodeType[node_num];
		m_node_num = new index_type[node_num];
		m_node_offset = new index_type[node_num];
		m_face_num = new index_type[node_num];
		m_face_idx = new index_type[node_num];
		m_cell_num = new index_type[node_num];
		m_neighor_cell_index = new index_type[node_num];
		m_node_idx = nullptr;
		m_neighor_face = nullptr;
		m_neighor_cell = nullptr;
	}

	NodeFN::NodeFN(index_type node_num, const index_type*neighbor_node_num, const index_type*neighbor_face_num, const index_type *neighbor_cell_num)
	{
		SetCount(node_num);
		m_coordinate = new double[3 * node_num];
		m_type = new NodeType[node_num];
		m_node_num = new index_type[node_num];
		for (index_type i = 0; i < node_num; i++)
		{
			m_node_num[i] = neighbor_node_num[i];
		}
		m_node_offset = new index_type[node_num];
		m_face_num = new index_type[node_num];
		m_face_idx = new index_type[node_num];
		m_cell_num = new index_type[node_num];
		m_neighor_cell_index = new index_type[node_num];
		index_type sum_node, sum_face, sum_cell;
		sum_node = sum_face = sum_cell = 0;
		for (index_type i = 0; i < node_num; i++)
		{
			m_node_offset[i] = sum_node;
			m_face_idx[i] = sum_face;
			m_neighor_cell_index[i] = sum_cell;
			sum_node += neighbor_node_num[i];
			sum_face += neighbor_face_num[i];
			sum_cell += neighbor_cell_num[i];
		}
		m_node_idx = new index_type[sum_node];
		m_neighor_face = new index_type[sum_face];
		m_neighor_cell = new index_type[sum_cell];
	}
	NodeFN::~NodeFN()
	{
		delete[] m_coordinate;
		delete[] m_type;
		delete[] m_node_num;
		delete[] m_node_offset;
		delete[] m_face_num;
		delete[] m_face_idx;
		delete[] m_cell_num;
		delete[] m_neighor_cell_index;
		delete[] m_node_idx;
		delete[] m_neighor_face;
		delete[] m_neighor_cell;
	}
	void NodeFN::SetCoord(index_type idx, const double *coord)
	{
		for (index_type i = 0; i < 3; i++)
		{
			m_coordinate[idx * 3 + i] = coord[i];
		}
	}
	void NodeFN::SetType(index_type idx, NodeType type)
	{
		m_type[idx] = type;
	}
	void NodeFN::SetNeighborNode(index_type node_num, const index_type *neighbor_node_num, const index_type *neighbor_node)
	{
		if (node_num != GetCount())
		{
			Log::error("NodeFN::SetNeighborNode: node_num != GetCount()");
			exit(0);
		}
		if (m_node_num != nullptr)
		{
			delete[] m_node_num;
		}
		m_node_num = new index_type[node_num];
		for (index_type i = 0; i < node_num; i++)
		{
			m_node_num[i] = neighbor_node_num[i];
		}
		index_type sum_node = 0;
		for (index_type i = 0; i < node_num; i++)
		{
			m_node_offset[i] = sum_node;
			sum_node += neighbor_node_num[i];
		}
		if (m_node_idx != nullptr)
		{
			delete[] m_node_idx;
		}
		m_node_idx = new index_type[sum_node];
		for (index_type i = 0; i < sum_node; i++)
		{
			m_node_idx[i] = neighbor_node[i];
		}
	}
	void NodeFN::SetNeighborNode(index_type idx, index_type neighbor_num, const index_type *neighbor_cloud)
	{
		if (neighbor_num != m_node_num[idx])
		{
			index_type delta = neighbor_num - m_node_num[idx];
			m_node_num[idx] = neighbor_num;
			for (index_type i = idx + 1; i < GetCount(); i++)
			{
				m_node_offset[i] += delta;
			}
			index_type sum = 0;
			for (index_type i = 0; i < GetCount(); i++)
			{
				sum += m_node_num[i];
			}
			// 重新分配内存
			index_type *temp = new index_type[sum];
			// 拷贝数据 0-->idx
			for (index_type i = 0; i < idx; i++)
			{
				for (index_type j = 0; j < m_node_num[i]; j++)
				{
					temp[m_node_offset[i] + j] = m_node_idx[m_node_offset[i] + j];
				}
			}
			// 拷贝数据 idx
			for (index_type i = 0; i < neighbor_num; i++)
			{
				temp[m_node_offset[idx] + i] = neighbor_cloud[i];
			}
			// 拷贝数据 idx+1-->node_num
			for (index_type i = idx + 1; i < GetCount(); i++)
			{
				for (index_type j = 0; j < m_node_num[i]; j++)
				{
					temp[m_node_offset[i] + j] = m_node_idx[m_node_offset[i] - delta + j];
				}
			}
			delete[] m_node_idx;
			m_node_idx = temp;
		}
		else
		{
			for (index_type i = 0; i < neighbor_num; i++)
			{
				m_node_idx[m_node_offset[idx] + i] = neighbor_cloud[i];
			}
		}
	}
	void NodeFN::SetNeighborFace(index_type idx, index_type neighbor_num, const index_type *neighbor_face)
	{
		for (index_type i = 0; i < neighbor_num; i++)
		{
			m_neighor_face[m_face_idx[idx] + i] = neighbor_face[i];
		}
	}
	void NodeFN::SetNeighborCell(index_type idx, index_type neighbor_num, const index_type *neighbor_cell)
	{
		for (index_type i = 0; i < neighbor_num; i++)
		{
			m_neighor_cell[m_neighor_cell_index[idx] + i] = neighbor_cell[i];
		}
	}
	const NodeType &NodeFN::GetType(const index_type &index) const
	{
		return m_type[index];
	}
	const double *NodeFN::GetCoord(const index_type &index) const
	{
		return m_coordinate + index * 3;
	}
	const index_type *NodeFN::GetNeighborNode(const index_type &idx) const
	{
		return m_node_idx + m_node_offset[idx];
	}
	index_type NodeFN::GetNeighborNodeNum(const index_type &idx) const
	{
		return m_node_num[idx];
	}
	const index_type *NodeFN::GetNeighborFace(const index_type &idx) const
	{
		return m_neighor_face + m_face_idx[idx];
	}
	index_type NodeFN::GetNeighborFaceNum(const index_type &idx) const
	{
		return m_face_num[idx];
	}
	const index_type *NodeFN::GetNeighborCell(const index_type &idx) const
	{
		return m_neighor_cell + m_neighor_cell_index[idx];
	}
	index_type NodeFN::GetNeighborCellNum(const index_type &idx) const
	{
		return m_cell_num[idx];
	}

}