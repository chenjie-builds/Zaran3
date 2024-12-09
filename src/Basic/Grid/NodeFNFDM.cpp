#include "NodeFNFDM.h"
#include "Log.h"
namespace zaran
{
	NodeFN::NodeFN(Id node_num)
	{
		SetCount(node_num);
		m_coordinate = new double[3 * node_num];
		m_type = new NodeType[node_num];
		m_node_num = new Id[node_num];
		m_node_offset = new Id[node_num];
		m_face_num = new Id[node_num];
		m_face_idx = new Id[node_num];
		m_cell_num = new Id[node_num];
		m_neighor_cell_index = new Id[node_num];
		m_node_idx = nullptr;
		m_neighor_face = nullptr;
		m_neighor_cell = nullptr;
	}

	NodeFN::NodeFN(Id node_num, const Id*neighbor_node_num, const Id*neighbor_face_num, const Id *neighbor_cell_num)
	{
		SetCount(node_num);
		m_coordinate = new double[3 * node_num];
		m_type = new NodeType[node_num];
		m_node_num = new Id[node_num];
		for (Id i = 0; i < node_num; i++)
		{
			m_node_num[i] = neighbor_node_num[i];
		}
		m_node_offset = new Id[node_num];
		m_face_num = new Id[node_num];
		m_face_idx = new Id[node_num];
		m_cell_num = new Id[node_num];
		m_neighor_cell_index = new Id[node_num];
		Id sum_node, sum_face, sum_cell;
		sum_node = sum_face = sum_cell = 0;
		for (Id i = 0; i < node_num; i++)
		{
			m_node_offset[i] = sum_node;
			m_face_idx[i] = sum_face;
			m_neighor_cell_index[i] = sum_cell;
			sum_node += neighbor_node_num[i];
			sum_face += neighbor_face_num[i];
			sum_cell += neighbor_cell_num[i];
		}
		m_node_idx = new Id[sum_node];
		m_neighor_face = new Id[sum_face];
		m_neighor_cell = new Id[sum_cell];
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
	void NodeFN::SetCoord(Id idx, const double *coord)
	{
		for (Id i = 0; i < 3; i++)
		{
			m_coordinate[idx * 3 + i] = coord[i];
		}
	}
	void NodeFN::SetType(Id idx, NodeType type)
	{
		m_type[idx] = type;
	}
	void NodeFN::SetNeighborNode(Id node_num, const Id *neighbor_node_num, const Id *neighbor_node)
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
		m_node_num = new Id[node_num];
		for (Id i = 0; i < node_num; i++)
		{
			m_node_num[i] = neighbor_node_num[i];
		}
		Id sum_node = 0;
		for (Id i = 0; i < node_num; i++)
		{
			m_node_offset[i] = sum_node;
			sum_node += neighbor_node_num[i];
		}
		if (m_node_idx != nullptr)
		{
			delete[] m_node_idx;
		}
		m_node_idx = new Id[sum_node];
		for (Id i = 0; i < sum_node; i++)
		{
			m_node_idx[i] = neighbor_node[i];
		}
	}
	void NodeFN::SetNeighborNode(Id idx, Id neighbor_num, const Id *neighbor_cloud)
	{
		if (neighbor_num != m_node_num[idx])
		{
			Id delta = neighbor_num - m_node_num[idx];
			m_node_num[idx] = neighbor_num;
			for (Id i = idx + 1; i < GetCount(); i++)
			{
				m_node_offset[i] += delta;
			}
			Id sum = 0;
			for (Id i = 0; i < GetCount(); i++)
			{
				sum += m_node_num[i];
			}
			// 重新分配内存
			Id *temp = new Id[sum];
			// 拷贝数据 0-->idx
			for (Id i = 0; i < idx; i++)
			{
				for (Id j = 0; j < m_node_num[i]; j++)
				{
					temp[m_node_offset[i] + j] = m_node_idx[m_node_offset[i] + j];
				}
			}
			// 拷贝数据 idx
			for (Id i = 0; i < neighbor_num; i++)
			{
				temp[m_node_offset[idx] + i] = neighbor_cloud[i];
			}
			// 拷贝数据 idx+1-->node_num
			for (Id i = idx + 1; i < GetCount(); i++)
			{
				for (Id j = 0; j < m_node_num[i]; j++)
				{
					temp[m_node_offset[i] + j] = m_node_idx[m_node_offset[i] - delta + j];
				}
			}
			delete[] m_node_idx;
			m_node_idx = temp;
		}
		else
		{
			for (Id i = 0; i < neighbor_num; i++)
			{
				m_node_idx[m_node_offset[idx] + i] = neighbor_cloud[i];
			}
		}
	}
	void NodeFN::SetNeighborFace(Id idx, Id neighbor_num, const Id *neighbor_face)
	{
		for (Id i = 0; i < neighbor_num; i++)
		{
			m_neighor_face[m_face_idx[idx] + i] = neighbor_face[i];
		}
	}
	void NodeFN::SetNeighborCell(Id idx, Id neighbor_num, const Id *neighbor_cell)
	{
		for (Id i = 0; i < neighbor_num; i++)
		{
			m_neighor_cell[m_neighor_cell_index[idx] + i] = neighbor_cell[i];
		}
	}
	const NodeType &NodeFN::GetType(const Id &index) const
	{
		return m_type[index];
	}
	const double *NodeFN::GetCoord(const Id &index) const
	{
		return m_coordinate + index * 3;
	}
	const Id *NodeFN::GetNeighborNode(const Id &idx) const
	{
		return m_node_idx + m_node_offset[idx];
	}
	Id NodeFN::GetNeighborNodeNum(const Id &idx) const
	{
		return m_node_num[idx];
	}
	const Id *NodeFN::GetNeighborFace(const Id &idx) const
	{
		return m_neighor_face + m_face_idx[idx];
	}
	Id NodeFN::GetNeighborFaceNum(const Id &idx) const
	{
		return m_face_num[idx];
	}
	const Id *NodeFN::GetNeighborCell(const Id &idx) const
	{
		return m_neighor_cell + m_neighor_cell_index[idx];
	}
	Id NodeFN::GetNeighborCellNum(const Id &idx) const
	{
		return m_cell_num[idx];
	}

}