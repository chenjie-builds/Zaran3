#include "NodeFNFDM.h"
#include "Log.h"
namespace zaran
{
	NodeFN::NodeFN(index_type node_num)
	{
		SetNodeNum(node_num);
	}

	NodeFN::NodeFN(index_type node_num, const index_type* neighbor_node_num, const index_type* neighbor_face_num, const index_type* neighbor_cell_num)
	{
		SetNodeNum(node_num);
		m_node_num.resize(node_num);
		for (index_type i = 0; i < node_num; i++)
		{
			m_node_num[i] = neighbor_node_num[i];
		}
		m_node_offset.resize(node_num);
		m_face_num.resize(node_num);
		m_face_idx.resize(node_num);
		m_cell_num.resize(node_num);
		m_neighor_cell_index.resize(node_num);
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
		m_node_idx.resize(node_num);
		m_neighor_face.resize(node_num);
		m_neighor_cell.resize(node_num);
	}

	NodeFN::NodeFN()
	{
		SetNodeNum(0);
	}

	void NodeFN::SetNodeNum(index_type node_num)
	{
		SetCount(node_num);
		m_coordinate.resize(node_num * 3);
		m_type.resize(node_num);
		m_node_num.resize(node_num);
		m_node_offset.resize(node_num);
		m_face_num.resize(node_num);
		m_face_idx.resize(node_num);
		m_cell_num.resize(node_num);
		m_neighor_cell_index.resize(node_num);
	}

	NodeFN::~NodeFN()
	{
	}
	void NodeFN::SetCoord(index_type idx, const double* coord)
	{
		for (index_type i = 0; i < 3; i++)
		{
			m_coordinate[idx * 3 + i] = coord[i];
		}
	}

	void NodeFN::SetCoord(const dynamic_array<dynamic_array<double>>& coord)
	{
#ifdef DEBUG_MODE
		if (coord.size() != m_coordinate.size())
		{
			Log::error("NodeFN::SetCoord: coord.size() != m_coordinate.size()");
			Log::error("coord.size() = {}, m_coordinate.size() = {}", coord.size(), m_coordinate.size());
			exit(0);
		}
#endif
		for (index_type i = 0; i < coord.size(); i++)
		{
			for (index_type j = 0; j < 3; j++)
			{
				m_coordinate[i * 3 + j] = coord[i][j];
			}
		}
	}

	void NodeFN::SetType(index_type idx, NodeType type)
	{
		m_type[idx] = type;
	}

	void NodeFN::SetType(const dynamic_array<NodeType>& type)
	{
#ifdef DEBUG_MODE
		if (type.size() != m_type.size())
		{
			Log::error("NodeFN::SetType: type.size() != m_type.size()");
			Log::error("type.size() = {}, m_type.size() = {}", type.size(), m_type.size());
			exit(0);
		}
#endif
		m_type = type;
	}

	void NodeFN::SetNeighborNode(index_type idx, index_type neighbor_num, const index_type* neighbor_cloud)
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
			dynamic_array<index_type> temp(sum);
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

	void NodeFN::SetNeighborNode(dynamic_array<dynamic_array<index_type>>& neighbor_node)
	{
		index_type node_num = neighbor_node.size();
#ifdef DEBUG_MODE
		if (node_num != GetCount())
		{
			Log::error("NodeFN::SetNeighborNode: node_num != GetCount()");
			Log::error("node_num = {}, GetCount() = {}", node_num, GetCount());
			exit(0);
		}
#endif
		index_type sum_node = 0;
		for (index_type i = 0; i < node_num; i++)
		{
			m_node_num[i] = neighbor_node[i].size();
			m_node_offset[i] = sum_node;
			sum_node += m_node_num[i];
		}
		m_node_idx.resize(sum_node);
		for (index_type i = 0; i < node_num; i++)
		{
			for (index_type j = 0; j < m_node_num[i]; j++)
			{
				m_node_idx[m_node_offset[i] + j] = neighbor_node[i][j];
			}
		}
	}

	void NodeFN::SetNeighborFace(index_type idx, index_type neighbor_num, const index_type* neighbor_face)
	{
		for (index_type i = 0; i < neighbor_num; i++)
		{
			m_neighor_face[m_face_idx[idx] + i] = neighbor_face[i];
		}
	}

	void NodeFN::SetNeighborFace(dynamic_array<dynamic_array<index_type>>& neighbor_face)
	{
		index_type node_num = neighbor_face.size();
#ifdef DEBUG_MODE
		if (node_num != GetCount())
		{
			Log::error("NodeFN::SetNeighborFace: node_num != GetCount()");
			Log::error("node_num = {}, GetCount() = {}", node_num, GetCount());
			exit(0);
		}
#endif
		index_type sum_face = 0;
		for (index_type i = 0; i < node_num; i++)
		{
			m_face_num[i] = neighbor_face[i].size();
			m_face_idx[i] = sum_face;
			sum_face += m_face_num[i];
		}
		m_neighor_face.resize(sum_face);
		for (index_type i = 0; i < node_num; i++)
		{
			for (index_type j = 0; j < m_face_num[i]; j++)
			{
				m_neighor_face[m_face_idx[i] + j] = neighbor_face[i][j];
			}
		}

	}

	void NodeFN::SetNeighborCell(index_type idx, index_type neighbor_num, const index_type* neighbor_cell)
	{
		for (index_type i = 0; i < neighbor_num; i++)
		{
			m_neighor_cell[m_neighor_cell_index[idx] + i] = neighbor_cell[i];
		}
	}

	void NodeFN::SetNeighborCell(dynamic_array<dynamic_array<index_type>>& neighbor_cell)
	{
		index_type node_num = neighbor_cell.size();
#ifdef DEBUG_MODE
		if (node_num != GetCount())
		{
			Log::error("NodeFN::SetNeighborCell: node_num != GetCount()");
			Log::error("node_num = {}, GetCount() = {}", node_num, GetCount());
			exit(0);
		}
#endif
		index_type sum_cell = 0;
		for (index_type i = 0; i < node_num; i++)
		{
			m_cell_num[i] = neighbor_cell[i].size();
			m_neighor_cell_index[i] = sum_cell;
			sum_cell += m_cell_num[i];
		}
		m_neighor_cell.resize(sum_cell);
		for (index_type i = 0; i < node_num; i++)
		{
			for (index_type j = 0; j < m_cell_num[i]; j++)
			{
				m_neighor_cell[m_neighor_cell_index[i] + j] = neighbor_cell[i][j];
			}
		}
	}

	const NodeType& NodeFN::GetType(const index_type& index) const
	{
		return m_type[index];
	}
	span<double> NodeFN::GetCoord(const index_type& index)
	{
		return span<double>(m_coordinate.data() + index * 3, 3);
	}
	span<index_type> NodeFN::GetNeighborNode(index_type idx)
	{
		return span<index_type>(m_node_idx.data() + m_node_offset[idx], m_node_num[idx]);
	}

	index_type NodeFN::GetNeighborNodeNum(const index_type& idx) const
	{
		return m_node_num[idx];
	}

	span<index_type> NodeFN::GetNeighborFace(index_type idx)
	{
		return span<index_type>(m_neighor_face.data() + m_face_idx[idx], m_face_num[idx]);
	}

	index_type NodeFN::GetNeighborFaceNum(const index_type& idx) const
	{
		return m_face_num[idx];
	}
	span<index_type> NodeFN::GetNeighborCell(index_type idx)
	{
		return span<index_type>(m_neighor_cell.data() + m_neighor_cell_index[idx], m_cell_num[idx]);
	}

	index_type NodeFN::GetNeighborCellNum(const index_type& idx) const
	{
		return m_cell_num[idx];
	}

}