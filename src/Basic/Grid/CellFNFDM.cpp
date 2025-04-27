#include "CellFNFDM.h"
#include "Log.h"
namespace zaran
{

	CellFN::CellFN(int cell_num) :CellBase()
	{
		SetCellNum(cell_num);
	}

	CellFN::~CellFN()
	{

	}

	void CellFN::SetCellNum(int cell_num)
	{
		if (cell_num <= 0)
		{
			Log::error("CellTopo::SetCellNum: The number of cells is less than 0");
			Log::error(" cell_num = {}", cell_num);
			exit(0);
		}
		CellBase::SetCellNum(cell_num);
		m_node_num.resize(cell_num);
		m_node_id.resize(cell_num);
		m_cell_face_num.resize(cell_num);
		m_face_id.resize(cell_num);
		m_center.resize(cell_num * 3);
	}

	void CellFN::SetNode(dynamic_array<dynamic_array<index_type>>& node_index)
	{
		if (GetCellNum() != node_index.size())
		{
			Log::error("CellTopo::SetNode: The number of cells is not equal to the number of nodes");
			Log::error(" cell_num = {}, node_index.size() = {}", GetCellNum(), node_index.size());
			exit(0);
		}
		int node_sum = 0;
		for (int i = 0; i < GetCellNum(); i++)
		{
			m_node_id[i] = node_sum;
			node_sum += node_index[i].size();
		}
		m_node.resize(node_sum);
		for (int i = 0; i < GetCellNum(); i++)
		{
			m_node_num[i] = node_index[i].size();
			for (int j = 0; j < m_node_num[i]; j++)
			{
				m_node[m_node_id[i] + j] = node_index[i][j];
			}
		}
	}
	void CellFN::SetFace(dynamic_array<dynamic_array<index_type>>& face_index)
	{
		if (GetCellNum() != face_index.size())
		{
			Log::error("CellTopo::SetFace: The number of cells is not equal to the number of faces");
			Log::error(" cell_num = {}, face_index.size() = {}", GetCellNum(), face_index.size());
			exit(0);
		}
		int face_sum = 0;
		for (size_t i = 0; i < GetCellNum(); i++)
		{
			m_face_id[i] = face_sum;
			face_sum += face_index[i].size();
		}
		m_face.resize(face_sum);
		for (int i = 0; i < GetCellNum(); i++)
		{
			m_cell_face_num[i] = face_index[i].size();
			for (int j = 0; j < m_cell_face_num[i]; j++)
			{
				m_face[m_face_id[i] + j] = face_index[i][j];
			}
		}
	}
	void CellFN::SetCenterCoord(index_type iCell, const double* center_coord)
	{
		if (iCell >= GetCellNum())
		{
			Log::error("CellTopo::SetCenter: The cell index is out of range");
			Log::error(" iCell = {}, m_cell_num = {}", iCell, GetCellNum());
			exit(0);
		}
		for (size_t i = 0; i < 3; i++)
		{
			m_center[iCell * 3 + i] = center_coord[i];
		}
	}

	index_type CellFN::GetNodeNum(index_type iCell) const
	{
		return m_node_num[iCell];
	}
	 span<index_type>  CellFN::GetNode(index_type iCell) 
	{
		return  span<index_type>(m_node.data() + m_node_id[iCell], m_node_num[iCell]);
	}
	index_type CellFN::GetFaceNum(index_type iCell) const
	{
		return m_cell_face_num[iCell];
	}
	span<index_type> CellFN::GetFace(index_type iCell)
	{
		return span<index_type>(m_face.data() + m_face_id[iCell], m_cell_face_num[iCell]);
	}
	span<double> CellFN::GetCenterCoord(index_type iCell) 
	{
		return span<double>(m_center.data() + iCell * 3, 3);
	}
} // namespace zaran