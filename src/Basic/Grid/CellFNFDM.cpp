#include "CellFNFDM.h"
#include "Log.h"
namespace zaran
{

CellFN::CellFN( int cell_num):CellBase()
{
    SetCellNum(cell_num);
    m_node_num = new index_type[cell_num];
    m_node_id = new index_type[cell_num];
    m_cell_face_num = new index_type[cell_num];
    m_face_id = new index_type[cell_num];
    m_center = new double[3 * cell_num];
    m_node = nullptr;
    m_face = nullptr;
}

CellFN::~CellFN()
{
    delete[] m_node_num;
    delete[] m_node_id;
    delete[] m_node;
    delete[] m_cell_face_num;
    delete[] m_face_id;
    delete[] m_face;
    delete[] m_center;
}
void CellFN::SetNode(dynamic_array<dynamic_array<index_type>>&node_index)
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
    if (m_node != nullptr)
    {
        delete[] m_node;
        m_node = nullptr;
    }
    else
    {
        m_node = new index_type[node_sum];
    }
    for (int i = 0; i < GetCellNum(); i++)
    {
        m_node_num[i] = node_index[i].size();
        for (int j = 0; j < m_node_num[i]; j++)
        {
            m_node[m_node_id[i] + j] = node_index[i][j];
        }
    }
}
void CellFN::SetFace(dynamic_array<dynamic_array<index_type>>&face_index)
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
    if (m_face != nullptr)
    {
        delete[] m_face;
        m_face = nullptr;
    }
    else
    {
        m_face = new index_type[face_sum];
    }
    for (int i = 0; i < GetCellNum(); i++)
    {
        m_cell_face_num[i] = face_index[i].size();
        for (int j = 0; j < m_cell_face_num[i]; j++)
        {
            m_face[m_face_id[i] + j] = face_index[i][j];
        }
    }
}
void CellFN::SetCenter(index_type iCell, const double *center)
{
	if(iCell >= GetCellNum())
	{
		Log::error("CellTopo::SetCenter: The cell index is out of range");
		Log::error(" iCell = {}, m_cell_num = {}", iCell, GetCellNum());
		exit(0);
	}
    for (size_t i = 0; i < 3; i++)
    {
        m_center[iCell * 3 + i] = center[i];
    }
}

index_type CellFN::GetNodeNum(index_type iCell) const
{
    return m_node_num[iCell];
}
const index_type*CellFN::GetNode(index_type iCell)
{
    return m_node + m_node_id[iCell];
}
index_type CellFN::GetFaceNum(index_type iCell) const
{
	return m_cell_face_num[iCell];
}
const index_type*CellFN::GetFace(index_type iCell)
{
	return m_face + m_face_id[iCell];
}
const double *CellFN::GetCenterCoord(index_type iCell) const
{
    return m_center + iCell * 3;
}
} // namespace zaran