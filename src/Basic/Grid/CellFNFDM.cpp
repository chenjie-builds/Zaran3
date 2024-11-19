#include "CellFNFDM.h"
#include "Log.h"
namespace zaran
{

CellFN::CellFN( int cell_num):CellBase()
{
    SetCellNum(cell_num);
    m_node_num = new int[cell_num];
    m_node_id = new int[cell_num];
    m_cell_face_num = new int[cell_num];
    m_face_id = new int[cell_num];
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
void CellFN::SetNode(Array<IArray> &node_index)
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
        m_node = new int[node_sum];
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
void CellFN::SetFace(Array<IArray> &face_index)
{
    if (GetCellNum() != face_index.size())
    {
        Log::error("CellTopo::SetFace: The number of cells is not equal to the number of faces");
        Log::error(" cell_num = {}, face_index.size() = {}", GetCellNum(), face_index.size());
        exit(0);
    }
    int face_sum = 0;
    for (int i = 0; i < GetCellNum(); i++)
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
        m_face = new int[face_sum];
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
void CellFN::SetCenter(int iCell, const double *center)
{
	if(iCell >= GetCellNum())
	{
		Log::error("CellTopo::SetCenter: The cell index is out of range");
		Log::error(" iCell = {}, m_cell_num = {}", iCell, GetCellNum());
		exit(0);
	}
    for (int i = 0; i < 3; i++)
    {
        m_center[iCell * 3 + i] = center[i];
    }
}

int CellFN::GetNodeNum(int iCell) const
{
    return m_node_num[iCell];
}
const int *CellFN::GetNode(int iCell)
{
    return m_node + m_node_id[iCell];
}
int CellFN::GetFaceNum(int iCell) const
{
	return m_cell_face_num[iCell];
}
const int *CellFN::GetFace(int iCell)
{
	return m_face + m_face_id[iCell];
}
const double *CellFN::GetCenterCoord(int iCell) const
{
    return m_center + iCell * 3;
}
} // namespace zaran