#include "CellStruct.h"
namespace zaran
{
CellStruct::CellStruct() 
{
    m_center_coord = nullptr;
    m_i_num = 0;
    m_j_num = 0;
    m_k_num = 0;
}
CellStruct::~CellStruct()
{
    if (m_center_coord)
        delete[] m_center_coord;
}
void CellStruct::Allocate(int i_num, int j_num, int k_num)
{
    m_i_num = i_num;
    m_j_num = j_num;
    m_k_num = k_num;
    if (m_center_coord)
    {
        delete[] m_center_coord;
    }
    m_center_coord = new double[3 * i_num * j_num * k_num];
}
void CellStruct::SetCenterCoord(int idx_i, int idx_j, int idx_k, double *center)
{
    int index = GetIdx(idx_i, idx_j, idx_k) * 3;
    for (int i_dim = 0; i_dim < 3; i_dim++)
    {
        m_center_coord[index + i_dim] = center[i_dim];
    }
}
const double *CellStruct::GetCenterCoord(int idx_i, int idx_j, int idx_k) const
{
    return m_center_coord + 3 * GetIdx(idx_i, idx_j, idx_k);
}
int CellStruct::GetIdx(int idx_i, int idx_j, int idx_k) const
{
    return idx_i + idx_j * m_i_num + idx_k * m_i_num * m_j_num;
}
} // namespace zaran