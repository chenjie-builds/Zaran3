#include "CellStruct.h"
namespace zaran
{
CellStruct::CellStruct(int ni, int nj, int nk) : m_ni(ni), m_nj(nj), m_nk(nk)
{
    SetCellNum(ni * nj * nk);
}
CellStruct::~CellStruct()
{
    if (m_center)
        delete[] m_center;
}
void CellStruct::SetCenterCoord(int i, int j, int k, double *center)
{
    int index = GetIndex(i, j, k) * 3;
    for (int i = 0; i < 3; i++)
    {
        m_center[index + i] = center[i];
    }
}
const double *CellStruct::GetCenterCoord(int i, int j, int k) const
{
    return m_center + 3 * GetIndex(i, j, k);
}
int CellStruct::GetIndex(int i, int j, int k) const
{
    return i + j * m_ni + k * m_ni * m_nj;
}
} // namespace zaran