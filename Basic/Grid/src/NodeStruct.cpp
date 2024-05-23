#include "NodeStruct.h"
using namespace zaran;
NodeStruct::NodeStruct(int ni, int nj, int nk) : m_ni(ni), m_nj(nj), m_nk(nk)
{
    m_coord = new double[3 * ni * nj * nk];
    SetNodeNum(ni * nj * nk);
}
NodeStruct::~NodeStruct()
{
    delete[] m_coord;
}

void NodeStruct::SetCoord(int i, int j, int k, const double *coord)
{
    int index = 3 * GetIndex(i, j, k);
    m_coord[index] = coord[0];
    m_coord[index + 1] = coord[1];
    m_coord[index + 2] = coord[2];
}
const double *NodeStruct::GetCoord(int i, int j, int k) const
{
    return m_coord + 3 * GetIndex(i, j, k);
}
int NodeStruct::GetIndex(int i, int j, int k) const
{
    return i + j * m_ni + k * m_ni * m_nj;
}
