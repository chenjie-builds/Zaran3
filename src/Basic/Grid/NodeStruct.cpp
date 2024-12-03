#include "NodeStruct.h"
using namespace zaran;
NodeStruct::NodeStruct() 
{
    m_coord = nullptr;
    m_i_num = 0;
    m_j_num = 0;
    m_k_num = 0;

}
NodeStruct::~NodeStruct()
{
    delete[] m_coord;
}

void zaran::NodeStruct::Allocate(int i_num, int j_num, int k_num)
{
    m_i_num = i_num;
    m_j_num = j_num;
    m_k_num = k_num;
    if(m_coord != nullptr)
    {
        delete[] m_coord;
    }
    m_coord = new double[3 * i_num * j_num * k_num];
    SetCount(i_num * j_num * k_num);
}

void NodeStruct::SetCoord(int idx_i, int idx_j, int idx_k, const double *coord)
{
    int index = 3 * GetIdx(idx_i, idx_j, idx_k);
    for(int i_dim = 0; i_dim < 3; i_dim++)
    { 
        m_coord[index + i_dim] = coord[i_dim];
    }
}
const double *NodeStruct::GetCoord(int idx_i, int idx_j, int idx_k) const
{
    return m_coord + 3 * GetIdx(idx_i, idx_j, idx_k);
}
const double* NodeStruct::GetCoord(IdxStruct* idx_proxy) const
{
	return m_coord + 3 * idx_proxy->GetIdx();
}
int NodeStruct::GetIdx(int idx_i, int idx_j, int idx_k) const
{
    return idx_i + idx_j * m_i_num + idx_k * m_i_num * m_j_num;
}
