#include "NodeStruct.h"
using namespace zaran;
NodeStruct::NodeStruct(index_type i_num, index_type j_num, index_type k_num)
{
	m_i_num = i_num;
	m_j_num = j_num;
	m_k_num = k_num;
	m_coord.resize(3 * i_num * j_num * k_num);
	SetCount(i_num * j_num * k_num);
}
NodeStruct::~NodeStruct()
{
}
zaran::Coord NodeStruct::GetCoord(index_type idx) const
{
	return &m_coord[3 * idx];
}

void NodeStruct::SetCoord(index_type idx_i, index_type idx_j, index_type idx_k, const double* coord)
{
	index_type index = 3 * GetIdx(idx_i, idx_j, idx_k);
	for (int i_dim = 0; i_dim < 3; i_dim++)
	{
		m_coord[index + i_dim] = coord[i_dim];
	}
}
const double* NodeStruct::GetCoord(index_type idx_i, index_type idx_j, index_type idx_k) const
{
	return GetCoord(GetIdx(idx_i, idx_j, idx_k));
}
const double* NodeStruct::GetCoord(IdProxyStruct* idx_proxy) const
{
	return GetCoord(idx_proxy->GetIdx());
}
index_type NodeStruct::GetIdx(index_type idx_i, index_type idx_j, index_type idx_k) const
{
	return idx_i + idx_j * m_i_num + idx_k * m_i_num * m_j_num;
}
