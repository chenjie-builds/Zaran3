#include "NodeStruct.h"
using namespace zaran;
NodeStruct::NodeStruct()
{
	m_coord.resize(0);
	m_i_num = 0;
	m_j_num = 0;
	m_k_num = 0;

}
NodeStruct::~NodeStruct()
{
}

void zaran::NodeStruct::Allocate(Id i_num, Id j_num, Id k_num)
{
	m_i_num = i_num;
	m_j_num = j_num;
	m_k_num = k_num;
	m_coord.resize(3 * i_num * j_num * k_num);
	SetCount(i_num * j_num * k_num);
}

zaran::Coord NodeStruct::GetCoord(Id idx) const
{
	return &m_coord[3 * idx];
}

void NodeStruct::SetCoord(Id idx_i, Id idx_j, Id idx_k, const double* coord)
{
	Id index = 3 * GetIdx(idx_i, idx_j, idx_k);
	for (int i_dim = 0; i_dim < 3; i_dim++)
	{
		m_coord[index + i_dim] = coord[i_dim];
	}
}
const double* NodeStruct::GetCoord(Id idx_i, Id idx_j, Id idx_k) const
{
	return GetCoord(GetIdx(idx_i, idx_j, idx_k));
}
const double* NodeStruct::GetCoord(IdProxyStruct* idx_proxy) const
{
	return GetCoord(idx_proxy->GetIdx());
}
Id NodeStruct::GetIdx(Id idx_i, Id idx_j, Id idx_k) const
{
	return idx_i + idx_j * m_i_num + idx_k * m_i_num * m_j_num;
}
