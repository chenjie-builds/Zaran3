#include "FieldDataCommInfo.h"
namespace zaran
{

	FieldDataCommInfo::FieldDataCommInfo(index_type recv_node_num, const dynamic_array<std::string>& recv_data_name, const index_type* recv_node_idx_local, const index_type* recv_field_idx_global, const index_type* recv_node_idx_global) : m_recv_data_name(recv_data_name)
	{
		m_idx_recv_node.resize(recv_node_num);
		m_idx_send_field.resize(recv_node_num);
		m_idx_send_node.resize(recv_node_num);
		for (int i = 0; i < recv_node_num; i++)
		{
			m_idx_recv_node[i] = recv_node_idx_local[i];
			m_idx_send_field[i] = recv_field_idx_global[i];
			m_idx_send_node[i] = recv_node_idx_global[i];
		}
		m_recv_data_cache.resize(recv_data_name.size());
		for (int i_recv_name = 0; i_recv_name < recv_data_name.size(); i_recv_name++)
		{
			m_recv_data_cache[i_recv_name].resize(recv_node_num);
		}
	}

	FieldDataCommInfo::~FieldDataCommInfo()
	{
	}
}