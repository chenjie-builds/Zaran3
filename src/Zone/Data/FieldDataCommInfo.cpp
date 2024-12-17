#include "FieldDataCommInfo.h"
namespace zaran
{

	FieldDataCommInfo::FieldDataCommInfo(index_type recv_node_num, const dynamic_array<std::string>& recv_data_name, const index_type* recv_node_idx_local, const index_type* recv_field_idx_global, const index_type* recv_node_idx_global) : m_data_name(recv_data_name)
	{
		m_data_idx_src.resize(recv_node_num);
		m_field_idx_tgt.resize(recv_node_num);
		m_data_idx_tgt.resize(recv_node_num);
		for (int i = 0; i < recv_node_num; i++)
		{
			m_data_idx_src[i] = recv_node_idx_local[i];
			m_field_idx_tgt[i] = recv_field_idx_global[i];
			m_data_idx_tgt[i] = recv_node_idx_global[i];
		}
		m_data_buffer_src.resize(recv_data_name.size());
		for (int i_recv_name = 0; i_recv_name < recv_data_name.size(); i_recv_name++)
		{
			m_data_buffer_src[i_recv_name].resize(recv_node_num);
		}
	}

	FieldDataCommInfo::~FieldDataCommInfo()
	{
	}
}