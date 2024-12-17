#pragma once
#include <vector>
#include <string>
#include "BasicType.h"
namespace zaran
{
	class FieldDataCommInfo
	{
	public:
		FieldDataCommInfo(index_type recv_node_num, const dynamic_array<std::string>& recv_data_name, const index_type* recv_node_idx_local, const index_type* recv_field_idx_global, const index_type* recv_node_idx_global);
		~FieldDataCommInfo();
		index_type GetRecvNodeNum() { return m_data_idx_src.size(); }
		const dynamic_array<std::string>& GetRecvDataName() const { return m_data_name; }
		/// @brief 获取接收数据缓存
		index_type* GetSrcDataIdx() { return m_data_idx_src.data(); }
		index_type* GetTgtFieldIdx() { return m_field_idx_tgt.data(); }
		index_type* GetTgtNodeIdx() { return m_data_idx_tgt.data(); }
		/// @brief 修改接收到的数据缓存
		/// @param idx_data_name 数据名称编号,即m_recv_data_name的编号
		/// @param idx_recv_node 需要接收数据的node编号,即m_recv_node_idx_src的编号,而不是全局编号
		void UpdateSrcDataBuffer(index_type idx_data_name, index_type idx_recv_node, double data) { m_data_buffer_src[idx_data_name][idx_recv_node] = data; }
		/// @brief 获取接收到的数据缓存
		/// @param idx_data_name 数据名称编号,即m_recv_data_name的编号
		/// @param idx_recv_node 需要接收数据的node编号,即m_recv_node_idx_src的编号,而不是全局编号
		double GetSrcDataBuffer(index_type idx_data_name, index_type idx_recv_node) { return m_data_buffer_src[idx_data_name][idx_recv_node]; }

	private:
		/// @brief 接收数据的名称
		dynamic_array<std::string> m_data_name;
		/// @brief 需要接收的数据编号
		dynamic_array<index_type> m_data_idx_src;
		/// @brief 接收数据来源的field编号
		dynamic_array<index_type> m_field_idx_tgt;
		/// @brief 接收数据来源的node编号
		dynamic_array<index_type> m_data_idx_tgt;
		/// @brief 接收到的数据缓存
		dynamic_array<dynamic_array<double>> m_data_buffer_src;
	};
}