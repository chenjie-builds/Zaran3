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
		index_type GetRecvNodeNum() { return m_idx_recv_node.size(); }
		const dynamic_array<std::string>& GetRecvDataName() const { return m_recv_data_name; }
		/// @brief 获取接收数据缓存
		index_type* GetIdxRecvNode() { return m_idx_recv_node.data(); }
		index_type* GetIdxSendField() { return m_idx_send_field.data(); }
		index_type* GetIdxSendNode() { return m_idx_send_node.data(); }
		/// @brief 修改接收到的数据缓存
		/// @param idx_data_name 数据名称编号,即m_recv_data_name的编号
		/// @param idx_recv_node 需要接收数据的node编号,即m_recv_node_idx_src的编号,而不是全局编号
		void SetRecvDataCache(int idx_data_name, int idx_recv_node, double data) { m_recv_data_cache[idx_data_name][idx_recv_node] = data; }
		/// @brief 获取接收到的数据缓存
		/// @param idx_data_name 数据名称编号,即m_recv_data_name的编号
		/// @param idx_recv_node 需要接收数据的node编号,即m_recv_node_idx_src的编号,而不是全局编号
		double GetRecvDataCache(int idx_data_name, int idx_recv_node) { return m_recv_data_cache[idx_data_name][idx_recv_node]; }

	private:
		/// @brief 接收数据的名称
		dynamic_array<std::string> m_recv_data_name;
		/// @brief 需要接收的数据编号
		dynamic_array<index_type> m_idx_recv_node;
		/// @brief 接收数据来源的field编号
		dynamic_array<index_type> m_idx_send_field;
		/// @brief 接收数据来源的node编号
		dynamic_array<index_type> m_idx_send_node;
		/// @brief 接收到的数据缓存
		dynamic_array<dynamic_array<double>> m_recv_data_cache;
	};
}