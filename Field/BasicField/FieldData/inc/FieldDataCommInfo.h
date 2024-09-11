#pragma once
#include <vector>
#include <string>
namespace zaran
{
    class FieldDataCommInfo
    {
    public:
        FieldDataCommInfo(int recv_per_node,int recv_node_num, const std::vector<std::string> &recv_data_name, const int *recv_node_idx_local, const int *recv_field_idx_global, const int *recv_node_idx_global);
        ~FieldDataCommInfo();
        int GetRecvNodeNum() { return m_recv_data_name.size(); }
        const std::vector<std::string> &GetRecvDataName() const { return m_recv_data_name; }
        int GetRecvPerNode() { return m_recv_per_node; }
        /// @brief 获取接收数据缓存
        double *GetRecvDataCache() { return m_recv_data_cache; }
        int *GetRecvNodeIdxLocal() { return m_recv_node_idx_local; }
        int *GetRecvFieldIdxGlobal() { return m_recv_field_idx_global; }
        int *GetRecvNodeIdxGlobal() { return m_recv_node_idx_global; }

    private:
        /// @brief 每个node需要发送的数据数量
        int m_recv_per_node;
        /// @brief 接收数据的数量
        std::vector<std::string> m_recv_data_name;
        /// @brief 接收数据的node数量
        int m_recv_node_num;
        /// @brief 需要接收的数据编号
        int *m_recv_node_idx_local;
        /// @brief 接收数据来源的field编号
        int *m_recv_field_idx_global;
        /// @brief 接收数据来源的node编号
        int *m_recv_node_idx_global;
        /// @brief 接收到的数据缓存
        double *m_recv_data_cache;
    };
}