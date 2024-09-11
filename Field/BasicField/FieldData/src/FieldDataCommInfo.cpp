#include "FieldDataCommInfo.h"
namespace zaran
{

    FieldDataCommInfo::FieldDataCommInfo(int recv_per_node, int recv_node_num, const std::vector<std::string> &recv_data_name, const int *recv_node_idx_local, const int *recv_field_idx_global, const int *recv_node_idx_global) : m_recv_per_node(0), m_recv_data_name(std::vector<std::string>())
    {
        m_recv_per_node = recv_per_node;
        m_recv_node_num = recv_node_num;
        m_recv_data_name = recv_data_name;
        m_recv_node_idx_local = new int[recv_node_num];
        m_recv_field_idx_global = new int[recv_node_num];
        m_recv_node_idx_global = new int[recv_node_num];
        for (int i = 0; i < recv_node_num; i++)
        {
            m_recv_node_idx_local[i] = recv_node_idx_local[i];
            m_recv_field_idx_global[i] = recv_field_idx_global[i];
            m_recv_node_idx_global[i] = recv_node_idx_global[i];
        }
        m_recv_data_cache = new double[recv_per_node * recv_node_num];
    }
    FieldDataCommInfo::~FieldDataCommInfo()
    {
        delete[] m_recv_node_idx_local;
        delete[] m_recv_field_idx_global;
        delete[] m_recv_node_idx_global;
        delete[] m_recv_data_cache;
    }
}