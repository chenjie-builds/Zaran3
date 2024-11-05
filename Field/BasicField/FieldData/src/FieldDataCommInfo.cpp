#include "FieldDataCommInfo.h"
namespace zaran
{

    FieldDataCommInfo::FieldDataCommInfo(int recv_node_num, const std::vector<std::string> &recv_data_name, const int *recv_node_idx_local, const int *recv_field_idx_global, const int *recv_node_idx_global) : m_recv_node_num(recv_node_num), m_recv_data_name(recv_data_name)
    {
        m_idx_recv_node = new int[recv_node_num];
        m_idx_send_field = new int[recv_node_num];
        m_idx_send_node = new int[recv_node_num];
        for (int i = 0; i < recv_node_num; i++)
        {
            m_idx_recv_node[i] = recv_node_idx_local[i];
            m_idx_send_field[i] = recv_field_idx_global[i];
            m_idx_send_node[i] = recv_node_idx_global[i];
        }
        m_recv_data_cache = new double *[recv_data_name.size()];
        for (int i_recv_name = 0; i_recv_name < recv_data_name.size(); i_recv_name++)
        {
            m_recv_data_cache[i_recv_name] = new double[recv_node_num];
        }
    }

    FieldDataCommInfo::~FieldDataCommInfo()
    {
        delete[] m_idx_recv_node;
        delete[] m_idx_send_field;
        delete[] m_idx_send_node;
        for (int i_recv_name = 0; i_recv_name < m_recv_data_name.size(); i_recv_name++)
        {
            delete[] m_recv_data_cache[i_recv_name];
        }
    }
}