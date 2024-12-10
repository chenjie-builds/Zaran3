#include "StructIdxProxy.h"
namespace zaran
{
    IdProxyStruct::IdProxyStruct(index_type ni, index_type nj, index_type nk)
    {
        m_ni = ni;
        m_nj = nj;
        m_nk = nk;
    }

    void IdProxyStruct::SetIdx(index_type i, index_type j, index_type k)
    {
        m_idx_i = i;
        m_idx_j = j;
        m_idx_k = k;
        IdProxy::SetIdx(i + j * m_ni + k * m_ni * m_nj);
    }
    void IdProxyStruct::SetIdx(index_type idx)
    {
        m_idx_i = idx % m_ni;
        m_idx_j = (idx / m_ni) % m_nj;
        m_idx_k = idx / (m_ni * m_nj);
        IdProxy::SetIdx(idx);
    }
    IdProxyStruct::~IdProxyStruct()
    {
    }

    void IdProxyStruct::GetIdxStruct(index_type &i, index_type&j, index_type&k) const
    {
        i = m_idx_i;
        j = m_idx_j;
        k = m_idx_k;
    }
}