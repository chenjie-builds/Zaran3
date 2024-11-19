#include "StructIdxProxy.h"
namespace zaran
{
    StructIdxProxy::StructIdxProxy(int ni, int nj, int nk)
    {
        m_ni = ni;
        m_nj = nj;
        m_nk = nk;
    }

    void StructIdxProxy::SetIdx(int i, int j, int k)
    {
        m_idx_i = i;
        m_idx_j = j;
        m_idx_k = k;
        IdxProxy::SetIdx(i + j * m_ni + k * m_ni * m_nj);
    }
    void StructIdxProxy::SetIdx(int idx)
    {
        m_idx_i = idx % m_ni;
        m_idx_j = (idx / m_ni) % m_nj;
        m_idx_k = idx / (m_ni * m_nj);
        IdxProxy::SetIdx(idx);
    }
    void StructIdxProxy::SetI(int i)
    {
        m_idx_i = i;
        IdxProxy::SetIdx(i + m_idx_j * m_ni + m_idx_k * m_ni * m_nj);
    }
    void StructIdxProxy::SetJ(int j)
    {
        m_idx_j = j;
        IdxProxy::SetIdx(m_idx_i + j * m_ni + m_idx_k * m_ni * m_nj);
    }
    void StructIdxProxy::SetK(int k)
    {
        m_idx_k = k;
        IdxProxy::SetIdx(m_idx_i + m_idx_j * m_ni + k * m_ni * m_nj);
    }
    StructIdxProxy::~StructIdxProxy()
    {
    }
    int StructIdxProxy::GetNi() const
    {
        return m_ni;
    }
    int StructIdxProxy::GetNj() const
    {
        return m_nj;
    }
    int StructIdxProxy::GetNk() const
    {
        return m_nk;
    }
    void StructIdxProxy::GetIdxStruct(int &i, int &j, int &k) 
    {
        i = m_idx_i;
        j = m_idx_j;
        k = m_idx_k;
    }
}