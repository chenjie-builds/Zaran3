#include "StructIdxProxy.h"
namespace zaran
{
    IdxStruct::IdxStruct(int ni, int nj, int nk)
    {
        m_ni = ni;
        m_nj = nj;
        m_nk = nk;
    }

    void IdxStruct::SetIdx(int i, int j, int k)
    {
        m_idx_i = i;
        m_idx_j = j;
        m_idx_k = k;
        Idx::SetIdx(i + j * m_ni + k * m_ni * m_nj);
    }
    void IdxStruct::SetIdx(int idx)
    {
        m_idx_i = idx % m_ni;
        m_idx_j = (idx / m_ni) % m_nj;
        m_idx_k = idx / (m_ni * m_nj);
        Idx::SetIdx(idx);
    }
    void IdxStruct::SetI(int i)
    {
        m_idx_i = i;
        Idx::SetIdx(i + m_idx_j * m_ni + m_idx_k * m_ni * m_nj);
    }
    void IdxStruct::SetJ(int j)
    {
        m_idx_j = j;
        Idx::SetIdx(m_idx_i + j * m_ni + m_idx_k * m_ni * m_nj);
    }
    void IdxStruct::SetK(int k)
    {
        m_idx_k = k;
        Idx::SetIdx(m_idx_i + m_idx_j * m_ni + k * m_ni * m_nj);
    }
    IdxStruct::~IdxStruct()
    {
    }
    int IdxStruct::GetNi() const
    {
        return m_ni;
    }
    int IdxStruct::GetNj() const
    {
        return m_nj;
    }
    int IdxStruct::GetNk() const
    {
        return m_nk;
    }
    void IdxStruct::GetIdxStruct(int &i, int &j, int &k) 
    {
        i = m_idx_i;
        j = m_idx_j;
        k = m_idx_k;
    }
}