#include "BoundStruct.h"
namespace zaran
{
    BoundStruct::BoundStruct(int i_bound, int j_bound, int k_bound, int i_ref, int j_ref, int k_ref, int i_ghost, int j_ghost, int k_ghost, const double* norm_bound):m_i_bound(i_bound),m_j_bound(j_bound),m_k_bound(k_bound),m_i_ref(i_ref),m_j_ref(j_ref),m_k_ref(k_ref),m_i_ghost(i_ghost),m_j_ghost(j_ghost),m_k_ghost(k_ghost)
    {
        for(int i=0;i<3;i++)
        {
            m_norm_bound[i] = norm_bound[i];
        }
    }

    BoundStruct::~BoundStruct()
    {
    }

    void BoundStruct::SetIdxBound(const int& i_bound, const int& j_bound, const int& k_bound)
    {
        m_i_bound = i_bound;
        m_j_bound = j_bound;
        m_k_bound = k_bound;
    }

    void BoundStruct::SetIdxRef(const int& i_ref, const int& j_ref, const int& k_ref)
    {
        m_i_ref = i_ref;
        m_j_ref = j_ref;
        m_k_ref = k_ref;
    }

    void BoundStruct::SetIdxGhost(const int& i_ghost, const int& j_ghost, const int& k_ghost)
    {
        m_i_ghost = i_ghost;
        m_j_ghost = j_ghost;
        m_k_ghost = k_ghost;
    }

    void BoundStruct::SetNormBound(const double* norm_bound)
    {
        for(int i=0;i<3;i++)
        {
            m_norm_bound[i] = norm_bound[i];
        }
    }

    void BoundStruct::GetIdxBound(int& i, int& j, int& k)
    {
        i = m_i_bound;
        j = m_j_bound;
        k = m_k_bound;
    }
    void BoundStruct::GetIdxRef(int& i, int& j, int& k)
    {
        i = m_i_ref;
        j = m_j_ref;
        k = m_k_ref;
    }
    void BoundStruct::GetIdxGhost(int& i, int& j, int& k)
    {
        i = m_i_ghost;
        j = m_j_ghost;
        k = m_k_ghost;
    }
const double* BoundStruct::GetNormBound()
    {
        return m_norm_bound;
    }
bool BoundStruct::operator==(const BoundStruct& bound)
{
    return m_i_bound == bound.m_i_bound && m_j_bound == bound.m_j_bound && m_k_bound == bound.m_k_bound && m_i_ref == bound.m_i_ref && m_j_ref == bound.m_j_ref && m_k_ref == bound.m_k_ref && m_i_ghost == bound.m_i_ghost && m_j_ghost == bound.m_j_ghost && m_k_ghost == bound.m_k_ghost && m_norm_bound == bound.m_norm_bound;
}
}