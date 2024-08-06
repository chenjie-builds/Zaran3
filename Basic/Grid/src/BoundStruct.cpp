#include "BoundStruct.h"
namespace zaran
{


    BoundStruct::BoundStruct(int i_bound, int j_bound, int k_bound, const int* direction, const double* norm_bound) :m_i_bound(i_bound), m_j_bound(j_bound), m_k_bound(k_bound)
    {
        for (int i = 0;i < 3;i++)
        {
            m_direction[i] = direction[i];
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

    void BoundStruct::SetDirection(const int* direction)
    {
        for (int i = 0;i < 3;i++)
        {
            m_direction[i] = direction[i];
        }
    }




    void BoundStruct::SetNormBound(const double* norm_bound)
    {
        for (int i = 0;i < 3;i++)
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
    const int* BoundStruct::GetDirection()
    {
        return m_direction;
    }
    const double* BoundStruct::GetNormBound()
    {
        return m_norm_bound;
    }
    bool BoundStruct::operator==(const BoundStruct& bound)
    {
        return m_i_bound == bound.m_i_bound && m_j_bound == bound.m_j_bound && m_k_bound == bound.m_k_bound && m_direction == bound.m_direction && m_norm_bound == bound.m_norm_bound;
    }
}