#pragma once
#include "BasicType.h"
namespace zaran
{
    class BoundStruct
    {
    public:
        BoundStruct(int i_bound, int j_bound, int k_bound, int i_ref, int j_ref, int k_ref, int i_ghost, int j_ghost, int k_ghost, const double* norm_bound);
        ~BoundStruct();
        void SetIdxBound(const int& i_bound, const int& j_bound, const int& k_bound);
        void SetIdxRef(const int& i_ref, const int& j_ref, const int& k_ref);
        void SetIdxGhost(const int& i_ghost, const int& j_ghost, const int& k_ghost);
        void SetNormBound(const double* norm_bound);
        void GetIdxBound(int& i, int& j, int& k);
        void GetIdxRef(int& i, int& j, int& k);
        void GetIdxGhost(int& i, int& j, int& k);
        const double* GetNormBound();
        bool operator==(const BoundStruct& bound);
    private:
        int m_i_bound, m_j_bound, m_k_bound;
        int m_i_ref, m_j_ref, m_k_ref;
        int m_i_ghost, m_j_ghost, m_k_ghost;
        double m_norm_bound[3];
    };
}