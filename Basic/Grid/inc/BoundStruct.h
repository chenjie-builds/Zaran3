#pragma once
#include "BasicType.h"
namespace zaran
{
    class BoundStruct
    {
    public:
        BoundStruct(int i_bound, int j_bound, int k_bound, const int* direction, const double* norm_bound);
        ~BoundStruct();
        void SetIdxBound(const int& i_bound, const int& j_bound, const int& k_bound);
        void SetDirection(const int* direction);
        void SetNormBound(const double* norm_bound);
        void GetIdxBound(int& i, int& j, int& k);
        const int* GetDirection();
        const double* GetNormBound();
        bool operator==(const BoundStruct& bound);
    private:
        int m_i_bound, m_j_bound, m_k_bound;
        // 边界的方向，用于标记边界的方向，只有一个方向为+1或-1，其余为0
        //  -1：左侧边界
        //  +1：右侧边界
        // example: {1,0,0}表示i方向上的右侧边界
        int m_direction[3];
        double m_norm_bound[3];
    };
}