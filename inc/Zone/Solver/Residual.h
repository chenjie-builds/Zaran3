#pragma once
#include"GridFNFDM.h"
namespace zaran
{
    /// @brief 记录残差信息
    /// @details 记录每个方程的最大残差和平均残差，以及最大残差的坐标
    /// @note 并不计算残差，只是记录残差信息，由Field类的派生类计算残差
    class ResInfo
    {
    public:
        ResInfo(index_type equ_num);
        ~ResInfo();
        double GetInfNorm(index_type iVar);
        double GetL2Norm(index_type iVar);
        const double* GetInfNormCoord(index_type Id)const;
        int GetInfNormIdx(index_type iVar)const;
        void SetInfNorm(index_type iVar, double value);
        void SetL2Norm(index_type iVar, double value);
        void SetInfNormCoord(index_type iVar, const double* coord);
        void SetInfNormIdx(index_type iVar, index_type idx);

    private:
        int m_equ_num;
        // 残差的L2范数
		dynamic_array<double> m_ref_L2;
        // 残差的无穷范数
		dynamic_array<double> m_res_inf;
        // 残差的无穷范数的坐标
		dynamic_array<double> m_res_inf_coord;
        // 残差的无穷范数的索引
        // 对于结构网格，还需要通过索引代理算出结构网格的节点编号
		dynamic_array<int> m_res_inf_idx;
    };
}