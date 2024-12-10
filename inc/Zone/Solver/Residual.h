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
        ResInfo(Id equ_num);
        ~ResInfo();
        double GetInfNorm(Id iVar);
        double GetL2Norm(Id iVar);
        const double* GetInfNormCoord(Id Id)const;
        int GetInfNormIdx(Id iVar)const;
        void SetInfNorm(Id iVar, double value);
        void SetL2Norm(Id iVar, double value);
        void SetInfNormCoord(Id iVar, const double* coord);
        void SetInfNormIdx(Id iVar, Id idx);

    private:
        int m_equ_num;
        // 残差的L2范数
		Array<double> m_ref_L2;
        // 残差的无穷范数
		Array<double> m_res_inf;
        // 残差的无穷范数的坐标
		Array<double> m_res_inf_coord;
        // 残差的无穷范数的索引
        // 对于结构网格，还需要通过索引代理算出结构网格的节点编号
		Array<int> m_res_inf_idx;
    };
}