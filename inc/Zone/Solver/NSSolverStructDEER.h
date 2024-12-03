#pragma once
#include "NSSolverStruct.h"
namespace zaran
{
    // 使用离散等价方程及其离散等价准则求解结构网格流动问题
    // 该方法度量系数存储在节点上
    // 计算半点数值通量时，使用离散点的度量系数
    class NSSolverStructDEER : public NSSolverStruct
    {
    public:
        NSSolverStructDEER(int index, string name, FlowSolverPara *para, GridStruct *grid, DataManagerNSStruct *data_manager);
        ~NSSolverStructDEER();

    protected:
        //==================================================================
        //@brief 通量差分：二阶中心差分
        void FluxDifference2nd() override;
        //@brief 通量差分：二阶中心差分,移动中间节点,用于梯度插值得到中间节点值
        void FlueDifference2nd_MoveMidNode();
		void FluxDifference4th() override;
        //@biref 通量差分：六阶中心差分
        void FluxDifference6th() override;
        /// @brief mid node value interpolate
        /// @param idx_left left node index
        /// @param idx_right right node index
        /// @param lef_coord coordinate of left point
        /// @param mid_coord coordinate of mid point
        /// @param right_coord coordinate of right point
        /// @param value_left interpolated value of left side at mid point
        /// @param value_right interpolated value of right side at mid point
       virtual void MidNodeGrad(int idx_left, int idx_right, const double *lef_coord, const double *mid_coord, const double *right_coord, double *value_left, double *value_right);
    };
}