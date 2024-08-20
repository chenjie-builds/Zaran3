#pragma once
#include "NSSolverStruct.h"
namespace zaran
{
    // 使用离散等价方程及其离散等价准则求解结构网格流动问题
    // 该方法度量系数存储在节点上
    // 计算半点数值通量时，使用离散点的度量系数
    class NSSolverStructDEER :public NSSolverStruct
    {
    public:
        NSSolverStructDEER(int index, string name, FlowSolverPara* para, GridStruct* grid, DataManagerNSStruct* data_manager);
        ~NSSolverStructDEER();
    protected:
        //==================================================================
        //@brief 通量差分：二阶中心差分
        void FluxDifference2nd()override;
        //@biref 通量差分：六阶中心差分
        void FluxDifference6th()override;
    };
}