#pragma once
#include "NSSolverStruct.h"
#include "Metric.h"
namespace zaran
{
  // 使用CMM方法求解结构网格流动问题
  // 该方法度量系数存储在各个半点上
  // 计算半点数值通量时，直接使用半点的度量系数
  class NSSolverStructCMM : public NSSolverStruct
  {
  public:
    NSSolverStructCMM(Id index, string name, std::shared_ptr<FlowSolverParam> para, std::shared_ptr < GridStruct> grid, std::shared_ptr < DataManagerNSStruct>data_manager);
    ~NSSolverStructCMM();

  protected:
    //==================================================================
    //! vvvv 以下是计算度量系数的函数 vvvv
    void CalcCoordTransCoef() override;
    //! @brief插值得到半点处的度量系数
    void CalcMidNodeMetrics();
    void CalcMidNodeMetrics2nd();
    void CalcMidNodeMetrics6th();
    //! ^^^^ 以上是计算度量系数的函数 ^^^^
    //==================================================================
    void CalcConvectionRes_1st() override;
    void CalcConvectionRes_MUSCL() override;
    void CalcConvectionRes_WCNS5() override;
    //@brief 根据半点左右值计算数值通量
    void CalcConvectionFluxMidNode();

    //@brief 通量差分：二阶中心差分
    void FluxDifference2nd() override;
    //@biref 通量差分：六阶中心差分
    void FluxDifference6th() override;
  };
} // namespace zaran