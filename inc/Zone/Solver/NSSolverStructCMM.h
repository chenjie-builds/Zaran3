/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file NSSolverStructCMM.h
 * \brief NSSolverStruct class, used to solve the NS equations on structured grid, using CMM method.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
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
    NSSolverStructCMM(index_type index, string name, shared_ptr<FlowSolverParamStruct> para, shared_ptr < GridStruct> grid, shared_ptr < DataManagerNSStruct>data_manager);
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
	void CalcConvectionRes_Grad() override;


    //@brief 根据半点左右值计算数值通量
    void CalcConvectionFluxMidNode();

    //@brief 通量差分：二阶中心差分
    void FluxDifference2nd() override;
    //@biref 通量差分：六阶中心差分
    void FluxDifference6th() override;
  };
} // namespace zaran