//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	FlowSolver.h														||
//*	@brief	流场求解器															||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "FieldSolver.h"
#include "flowsolverpara.h"
namespace zaran
{
/// @brief 流场求解器基类
class FlowSolver : public FieldSolver
{
  public:
    FlowSolver(int index, string name, FlowSolverPara *para, GridBase *grid, FieldData *fieldData);
    virtual ~FlowSolver();
  public:
    void InitField() override;
    void InitSolver() override;
    // 计算时间步长
    virtual void CalcTimeStep() = 0;
    void Solve() override;
    void Postprocess() override;
    // 同步时间步长为全局时间步
    virtual void ReduceTimeStep(double &dt) = 0;

  public:
    // 计算当前时刻的CFL数
    double ComputeCFL();

  protected:
    // 返回当前求解器的参数
    FlowSolverPara *GetPara() override;
  public:
    // 时间推进
    virtual void TimeAdvance() = 0;


  protected:
};
} // namespace zaran