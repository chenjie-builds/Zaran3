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
#include"FieldSolver.h"
#include "flowsolverpara.h"
namespace zaran
{
	/// @brief 流场求解器基类
	class FlowSolver :public FieldSolver
	{
	public:
		FlowSolver();
		virtual ~FlowSolver();
	public:
		void InitData() override;
		void InitSolver() override;
		// 计算时间步长
		virtual void CalcTimeStep() = 0;
		void Solve() override;
		void Post()override;
		// 同步时间步长为全局时间步
		virtual void SnycTimeStepWithGlobal(double& dt) = 0;
		virtual	double ComputeMaxResidual() = 0;
	public:
		//计算当前时刻的CFL数
		double ComputeCFL();
	protected:
		//将残差归零
		void ZeroResidual();
		//返回当前求解器的参数
		FlowSolverPara* GetPara();
		//从网格中取数据，为了防止每次直接使用数据名称
	public:
		// 时间推进
		virtual void TimeAdvance() = 0;
		// 计算变量梯度
		virtual void CalcPrimGrad() = 0;
	protected:
		double** m_prim;
		double**m_cons;
		double** m_metric;
		double* m_dt;
		double** m_residual;
		double** m_limiter;
		int* m_non_physical;
	};
}