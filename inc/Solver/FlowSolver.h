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
	using FlowSolverParaPtr = Ptr<FlowSolverPara>;
	class FlowSolver :public FieldSolver
	{
	public:
		FlowSolver() {};
		virtual ~FlowSolver() {};
	public:
		void InitField() override;
		void InitSolver() override;
		// 计算时间步长
		virtual void ComputeTimeStep() = 0;
		void Solve() override;
		void Post()override;
		// 同步时间步长为全局时间步
		virtual void SnycTimeStepWithGlobal(double& dt) = 0;
	public:
		//计算当前时刻的CFL数
		double ComputeCFL();
	protected:
		//将残差归零
		void ZeroResidual();
		//返回当前求解器的参数
		FlowSolverParaPtr GetPara();
		//从网格中取数据，为了防止每次直接使用数据名称
	public:
		// 时间推进
		virtual void TimeAdvance() = 0;
		// 计算变量梯度
		virtual void ComputePrimtiveGradient() = 0;
	protected:
		Array<DArray*> m_Primtive;
		Array<DArray*> m_Conservative;
		Array<DArray*> m_CoordTrans;
		DArray* m_TimeStep;
		Array<DArray*> m_Residual;
		Array<DArray*> m_LimiterCoef;
		Array<DArray*> m_ConservativeRK;// 用于存储Runge-Kutta的中间量
	};
}