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
		void InitField() override;
		void InitSolver() override;
		// 计算时间步长
		virtual void CalcTimeStep() = 0;
		void Solve() override;
		void Postprocess()override;
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

		//访问数据
	public:
		double* GetLimiter(int iNode) { return m_limiter + iNode * GetNumberOfEquations(); }
		double& GetLimiter(int iNode, int iEqu) { return m_limiter[iNode * GetNumberOfEquations() + iEqu]; }
		double* GetCons(int iNode) { return m_cons + iNode * GetNumberOfEquations(); }
		double& GetCons(int iNode, int iEqu) { return m_cons[iNode * GetNumberOfEquations() + iEqu]; }
		double* GetResidual(int iNode) { return m_residual + iNode * GetNumberOfEquations(); }
		double& GetResidual(int iNode, int iEqu) { return m_residual[iNode * GetNumberOfEquations() + iEqu]; }
		double& GetTemperture(int iNode, int iEqu) { return m_temperture[iNode]; }
		double* GetMetric(int iNode) { return m_metric + iNode * 17; }
		double* GetMetricXi(int iNode) { return m_metric + iNode * 17; }
		double* GetMetricEta(int iNode) { return m_metric + iNode * 17 + 4; }
		double* GetMetricZeta(int iNode) { return m_metric + iNode * 17 + 8; }
		double* GetMetricTau(int iNode) { return m_metric + iNode * 17 + 12; }
		double& GetMetricJacob(int iNode) { return m_metric[iNode * 17 + 16]; }
	protected:
		double** m_prim;
		//守恒变量
		double* m_cons;
		double* m_metric;
		double* m_dt;
		double* m_residual;
		double* m_limiter;
		double* m_temperture;
		int* m_non_physical;
	};
}