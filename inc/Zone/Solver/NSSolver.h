//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	NSSolver.h															||
//*	@brief	NS 方程求解器														||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "FlowSolver.h"
#include "RiemannSolver.h"
#include "Gas.h"
#include "DataManagerNS.h"
namespace zaran
{
	class NSSolver : public FlowFieldSolver
	{
	public:
		NSSolver(int index, string name, FlowSolverPara *para, GridBase *grid, DataManagerNS *data_manager);
		~NSSolver();
		void Init() override;
		void InitField() override;
		DataManagerNS *GetDataManager() override;

	protected:
		virtual void InitFieldFarfield() = 0;
		virtual void InitFieldFarFieldNoVelocity() = 0;
		virtual void InitFieldBackup() = 0;
		virtual void InitFieldExplosion() = 0;
		void InitSolver() override;
		void Solve() override;
		void UpdateField() override;
		void Preprocess() override;
		void Postprocess() override;
		// 计算坐标变换系数，度量系数，雅克比矩阵等
		virtual void CalcCoordTransCoef() = 0;

	protected:
		Gas *GetGas() { return m_gas; }
		/*
		 * 时间步计算与时间推进相关函数
		 */
	protected:
		void CalcTimeStep() override;
		void TimeAdvance() override;
		// 计算当地时间步
		virtual void CalcTimeStepLocal() = 0;
		// 统计全局最小时间步
		virtual void CalcMinTimeStep(double& min_dt) = 0;
		// 龙格库塔法
		virtual void RungeKutta() = 0;
		/**
		 * 空间流动参数计算
		 */
	protected:
		// 原始变量到守恒变量
		virtual void Prim2Cons() = 0;
		// 守恒变量到原始变量
		virtual void Cons2Prim() {}
		virtual void ZeroResidual() = 0;
		// 计算流场残差，即右端项
		virtual void CalcResidual();
		// 计算流动通量残差
		virtual void CalcConvectionResidual() = 0;
		// 计算粘性通量残差
		virtual void CalcViscousResidual() = 0;
		virtual void CalcViscousFlux() = 0;
		virtual void CalcViscousFluxGrad() = 0;
		// 计算源项
		virtual void CalcSourceResidual() = 0;
		/// @brief 检查原始变量
		virtual void CheckPrimtive() = 0;
		/// @brief 检查残差
		virtual void CheckResidual() = 0;
		/// 如果原始变量出现异常值，将其修正为合理值
		virtual void FixPrimtive() = 0;
		virtual void CalcForce() {};

	protected:
		Gas *m_gas;
		// 通量求解器
		RiemannSolver *m_riemann_solver;
	};
}