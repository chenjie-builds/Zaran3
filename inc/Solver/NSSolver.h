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
#include "flowsolver.h"
#include "RiemannSolver.h"
namespace zaran
{
	class NSSolver :public FlowSolver
	{
	public:
		NSSolver() {}
		~NSSolver() {}
		void Init()override;
		void InitField()override;
		void InitSolver()override;
		void Solve() override;
		void UpdateField()override;
		void CreateFieldData()override;
		void RegisterFieldData()override;
		virtual void ComputeCoordTrans();
		double ComputeMaxResidual()override;
		/*
		* 梯度计算函数
		* 目前包括最小二乘法、非结构有限差分法
		*/
	protected:
		void ComputePrimtiveGradient()override;
		// 使用最小二乘求梯度
		virtual void ComputeGradientWLS();
		// 使用非结构有限差分法求梯度
		virtual void ComputeGradientUFDM();
		void NoGradient();
		/**
		 * 时间步计算与时间推进相关函数
		 */
	protected:
		void ComputeTimeStep()override;
		void SnycTimeStepWithGlobal(double& dt)override;
		void TimeAdvance()override;
		// 计算当地时间步
		virtual void ComputeTimeStepLocal();
		//龙格库塔法
		virtual void RungeKutta();
		/**
		* 空间流动参数计算
		*/
	protected:
		// 原始变量到守恒变量
		void Primitive2Conservative();
		void Primitive2Conservative(double& rho, double& u, double& v, double& w, double& p, double& cons0, double& cons1, double& cons2, double& cons3, double& cons4);
		// 守恒变量到原始变量
		void Conservative2Primitive();
		void Conservative2Primitive(double& cons0, double& cons1, double& cons2, double& cons3, double& cons4, double& rho, double& u, double& v, double& w, double& p);
		// 计算流动通量
		virtual void InviscidFlux();
		//计算粘性通量
		virtual void ViscousFlux();
		//计算源项
		virtual void SourceFlux();
		// 计算流场残差，即右端项
		virtual void ComputeResidual();
		// 计算限制器系数
		virtual void ComputeLimiterCoef();
		// 创建限制器
		void CreateLimiter();

		/**
		* 边界条件处理
		*/
	protected:
		void BoundaryCondition()override;
		// 超声速入口边界条件
		void ComputeInletBC(Boundary& bound);
		// 超声速出口边界条件
		void ComputeOutletBC(Boundary& bound);
		// 壁面边界条件
		void  ComputeWallBC(Boundary& bound);
		// 原始变量梯度
		Array<DArray*> m_PrimGradX;
		Array<DArray*> m_PrimGradY;
		Array<DArray*> m_PrimGradZ;
		//通量求解器
		Ptr<RiemannSolver> riemannSolver_;
		//限制器函数指针
		double (*limiter)(const double&, const double&);
	};
}