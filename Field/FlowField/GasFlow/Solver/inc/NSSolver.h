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
#include "log.h"
#include "MathBasic.h"
#include "CommonPara.h"
#include "RiemannSolverFactory.h"
#include "Limiter.h"
#include "Gas.h"
#include"PerfectGas.h"

namespace zaran
{
	class NSSolver :public FlowSolver
	{
	public:
		NSSolver() :FlowSolver() {}
		~NSSolver() {}
		void Init()override;
		void InitField()override;
		void InitFieldFarFlow();
		void InitFieldFarFieldNoVelocity();
		void InitFieldBackup();
		void InitSolver()override;
		void Solve() override;
		void UpdateField()override;
		void CreateFieldData()override;
		void RegisterFieldData()override;
		void Preprocess()override;
		void Postprocess()override;
		virtual void CalcMetric() = 0;
		double ComputeMaxResidual()override;
		void BackupField(std::string&back_folder)override;
		/*
		* 梯度计算函数
		* 目前包括最小二乘法、非结构有限差分法
		*/
	protected:
		void CalcPrimGrad()override;
		void CalcPrimGradBound();
		// 使用最小二乘求梯度
		virtual void CalcGradWLS() = 0;
		// 使用非结构有限差分法求梯度
		virtual void CalcGradUFDM();
		void NoGradient();
		/**
		 * 时间步计算与时间推进相关函数
		 */
	protected:
		void CalcTimeStep()override;
		void SnycTimeStepWithGlobal(double& dt)override;
		void TimeAdvance()override;
		// 计算当地时间步
		virtual void CalcTimeStepLocal() = 0;
		//龙格库塔法
		virtual void RungeKutta();
		/**
		* 空间流动参数计算
		*/
	protected:
		// 原始变量到守恒变量
		void Prim2Cons();
		void Prim2Cons(double& rho, double& u, double& v, double& w, double& p, double& cons0, double& cons1, double& cons2, double& cons3, double& cons4);
		// 守恒变量到原始变量
		void Cons2Prim();
		void Cons2Prim(double& cons0, double& cons1, double& cons2, double& cons3, double& cons4, double& rho, double& u, double& v, double& w, double& p);
		// 计算流动通量
		virtual void ConvectiveResidual() = 0;
		//计算粘性通量
		virtual void ViscousResidual() = 0;
		//计算源项
		virtual void SourceTermResidual() = 0;
		// 计算流场残差，即右端项
		virtual void CalcResidual();
		// 计算限制器系数
		virtual void CalcLimiter();
		virtual void CalcLimiterVK();
		virtual void CalcLimiterBJ();
		virtual void CalcLimiterNone();
		virtual void CalcLimiterFirstOrder();
		virtual void CalcLimiterBound();
		/// @brief 检查原始变量
		virtual void CheckPrimtive();
		/// @brief 检查残差
		virtual void CheckResidual();
		/// 如果原始变量出现异常值，将其修正为合理值
		virtual void FixPrimtive();
		/// @brief mid point value reconstruct
		/// @param index_left left node index
		/// @param index_right right node index 
		/// @param value_rec_left reconstructed value of left side at mid point
		/// @param value_rec_right reconstructed value of right side at mid point
		void MidPointReconstruct(int index_left, int index_right, double* value_rec_left, double* value_rec_right);
		void MidPointReconstructFirstOrder(int index_left, int index_right, double* value_rec_left, double* value_rec_right);
		virtual	void CalcForce() {};
		double* GetPrimGrad(int iNode, int iVar) { return m_prim_grad + iNode * GetNumberOfEquations() * 3 + iVar * 3; }
		double& GetPrimGrad(int iNode, int iVar, int iDim) { return m_prim_grad[iNode * GetNumberOfEquations() * 3 + iVar * 3 + iDim]; }

	protected:

		void BoundaryCondition()override;
		// 超声速入口边界条件
		virtual void InletBC(Boundary& bound);
		// 超声速出口边界条件
		virtual void OutletBC(Boundary& bound);
		// 壁面边界条件
		virtual void  WallBC(Boundary& bound);
		// 原始变量梯度
		double* m_prim_grad;
		//通量求解器
		Ptr<RiemannSolver> riemannSolver_;
		//降阶标识，对于二阶精度的格式，降阶到一阶
		IArray m_reduce_order;

	};
}