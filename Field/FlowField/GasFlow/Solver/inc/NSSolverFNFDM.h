//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Solver_NS_3D.h															||
//*	@brief	三维NS 方程求解器														||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "NSSolver.h"
#include "GridFNFDM.h"
#include "FieldDataManagerNS_FNFDM.h"
#include"Residual.h"
#include"GradWLSQ.h"
#include"NodeMetric.h"
namespace zaran
{
	class NSSolverFNFDM :public NSSolver
	{
	public:
		NSSolverFNFDM(int index, string name, FlowSolverPara* para, GridFN* grid, FieldData* fieldData, DataManagerNS_FNFDM* data_manager);
		~NSSolverFNFDM();
	protected:
		void InitFieldFarFlow()override;
		void InitFieldFarFieldNoVelocity()override;
		void InitFieldBackup()override;
		void CalcMetric()override;
		void BackupField(std::string& back_folder)override;
	protected:
		void CalcTimeStepLocal() override;
		void SnycTimeStepWithGlobal(double& dt)override;
		void RungeKutta()override;
		void Prim2Cons()override;
		void Cons2Prim()override;

		/// @brief mid point value reconstruct
		/// @param index_left left node index
		/// @param index_right right node index 
		/// @param value_rec_left reconstructed value of left side at mid point
		/// @param value_rec_right reconstructed value of right side at mid point
		void MidPointReconstruct2ndOrder(int index_left, int index_right, double* value_rec_left, double* value_rec_right);
		void MidPointReconstruct1stOrder(int index_left, int index_right, double* value_rec_left, double* value_rec_right);

	protected:
		void CalcPrimGradBound()override;
		void CalcGradWLS()override;

		// 计算限制器系数
		void CalcLimiter()override;
		void CalcLimiterVK();
		void CalcLimiterBJ();
		void CalcLimiterNone();
		void CalcLimiterFirstOrder();
		void CalcLimiterBound();
		void CheckPrimtive()override;
		void CheckResidual()override;
		void FixPrimtive()override;

		void ZeroResidual() override;
		// 计算流动通量
		void ConvectiveResidual()override;
		//计算粘性通量
		void ViscousResidual() override;
		void CalcViscousFlux() override;
		void CalcViscousFluxGrad()override;
		//计算源项
		void SourceTermResidual() override;
		/// 计算气动力
		void CalcForce()override;
		void BoundaryCondition()override;
		// 超声速入口边界条件
		virtual void InletBC(Boundary& bound);
		// 超声速出口边界条件
		virtual void OutletBC(Boundary& bound);
		// 壁面边界条件
		virtual void  WallBC(Boundary& bound);
		// 黎曼边界条件
		virtual void RiemannBC(Boundary& bound);
	protected:
		virtual GridFN* GetGrid();
	private:
		DataManagerNS_FNFDM* m_data_manager;
		GradWLSQ* m_grad_wlsq;
		NodeMetric* m_node_metric;
	};
}