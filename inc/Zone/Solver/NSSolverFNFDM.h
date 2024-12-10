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
#include "DataManagerNS.h"
#include "Residual.h"
#include "GradWLSQ.h"
#include "Metric.h"
namespace zaran
{
	class NSSolverFNFDM : public NSSolver
	{
	public:
		NSSolverFNFDM(Id index, string name, std::shared_ptr<FlowSolverParam> para, std::shared_ptr < GridFN> grid, std::shared_ptr < DataManagerNS>data_manager);
		~NSSolverFNFDM();

	protected:
		void InitFieldFarfield() override;
		void InitFieldFarFieldZeroVel() override;
		void InitFieldBackup() override;
		void InitFieldExplosion() override;
		void CalcCoordTransCoef() override;
		void BackupField(std::string& back_folder) override;

	protected:
		void Preprocess() override;
		void CalcTimeStepLocal() override;
		void CalcMinTimeStep(double& dt) override;
		void ReduceTimeStep(double& dt) override;
		void RungeKutta() override;
		void Prim2Cons() override;
		void Cons2Prim() override;

		/// @brief mid point value reconstruct
		/// @param index_left left node index
		/// @param index_right right node index
		/// @param value_rec_left reconstructed value of left side at mid point
		/// @param value_rec_right reconstructed value of right side at mid point
		void MidPointReconstruct2ndOrder(int index_left, int index_right, double* value_rec_left, double* value_rec_right);
		void MidPointReconstruct1stOrder(int index_left, int index_right, double* value_rec_left, double* value_rec_right);

	protected:
		void CalcPrimGrad();
		void CalcPrimGradBound();
		void CalcGradWLS();
		void CalcGradFNFDM();
		void NoGradient();

		// 计算限制器系数
		void CalcLimiter();
		void CalcLimiterVK();
		void CalcLimiterBJ();
		void CalcLimiterNone();
		void CalcLimiterFirstOrder();
		void CalcLimiterBound();
		void CheckPrimtive() override;
		void CheckResidual() override;
		void FixPrimtive() override;

		void ZeroResidual() override;
		// 计算流动通量
		void CalcConvectionResidual() override;
		// 计算粘性通量
		void CalcViscousResidual() override;
		void CalcViscousFlux() override;
		void CalcViscousFluxGrad() override;
		// 计算源项
		void CalcSourceResidual() override;
		/// 计算气动力
		void CalcForce() override;
		void BoundaryCondition() override;
		// 超声速入口边界条件
		virtual void BCInlow(BoundFN& bound);
		// 超声速出口边界条件
		virtual void BCOutflow(BoundFN& bound);
		// 壁面边界条件
		virtual void BCWall(BoundFN& bound);
		// 黎曼边界条件
		virtual void BCFarfield(BoundFN& bound);
		std::shared_ptr<GridFN> GetGrid();
	private:
		std::shared_ptr<GradWLSQ> m_grad_wlsq;
		std::shared_ptr<Metric> m_node_metric;
		std::shared_ptr<GridFN> m_grid;
	};
}