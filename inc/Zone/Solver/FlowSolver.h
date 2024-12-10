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
#include "FlowSolverPara.h"
#include"FieldDataManager.h"
namespace zaran
{
	/// @brief 流场求解器基类
	class FlowFieldSolver : public FieldSolver
	{
	public:
		FlowFieldSolver(index_type index, string name, shared_ptr<FlowSolverParam> para, shared_ptr < GridBase> grid, shared_ptr < DataManager>data_manager);
		virtual ~FlowFieldSolver();
	public:
		void InitField() override;
		void InitSolver() override;
		// 计算时间步长
		virtual void CalcTimeStep() = 0;
		void Solve() override;
		void Postprocess() override;
		// 同步时间步长为全局时间步
		virtual void ReduceTimeStep(double& dt) = 0;

	public:
		// 计算当前时刻的CFL数
		double ComputeCFL();

	protected:
		// 返回当前求解器的参数
		FlowSolverParam* GetPara();
	public:
		// 时间推进
		virtual void TimeAdvance() = 0;

	private:
		// 流场参数
		shared_ptr<FlowSolverParam> m_para;

	};
} // namespace zaran