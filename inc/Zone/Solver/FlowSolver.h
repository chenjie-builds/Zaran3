/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file FlowSolver.h
 * \brief FlowFieldSolver class, used to solve the flow field.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
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