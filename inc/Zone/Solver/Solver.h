/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file Solver.h
 * \brief Solver class, used to control the simulation.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "SolverPara.h"
#include "GlobalData.h"
#include "GridBase.h"
#include <memory>
namespace zaran
{

	class Solver
	{
	public:
		Solver(index_type index, string name, shared_ptr<SolverParam>para, shared_ptr<GridBase>grid);
		virtual ~Solver();
	public:
		const string& GetName()const { return m_name; }
		GridBase* GetGrid() { return m_grid.get(); }
		const index_type& GetIndex()const { return m_id; }
		SolverParam* GetPara() { return m_para.get(); }
	public:
		virtual void Init();
		virtual void Solve() = 0;
		virtual void Preprocess() = 0;
		virtual void Postprocess() = 0;
		// 初始化自己
		virtual void InitSolver() = 0;
	private:
		// solver index
		index_type m_id;
		string m_name;
		shared_ptr<GridBase> m_grid;
		// solver 的参数
		shared_ptr<SolverParam> m_para;
	};
}