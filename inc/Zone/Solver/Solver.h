//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Solver.h															||
//*	@brief	求解器基类,通过一些算法修改网格数据										||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "SolverPara.h"
#include "GlobalData.h"
#include "GridList.h"
#include <memory>
namespace zaran
{

	//网格指针
	class Solver
	{
	public:
		Solver(Id index, string name, std::shared_ptr<SolverParam>para, std::shared_ptr<GridBase>grid);
		virtual ~Solver();
	public:
		const string& GetName()const { return m_name; }
		GridBase* GetGrid() { return m_grid.get(); }
		const Id& GetIndex()const { return m_id; }
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
		Id m_id;
		string m_name;
		std::shared_ptr<GridBase> m_grid;
		// solver 的参数
		std::shared_ptr<SolverParam> m_para;
	};
}