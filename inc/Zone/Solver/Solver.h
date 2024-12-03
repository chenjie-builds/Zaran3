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
namespace zaran
{

	//网格指针
	class Solver
	{
	public:
		Solver(int index, string name, SolverParam* para, GridBase* grid);
		virtual ~Solver();
	public:
		const string& GetName()const { return name_; }
		virtual GridBase* GetGrid() { return m_grid; }
		const int& GetIndex()const { return index_; }
		virtual SolverParam* GetPara();
	public:
		virtual void Init() ;
		virtual void Solve() = 0;
		virtual void Preprocess() = 0;
		virtual void Postprocess() = 0;
		// 初始化自己
		virtual void InitSolver() = 0;
	protected:
		// solver index
		int index_;
		string name_;
		// 全场网格数组
		GridBase* m_grid;
		// solver 的参数
		SolverParam* m_para;
	};
}