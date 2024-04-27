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
#include "grid.h"
#include "solverpara.h"
#include "GlobalData.h"
#include "GridList.h"
namespace zaran
{

	//网格指针
	class Solver
	{
	public:
		Solver();
		virtual ~Solver();
	public:
		void SetIndex(const int& index) { index_ = index; }
		void SetName(const string& name) { name_ = name; }
		void SetPara(SolverPara* para) { para_ = para; }
		void SetGrid(Grid* grid) { m_grid = grid; }
		const string& GetName()const { return name_; }
		virtual Grid* GetGrid() { return m_grid; }
		const int& GetIndex()const { return index_; }
		virtual SolverPara* GetPara();
	public:
		virtual void Init() = 0;
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
		Grid* m_grid;
		// solver 的参数
		SolverPara* para_;
	};
}