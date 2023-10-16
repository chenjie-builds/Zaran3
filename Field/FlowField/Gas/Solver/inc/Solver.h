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
	using GridPtr =Ptr<Grid>;
	using GridListPtr = Ptr<GridList>;
	using SolverParaPtr = Ptr<SolverPara>;
	class Solver
	{
	public:
		Solver();
		virtual ~Solver();
	public:
		void SetIndex(const int& index) { index_ = index; }
		void SetName(const string& name) { name_ = name; }
		void SetPara(SolverParaPtr para) { para_ = para; }
		void SetGridList(GridListPtr& gridList) { gridList_ = gridList; }
		void SetGridIndex(const int& gridIndex) { gridIndex_ = gridIndex; }
		const string& GetName()const { return name_; }
		GridPtr& GetGrid() { return gridList_->GetGrid(gridIndex_); }
		const int& GetIndex()const { return index_; }
		SolverParaPtr& GetPara();
	public:
		virtual void Init() = 0;
		virtual void Solve() = 0;
		virtual void Post() = 0;

		
		// 初始化自己
		virtual void InitSolver() = 0;


	protected:
		// solver index
		int index_;
		//! solver name
		string name_;
		// 全场网格数组
		GridListPtr gridList_;
		//! solver 作用的网格id
		int gridIndex_;
		// solver 的参数
		SolverParaPtr para_;
	};
}