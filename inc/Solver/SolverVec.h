//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	SolverVec.h															||
//*	@brief	求解器容器类，存储所有求解器											||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "solver.h"
#include <vector>
namespace zaran
{
	using std::shared_ptr;
	//网格的集合，所有的网格均存在此处
	class SolverVec
	{
	public:
		SolverVec();
		~SolverVec();
		//向网格集合中添加一个网格
		void AddSolver(shared_ptr<Solver>& solver);
		//检查网格的index是否和网格集合中匹配
		void CheckSolver();
		//返回一个网格
		shared_ptr<Solver>& GetSolverPtr(const int index);
		//返回网格的个数
		int GetSolverNumber() { return solverVec_.size(); }
	private:
		std::vector<shared_ptr<Solver>> solverVec_;
	};
}