//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Controller.h														||
//*	@brief	Control solver to modify field data in grid							||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include"BasicType.h"
#include "SolverVec.h"
#include "GridList.h"
#include "Visual.h"
#include <iostream>
namespace zaran
{
	class Controller
	{
	public:
		Controller(Ptr<GridList>& gridList, Ptr<SolverVec>& solverVec);
		~Controller();
	public:
		// 流场求解
		void SolveField();
	protected:
		// 前处理
		void PreSolve();
		// 计算一步
		void SolveFieldOneStep();
		// 后处理
		void PostSolve();
	protected:
		//初始化，包括求解器初始化和流场初始化
		void Initialize();
		//计算流场中的最大和最小残差
		double CalcMaxAveResidual();
		// 计算残差
		void CalcResidual();
		// 交互插值节点信息
		void CommInterNodeData();
		// 是否停止计算
		bool IsStopSolve();
	protected:
		//备份流场数据
		void SaveFieldData();
		// 备份边界节点数据
		void SaveWallNode();
		// 备份残差
		void SaveResidual();
		// 输出流场为VTK格式
		void SaveDataVTK(std::ostream& os);
		// 输出流场为Tecplot格式
		void SaveDataTecplot();
	private:
		Ptr<GridList> gridList_;
		Ptr<SolverVec> solverVec_;
		Ptr<Visual> visual_;
		double maxResidual_;
		double aveResidual_;
	};

}