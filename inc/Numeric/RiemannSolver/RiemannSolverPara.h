//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	RiemannSolverPara.h													||
//*	@brief	黎曼求解器参数, RiemannSolver计算所需参数, 用于防止每次新建临时变量		||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include <Eigen/Dense>
namespace zaran
{
	using namespace Eigen;
	// 通量求解器参数, 包括左右变量值
	struct RiemannSolverPara
	{
		//左右原始变量
		Vector<double ,5> primL, primR;
		//左右比热比
		double gammaL, gammaR;
		//左右能量
		double eL, eR;
		//左右声速
		double cL, cR;
		//左右速度大小
		double v2L, v2R;
		//法向向量
		Vector3d norm;

		double nt;
		//法向速度
		double vnL, vnR;
		//左右马赫数
		double machL, machR;
		//重构后的通量
		Vector<double, 5> flux;
	};
}