/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file RiemannSolverPara.h
 * \brief Riemann solver parameters, including left and right variable values.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */

#pragma once
#include "Eigen/Dense"
namespace zaran
{
	using namespace Eigen;
	// 通量求解器参数, 包括左右变量值
	struct RiemannSolverPara
	{
		//左右原始变量
		Vector<double ,5> prim_left, prim_right;
		//左右比热比
		double gamma_left, gamma_right;
		//左右能量
		double energy_left, energy_right;
		//左右焓
		double enthalpy_left, enthalpy_right;
		//左右声速
		double c_left, c_right;
		//左右速度大小
		double v2_left, v2_right;
		//法向向量
		Vector3d norm;

		double nt;
		//法向速度
		double vn_left, vn_right;
		//左右马赫数
		double mach_left, mach_right;
		//重构后的通量
		Vector<double, 5> flux;
	};
}