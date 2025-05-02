/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file LastSquare.h
 * \brief Least square solver.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include"BasicType.h"
namespace zaran
{
	//@brief 使用最小二乘法求解Ax=b
	class LastSquare
	{
	public:
		LastSquare(const Eigen::MatrixXd& A, const  Eigen::VectorXd& b, const int& solver_type = 0) :A_(A), b_(b), solver_type_(solver_type) {};
		Eigen::VectorXd Solver();
	private:
		Eigen::MatrixXd A_;
		Eigen::VectorXd b_;
		Eigen::VectorXd x_;
		int solver_type_;
	private:
		//几种最小二乘方法 https://blog.csdn.net/weixin_46581517/article/details/105178304
		void SolverNormalMatrix();
		void SolverSVD();
		void SolverColPivHouseholderQR();
	};
}