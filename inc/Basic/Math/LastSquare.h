//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	LastSquare.h														||
//*	@brief	最小二乘法求解														||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include"BasicType.h"
namespace zaran
{
	//@brief 使用最小二乘法求解Ax=b
	class LastSquare
	{
	public:
		LastSquare(const Matrix& A, const DVector& b, const int& solver_type = 0) :A_(A), b_(b), solver_type_(solver_type) {};
		DVector Solver();
	private:
		Matrix A_;
		DVector b_;
		DVector x_;
		int solver_type_;
	private:
		//几种最小二乘方法 https://blog.csdn.net/weixin_46581517/article/details/105178304
		void SolverNormalMatrix();
		void SolverSVD();
		void SolverColPivHouseholderQR();
	};
}