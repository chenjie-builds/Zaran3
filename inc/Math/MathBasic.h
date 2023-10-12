//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	MathBasic.h															||
//*	@brief	some basic math functions, such as min, max, sign					||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include"BasicType.h"
namespace zaran
{
	double DistanceOfTwoPoints(const double* pt1, const double* pt2);
	template < typename T >
	bool IsErrorData(T data);
	template < typename T >
	T Abs(const T& a) { return (a < 0) ? -a : a; }
	double GetRand(const double& range_low, const double& range_high);
	double Southerland(double T, double mu0 = 1.71608e-5, double T0 = 273.16,
		double Ts = 110.4);
	template<typename T1, typename T2>
	T1 Max(const T1& a, const T2& b) { return a > b ? a : b; }
	template<typename T>
	T Sign(const T& a, const T& b)
	{
		int s = b > 0 ? 1 : -1;
		return s * abs(a);
	}
	template<typename T1, typename T2>
	T1 Min(const T1& a, const T2& b) { return a < b ? a : b; }
	//三点拟合一个圆
	//根据三个二维空间点pt1,pt2,pt3拟合出其外接圆的半径radius,圆心 center
	//https://blog.csdn.net/liyuanbhu/article/details/52891868
	void CircleFrom3Point(const double* pt1, const double* pt2, const double* pt3, double& radius, double* center);
	//向量(x2,y2)逆时针到(x1,y1)的角度，0~2pi，默认(x2,y2)为x轴正方向
	double AngleOfTwoArray(const double& x1, const double& y1, const double& x2 = 1, const double& y2 = 0);
	//三维向量A和向量B的夹角，0~pi
	double AngleOfTwoArray3D(const double* A, const double* B);
}