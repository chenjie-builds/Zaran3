//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	CommonPara.h														||
//*	@brief	Define some	 constant parameter for global use						||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "BasicType.h"
namespace zaran
{
	//空间维数
	enum class Dimension
	{
		one = 1,
		two = 2,
		three = 3
	};
	const dimension_type one = 1;
	const dimension_type two = 2;
	const dimension_type three = 3;

	const index_type density = 0;
	const index_type velocity_x = 1;
	const index_type velocity_y = 2;
	const index_type velocity_z = 3;
	const index_type pressure = 4;
	const index_type temperature = 5;

	const double LARGE_NUMBER = 1.0e40;
	const double SMALL_NUMBER = 1.0e-40;
	const double TINY_NUMBER = 1.0e-40;
	const double EPSILON_NUMBER = 1.0e-12;

	const double PI = 3.141592653589793238462643383279502884197169399;//圆周率pi
	struct Box
	{
		double x_min, x_max, y_min, y_max, z_min, z_max;
	};
}