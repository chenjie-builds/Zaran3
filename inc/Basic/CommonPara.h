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
const double gama = 1.4;
//空间维数
enum class Dimension
{
	one = 1,
	two = 2,
	three = 3
};

const double LARGE_NUMBER = 1.0e40;
const double SMALL_NUMBER = 1.0e-40;
const double TINY_NUMBER = 1.0e-40;
const double EPSILON_NUMBER = 1.0e-12;

const double PI = 3.141592653589793238462643383279502884197169399;//圆周率pi



