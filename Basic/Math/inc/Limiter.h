//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Limiter.h															||
//*	@brief	限制器函数															||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
namespace zaran
{
	enum class LimiterType//限制器类型
	{
		none,//无限制器
		barth,//Barth限制器
		first_order,//一阶限制器
		vk//Venkatakrishnan限制器
	};
	double LimiterNone(const double& x, const double& y);
	double LimiterMinMod(const double& x, const double& y);
	double LimiterVanLeer(const double& x, const double& y);
	double LimiterMixMinModVanLeer(const double& x, const double& y);
	double LimiterVanAlbada(const double& x, const double& y);
	double LimiterFirstOrder(const double& x, const double& y);
	double LimiterBarth(const double& x, const double& y);
	double LimiterVK(const double& x, const double& y, const double& eps);

}