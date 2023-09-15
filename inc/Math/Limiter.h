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
		minmod,
		vanleer,
		vanalbada,
		mixminmodvanleer,
		nolimit,
		barth,
		oneorder
	};
	double NoLimiter(const double& x, const double& y);
	double MinMod(const double& x, const double& y);
	double VanLeer(const double& x, const double& y);
	double MixMinModVanLeer(const double& x, const double& y);
	double VanAlbada(const double& x, const double& y);
	double OneOrder(const double& x, const double& y);
	double Barth(const double& x, const double& y);
}