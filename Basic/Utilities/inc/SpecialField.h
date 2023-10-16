//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	SpecialField.h														||
//*	@brief	特殊流场生成, 用于测试生成特殊流场										||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include"BasicType.h"
namespace zaran
{
	void CalcSupersonicVortex(const double& x, const double& y, const double& x0, const double y0, const double& ci, const double& Mi, const double& rhoi, DVector& primitive);
	void CalcIsentropicVortex(const double& x, const double& y, const double& beta, DVector& primitive);
}