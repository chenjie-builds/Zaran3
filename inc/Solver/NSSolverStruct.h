//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	NSSolverStruct.h													||
//*	@brief	结构网格NS方程求解器													||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "NSSolver.h"
#include "StructGrid.h"
namespace zaran
{
	class NSSolverStruct :public NSSolver
	{
	public:
		virtual void ComputeCoordTrans()override;
		// 使用最小二乘求梯度
		virtual void ComputeGradientWLS()override;
	protected:
		Ptr<StructGrid> GetGrid();
		virtual void ComputeTimeStepLocal()override;
		// 计算流动通量
		virtual void InviscidFlux()override;
		// 计算限制器系数
		virtual void ComputeLimiterCoef()override;
	private:
	};
}