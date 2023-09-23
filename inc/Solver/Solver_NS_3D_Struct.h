//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Solver_NS_3D_Struct.h												||
//*	@brief	三维结构网格NS方程求解器												||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "Solver_NS_3D.h"
#include "Grid_Struct_3D.h"
namespace zaran
{
	class Solver_NS_3D_Struct :public Solver_NS_3D
	{
	public:
		Solver_NS_3D_Struct() {}
		~Solver_NS_3D_Struct() {}
		virtual void ComputeCoordTrans()override;
		// 使用最小二乘求梯度
		virtual void ComputeGradientWLS()override;
	protected:
		Ptr<Grid_Struct_3D> GetGrid();
		virtual void ComputeTimeStepLocal()override;
		// 计算流动通量
		virtual void InviscidFlux()override;
		// 计算限制器系数
		virtual void ComputeLimiterCoef()override;
	private:
	};
}