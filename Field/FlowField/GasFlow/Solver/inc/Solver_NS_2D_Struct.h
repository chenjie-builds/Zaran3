//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Solver_NS_2D_Struct.h												||
//*	@brief	二维结构网格NS方程求解器												||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "Solver_NS_2D.h"
#include "Grid_Struct_2D.h"
namespace zaran
{
	class Solver_NS_2D_Struct :public Solver_NS_2D
	{
	public:
		Solver_NS_2D_Struct() {}
		~Solver_NS_2D_Struct() {}
		virtual void InitField()override;
		virtual void ComputeCoordTrans()override;
		// 使用最小二乘求梯度
		virtual void ComputeGradientWLS()override;
	protected:
		Ptr<Grid_Struct_2D> GetGrid();
		virtual void ComputeTimeStepLocal()override;
		// 计算流动通量
		virtual void InviscidFlux()override;
		// 使用MUSCL格式计算流动通量
		void InviscidFluxMUSCL();
		// 计算限制器系数
		virtual void ComputeLimiterCoef()override;
		// 超声速入口边界条件
		virtual void ComputeInletBC(Boundary& bound);
		// 超声速出口边界条件
		virtual void ComputeOutletBC(Boundary& bound);
		// 壁面边界条件
		virtual void  ComputeWallBC(Boundary& bound) {};
	private:
		//限制器函数指针
		double (*limiter)(const double&, const double&);
	};
}