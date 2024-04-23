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
		virtual void InitData()override;
		virtual void CalcMetric()override;
		// 使用最小二乘求梯度
		virtual void CalcGradWLS()override;
	protected:
		Ptr<Grid_Struct_2D> GetGrid();
		virtual void CalcTimeStepLocal()override;
		// 计算流动通量
		virtual void ConvectiveResidual()override;
		// 使用MUSCL格式计算流动通量
		void InviscidFluxMUSCL();
		// 计算限制器系数
		virtual void CalcLimiter()override;
		// 超声速入口边界条件
		virtual void InletBC(Boundary& bound);
		// 超声速出口边界条件
		virtual void OutletBC(Boundary& bound);
		// 壁面边界条件
		virtual void  WallBC(Boundary& bound) {};
	private:
		//限制器函数指针
		double (*limiter)(const double&, const double&);
	};
}