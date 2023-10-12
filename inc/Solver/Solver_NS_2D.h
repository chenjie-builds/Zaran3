//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Solver_NS_2D.h															||
//*	@brief	二维NS 方程求解器														||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "NSSolver.h"
namespace zaran
{
	class Solver_NS_2D :public NSSolver
	{
	public:
		Solver_NS_2D() {}
		~Solver_NS_2D() {}
		void ComputeCoordTrans()override;
		void InitField()override;
	protected:
		void ComputeGradientWLS()override;
		void ComputeTimeStepLocal() override;
		// 计算流动通量
		void InviscidFlux()override;
		void SolveHoleNode();
		void SolveHoleNodeFNFDM();//使用自由节点有限差分法
		void SolveHoleNodeIDW();//逆距离加权法
		void BoundaryCondition()override;
		//计算粘性通量
		void ViscousFlux() override;
		//计算源项
		void SourceFlux() override;
	};
}