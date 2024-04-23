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
		void CalcMetric()override;
		void InitData()override;
		void InitFieldFarField();//初始化为来流
		void InitFieldNoFlow();//初始化为静止流场
		void InitFieldShockReflection();//初始化为激波反射
		void InitFieldIsentropicVortex();
	protected:
		void CalcGradWLS()override;
		void CalcTimeStepLocal() override;
		// 计算流动通量
		void ConvectiveResidual()override;
		void SolveHoleNode();
		void SolveHoleNodeFNFDM();//使用自由节点有限差分法
		void SolveHoleNodeIDW();//逆距离加权法
		//自定义边界
		void SolveUserDefinedBoundary();

		void BoundaryCondition()override;
		//计算粘性通量
		void ViscousResidual() override;
		//计算源项
		void SourceTermResidual() override;
	};
}