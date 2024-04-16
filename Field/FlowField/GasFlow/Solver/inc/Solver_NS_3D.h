//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Solver_NS_3D.h															||
//*	@brief	三维NS 方程求解器														||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "NSSolver.h"
namespace zaran
{
	class Solver_NS_3D :public NSSolver
	{
	public:
		Solver_NS_3D()  {}
		~Solver_NS_3D() {}
		void CalcMetric()override;
	protected:
		void CalcGradWLS()override;
		void CalcTimeStepLocal() override;
		// 计算流动通量
		void InviscidFlux()override;
		//计算粘性通量
		void ViscousFlux() override;
		//计算源项
		void SourceFlux() override;
		void InviscidFluxNew();
		/// 计算气动力
		void CalcForce();
		

	};
}