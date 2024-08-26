//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Roe.h																||
//*	@brief	Roe ÀèÂüÇó½âÆ÷														||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "RiemannSolver.h"
namespace zaran
{
	class Roe :public RiemannSolver
	{
	public:
		void Solver(RiemannSolverPara& para)override;
	private:
		double a1, a2, a3, a4, a5, a6, a7, a8;
		double delta_rho, delta_u, delta_v, delta_w, delta_p;
		double delta_vn;
		double rho_roe, u_roe, v_roe, w_roe;
		double energy_roe;
		double enthalpy_roe;
		double c_roe;
		double vn_roe;
		Vector<double, 5> flux_roe;

	
	};
}