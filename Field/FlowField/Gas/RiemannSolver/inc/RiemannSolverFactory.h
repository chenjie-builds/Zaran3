//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	.h															||
//*	@brief						||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "RiemannSolver.h"
namespace zaran
{
	class RiemannSolverFactory
	{
	public:
		void Create(Ptr<RiemannSolver>& riemannSolver,string& riemannSolverName);
	private:
		void CreateVanLeer(Ptr<RiemannSolver>& riemannSolver);
		void CreateHLLC(Ptr<RiemannSolver>& riemannSolver);
		void CreateRoe(Ptr<RiemannSolver>& riemannSolver);
		void CreateStegerWarming(Ptr<RiemannSolver>& riemannSolver);
		void CreateAusmpw(Ptr<RiemannSolver>& riemannSolver);
	};
}