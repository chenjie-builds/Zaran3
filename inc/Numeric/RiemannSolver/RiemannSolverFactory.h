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
#include <string>
#include <memory>
namespace zaran
{
	using std::string;
	class RiemannSolverFactory
	{
	public:
		void Create(std::shared_ptr<RiemannSolver>& riemannSolver,string& riemannSolverName);
	private:
		void CreateVanLeer(std::shared_ptr<RiemannSolver>& riemannSolver);
		void CreateHLLC(std::shared_ptr<RiemannSolver>& riemannSolver);
		void CreateRoe(std::shared_ptr<RiemannSolver>& riemannSolver);
		void CreateStegerWarming(std::shared_ptr<RiemannSolver>& riemannSolver);
		void CreateAusmpw(std::shared_ptr<RiemannSolver>& riemannSolver);
	};
}