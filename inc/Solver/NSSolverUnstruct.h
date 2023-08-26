//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	NSSolverUnstruct.h													||
//*	@brief	非结构网格NS方程求解器												||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "NSSolver.h"
#include "UnstructGrid.h"
namespace zaran
{
	class NSSolverUnstruct :public NSSolver
	{
	public:
		void InitField()override;
		void InitSolver()override;
		void Post()override;
		void ComputeTimeStep()override;
		void TimeAdvance()override;
		void BoundaryCondition()override;
		void ComputePrimtiveGradient()override;
	protected:
		using UnstructGridPtr = std::shared_ptr<UnstructGrid>;
		UnstructGridPtr GetGrid();
		void ComputeTimeStepLocal() override;
	private:
	};
}