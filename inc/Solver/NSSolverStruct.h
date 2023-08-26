//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	NSSolverStruct.h													||
//*	@brief	结构网格NS方程求解器													||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "NSSolver.h"
#include "StructGrid.h"
namespace zaran
{
	class NSSolverStruct :public NSSolver
	{
	public:
		void InitField()override;
		void InitSolver()override;
		void ComputeTimeStep()override;
		void TimeAdvance()override;
		void BoundaryCondition()override;
		void ComputePrimtiveGradient()override;
	protected:
		using StructGridPtr = std::shared_ptr<StructGrid>;
		StructGridPtr GetGrid();
		void ComputeTimeStepLocal() override;
	private:
	};
}