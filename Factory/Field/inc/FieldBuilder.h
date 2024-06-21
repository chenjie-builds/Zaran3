//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	SolverFactory.h														||
//*	@brief	求解器工厂类, 生成求解器												||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "BasicType.h"
#include "solver.h"
#include "flowsolver.h"
#include "NSSolver.h"
#include "GlobalField.h"
namespace zaran
{
    class FieldBuilder
    {
    public:
        FieldBuilder(GridType grid_type, FieldSolverType solver_type,Dimension dim) :m_grid_type(grid_type), m_solver_type(solver_type),m_dim(dim) {};
        FieldManager* Create();
    private:
        void CreateGrid();
        void CreateField();
        void CreateSolver();
    private:
        GridType m_grid_type;
        FieldSolverType m_solver_type;
        Dimension m_dim;
    };
}