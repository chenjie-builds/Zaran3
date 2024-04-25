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
#include "SolverVec.h"
#include "solver.h"
#include "flowsolver.h"
#include "NSSolver.h"
#include "Solver_NS_3D.h"
#include"Solver_NS_2D.h"
#include "Solver_NS_3D_Struct.h"
#include "Solver_NS_2D_Struct.h"
#include "Field.h"
namespace zaran
{
    class FieldFactory
    {
    public:
        FieldFactory(GridType grid_type, FieldSolverType solver_type,Dimension dim) :m_grid_type(grid_type), m_solver_type(solver_type),m_dim(dim) {};
        void Create();
       Field** GetField() { return m_field; }
    private:
        void CreateGrid();
        void CreateField();
        void CreateSolver();
    private:
        Field** m_field;
        GridType m_grid_type;
        FieldSolverType m_solver_type;
        Dimension m_dim;

    };
}