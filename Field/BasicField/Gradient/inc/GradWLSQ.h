//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Gradient.h															||
//*	@brief	根据流场信息，求出对应的梯度    										||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once

#include "GridFNFDM.h"
namespace zaran
{
    class GradWLSQ
    {
    public:
        GradWLSQ(GridFN* grid);
        ~GradWLSQ();
        /// @brief Calculate gradient of data using Weighted Least Square method
        /// @param grid  grid information
        /// @param data  data to calculate gradient
        /// @param grad_x  gradient in x direction
        /// @param grad_y  gradient in y direction
        /// @param grad_z  gradient in z direction
        void CalcGradient(GridFN* grid, const double* data, double* grad_x, double* grad_y, double* grad_z);
    private:
        GridFN* m_grid;
        double*** m_omega;
    };
}