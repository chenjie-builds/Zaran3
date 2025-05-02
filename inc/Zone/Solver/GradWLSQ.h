/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file GradWLSQ.h
 * \brief GradWLSQ class, used to calculate gradient using Weighted Least Square method.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once

#include "GridFNFDM.h"
namespace zaran
{
    class GradWLSQ
    {
    public:
        GradWLSQ(shared_ptr<GridFN> grid);
        ~GradWLSQ();
        /// @brief Calculate gradient of data using Weighted Least Square method
        /// @param grid  grid information
        /// @param data  data to calculate gradient
        /// @param grad_x  gradient in x direction
        /// @param grad_y  gradient in y direction
        /// @param grad_z  gradient in z direction
        void CalcGradient(shared_ptr<GridFN> grid, const double* data, double* grad_x, double* grad_y, double* grad_z);
    private:
        shared_ptr<GridFN> m_grid;
        double*** m_omega;
    };
}