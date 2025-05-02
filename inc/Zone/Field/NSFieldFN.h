/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file NSFieldFNFDM.h
 * \brief NS Field class, used to control the NS field simulation in flexible node grid.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "NSField.h"
#include "GridFNFDM.h"
#include "NSSolverFNFDM.h"
#include "DataManagerNS.h"
#include "Residual.h"
namespace zaran
{

    class NSFieldFNFDM :public FieldNS
    {
    public:
        NSFieldFNFDM(shared_ptr<GridFN> grid);
        ~NSFieldFNFDM();
    public:
		shared_ptr<GridFN> GetGrid();
		shared_ptr<NSSolverFNFDM> GetSolver();
        void CalcResidual() override;
    protected:
        void AllocateSolver() override;
        void AllocateDataManager() override;
    };
}