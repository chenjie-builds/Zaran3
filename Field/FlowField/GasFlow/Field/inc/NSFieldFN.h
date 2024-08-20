#pragma once
#include "NSField.h"
#include"GridFNFDM.h"
#include"NSSolverFNFDM.h"
#include "DataManagerNS.h"
#include "Residual.h"
namespace zaran
{

    class NSFieldFNFDM :public FieldNS
    {
    public:
        NSFieldFNFDM(GridBase* grid);
        ~NSFieldFNFDM();
    public:
        GridFN* GetGrid() override;
        NSSolverFNFDM* GetSolver() override;
        void CalcResidual() override;
    protected:
        void AllocateSolver() override;
        void AllocateDataManager() override;
    };
}