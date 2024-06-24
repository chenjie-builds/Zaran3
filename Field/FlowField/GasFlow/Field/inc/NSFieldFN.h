#pragma once
#include "NSField.h"
#include"GridFNFDM.h"
#include"NSSolverFNFDM.h"
#include "DataManagerNS.h"
#include "residual.h"
namespace zaran
{

    class NSFieldFNFDM :public FieldNS
    {
    public:
        NSFieldFNFDM(GridBase* grid);
        ~NSFieldFNFDM();
        GridFN* GetGrid() override;
        NSSolverFNFDM* GetSolver() override;
        void CalcResidual() override;
    protected:
        void Allocate() override;
    };
}