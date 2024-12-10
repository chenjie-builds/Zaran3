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