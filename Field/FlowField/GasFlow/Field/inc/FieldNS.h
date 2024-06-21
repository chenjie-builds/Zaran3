#pragma once
#include "Field.h"
#include"GridFNFDM.h"
#include"NSSolverFNFDM.h"
#include "DataManagerNS.h"
#include "residual.h"
namespace zaran
{

    class FieldNS_FNFDM :public Field
    {
    public:
        FieldNS_FNFDM(GridBase* grid);
        ~FieldNS_FNFDM();
        GridFN* GetGrid() override;
        NSSolverFNFDM* GetSolver() override;
        FlowSolverPara* GetSolverPara() override;
        DataManagerNS* GetDataManager() override;
        ResAnalyzerFN* GetResAnalyzer() { return m_res_analyze; }
    protected:
        void Allocate() override;
        ResAnalyzerFN* m_res_analyze;

    };
}