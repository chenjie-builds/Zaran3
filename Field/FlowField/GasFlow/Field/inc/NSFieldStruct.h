#pragma once
#include "NSField.h"
#include "GridStruct.h"
#include "NSSolverStruct.h"
#include "DataManagerNS.h"
#include "FlowSolverStructPara.h"
namespace zaran
{
    class NSFieldStruct : public FieldNS
    {
    public:
        NSFieldStruct(GridBase *grid);
        ~NSFieldStruct();
        void Allocate() override;

    public:
        GridStruct *GetGrid() override;
        FlowSolverStructPara *GetSolverPara() override;
        NSSolverStruct *GetSolver() override;
        DataManagerNSStruct *GetDataManager() override;
        void CalcResidual() override;
        StructIdxProxy& GetIdxProxy() { return *m_idx_proxy; }
    protected:
        void AllocateResInfo() override;
        void AllocateSolver() override;
        void AllocateDataManager() override;
        void AllocateIdxProxy();
        void AllocateSolverPara() override;
        StructIdxProxy *m_idx_proxy;
    };
}