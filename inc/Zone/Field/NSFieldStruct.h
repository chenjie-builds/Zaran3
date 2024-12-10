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
        NSFieldStruct(shared_ptr<GridStruct> grid);
        ~NSFieldStruct();
        void Allocate() override;

    public:
        shared_ptr<GridStruct>GetGrid();
		shared_ptr<FlowSolverParamStruct>GetSolverPara();
		shared_ptr<NSSolverStruct>GetSolver();
		shared_ptr<DataManagerNSStruct>GetDataManager();
        void CalcResidual() override;
        IdProxyStruct& GetIdxProxy() { return *(GetGrid()->GetIdxProxy()); }
    protected:
        void AllocateResInfo() override;
        void AllocateSolver() override;
        void AllocateDataManager() override;
        void AllocateSolverPara() override;
    };
}