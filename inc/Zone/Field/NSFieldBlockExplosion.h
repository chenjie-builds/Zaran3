#pragma once
#include "NSFieldStruct.h"
#include"NSSolverBlock.h"
#include "GridBlock.h"
#include "FlowSolverStructPara.h"
#include "DataManagerNSStruct.h"
#include "StructIdxProxy.h"
#include "GlobalField.h"
namespace zaran
{

    class NSFieldBlockExplosion : public NSFieldStruct
    {
    public:
        NSFieldBlockExplosion(shared_ptr<GridBlock>grid);
        ~NSFieldBlockExplosion();
        void Allocate() override;
        shared_ptr<GridBlock> GetGrid();
        shared_ptr<NSSolverBlock> GetSolver();
    protected:
        void AllocateSolver() override;
    };
} // namespace zaran