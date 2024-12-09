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
        NSFieldBlockExplosion(std::shared_ptr<GridBlock>grid);
        ~NSFieldBlockExplosion();
        void Allocate() override;
        std::shared_ptr<GridBlock> GetGrid();
        std::shared_ptr<NSSolverBlock> GetSolver();
    protected:
        void AllocateSolver() override;
    };
} // namespace zaran