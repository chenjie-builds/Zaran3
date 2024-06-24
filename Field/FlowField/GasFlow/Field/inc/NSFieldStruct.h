#pragma once
#include"NSField.h"
#include"GridStruct.h"
#include"NSSolverStruct.h"
#include"DataManagerNS.h"
namespace zaran
{
    class NSFieldStruct :public FieldNS
    {
    public:
        NSFieldStruct(GridBase* grid);
        ~NSFieldStruct();
        GridStruct* GetGrid() override;
        NSSolverStruct* GetSolver() override;
        void CalcResidual() override;
    protected:
        void Allocate() override;
        StructIdxProxy* m_idx_proxy;
    };
}