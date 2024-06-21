#pragma once
#include"Field.h"
#include"GridStruct.h"
#include"NSSolverStruct.h"
#include"DataManagerNS.h"
namespace zaran
{
    class FieldNS_Struct :public Field
    {
    public:
        FieldNS_Struct(GridBase* grid);
        ~FieldNS_Struct();
        GridStruct* GetGrid() override;
        NSSolverStruct* GetSolver() override;
        FlowSolverPara* GetSolverPara() override;
        DataManagerNS* GetDataManager() override;
    protected:
        void Allocate() override;
        
    };
}