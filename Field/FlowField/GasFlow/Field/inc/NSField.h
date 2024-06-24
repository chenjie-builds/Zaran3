#pragma once
#include "Field.h"
#include "DataManagerNS.h"
#include "residual.h"
#include "FlowSolverPara.h"
namespace zaran
{

    class FieldNS :public Field
    {
    public:
        FieldNS(GridBase* grid, FieldType fieldType);
        ~FieldNS();
        FlowSolverPara* GetSolverPara() override;
        DataManagerNS* GetDataManager() override;
        ResInfo* GetResInfo() { return m_res_info; }
        virtual void CalcResidual() = 0;
    protected:
        void Allocate() override;
        ResInfo* m_res_info;
    };
}