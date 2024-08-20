#pragma once
#include "Field.h"
#include "DataManagerNS.h"
#include "Residual.h"
#include "FlowSolverPara.h"
namespace zaran
{

    class FieldNS :public Field
    {
    public:
        FieldNS(GridBase* grid, FieldType fieldType);
        ~FieldNS();
        void Allocate() override;
    public:
        FlowSolverPara* GetSolverPara() override;
        DataManagerNS* GetDataManager() override;
        ResInfo* GetResInfo() { return m_res_info; }
        virtual void CalcResidual() = 0;
        virtual void AllocateSolverPara()override;
    protected:
        virtual void AllocateResInfo();
        /// @brief 残差信息，用于存储残差
        ResInfo* m_res_info;
    };
}