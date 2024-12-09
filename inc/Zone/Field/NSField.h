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
        FieldNS(std::shared_ptr<GridBase> grid=nullptr, FieldType fieldType=FieldType::Unset);
        ~FieldNS();
        void Allocate() override;
    public:
        std::shared_ptr<FlowSolverParam> GetSolverPara();
		std::shared_ptr<DataManagerNS> GetDataManager();
        ResInfo* GetResInfo() { return m_res_info; }
        virtual void CalcResidual() = 0;
        virtual void AllocateSolverPara()override;
    protected:
        virtual void AllocateResInfo();
        /// @brief 残差信息，用于存储残差
        ResInfo* m_res_info;
    };
}