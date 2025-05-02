/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file FieldNS.h
 * \brief NS Field class, used to control the NS field simulation.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
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
        FieldNS(shared_ptr<GridBase> grid=nullptr, FieldType fieldType=FieldType::Unset);
        ~FieldNS();
        void Allocate() override;
    public:
        shared_ptr<FlowSolverParam> GetSolverPara();
		shared_ptr<DataManagerNS> GetDataManager();
		ResInfo* GetResInfo() { return m_res_info.get(); }
        virtual void CalcResidual() = 0;
        virtual void AllocateSolverPara()override;
    protected:
        virtual void AllocateResInfo();
        /// @brief 残差信息，用于存储残差
		shared_ptr<ResInfo> m_res_info;
    };
}