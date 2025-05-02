/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file DataManager.h
 * \brief DataManagerNSStruct class, used to manage field data.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include"FieldData.h"
#include "GridBase.h"
namespace zaran
{
    class DataManager
    {
    public:
        DataManager(shared_ptr<FieldData> fieldData, int data_num);
        virtual ~DataManager();
    public:
        /// @brief 创建数据
        virtual void CreateData() = 0;
        /// @brief 注册数据
        virtual void RegisterData() = 0;
    protected:
        /// @brief 场数据
        shared_ptr<FieldData>m_data;
        /// @brief 数据个数
        int m_data_num;
    };
}