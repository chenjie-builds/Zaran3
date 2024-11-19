#pragma once
#include"FieldData.h"
#include "GridBase.h"
namespace zaran
{
    class DataManager
    {
    public:
        DataManager(FieldData* fieldData, int data_num);
        virtual ~DataManager();
    public:
        /// @brief 创建数据
        virtual void CreateData() = 0;
        /// @brief 注册数据
        virtual void RegisterData() = 0;
    protected:
        /// @brief 场数据
        FieldData* m_data;
        /// @brief 数据个数
        int m_data_num;
    };
}