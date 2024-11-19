#pragma once
#include "Field.h"
#include "FieldDataCommInfo.h"
#include <vector>
namespace zaran
{
    /// @brief Field管理类
    /// @details 用于管理全部的Field，提供Field的添加、删除、获取等功能
    class FieldManager
    {
    public:
        FieldManager();
        ~FieldManager();

    public:
        /// @brief 释放
        void FreeMemory();
        /// @brief 获取Field
        Field *GetField(int idx_field);
        /// @brief 获取field data communication information
        FieldDataCommInfo *GetFieldDataCommInfo(int idx_field);
        void SetFieldDataCommInfo(int idx_field, FieldDataCommInfo *field_data_comm_info);
        /// @brief 获取Field数量
        int GetFieldNum();
        /// @brief 添加Field
        void AddField(Field *field, FieldDataCommInfo *field_data_comm_info);
        /// @brief 移除Field
        void RemoveField(Field *field);
    private:
        /// @brief Field, store all fields
        std::vector<Field *> m_field;
        /// @brief  field data communication information
        std::vector<FieldDataCommInfo *> m_field_data_comm_info;
    };
} // namespace zaran
