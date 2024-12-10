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
        /// @brief 获取Field
        shared_ptr<Field> GetField(index_type idx_field);
        /// @brief 获取field data communication information
        shared_ptr<FieldDataCommInfo> GetFieldDataCommInfo(index_type idx_field);
        void SetFieldDataCommInfo(index_type idx_field, shared_ptr<FieldDataCommInfo> field_data_comm_info);
        /// @brief 获取Field数量
        index_type GetFieldNum();
        /// @brief 添加Field
        void AddField(shared_ptr<Field> field, shared_ptr<FieldDataCommInfo> field_data_comm_info);
        /// @brief 移除Field
        void RemoveField(shared_ptr<Field> field);
    private:
        /// @brief Field, store all fields
        dynamic_array<shared_ptr<Field>> m_field;
        /// @brief  field data communication 
        dynamic_array<shared_ptr<FieldDataCommInfo>> m_field_data_comm_info;
    };
} // namespace zaran
