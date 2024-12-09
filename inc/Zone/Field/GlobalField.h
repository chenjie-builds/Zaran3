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
        std::shared_ptr<Field> GetField(Id idx_field);
        /// @brief 获取field data communication information
        std::shared_ptr<FieldDataCommInfo> GetFieldDataCommInfo(Id idx_field);
        void SetFieldDataCommInfo(Id idx_field, std::shared_ptr<FieldDataCommInfo> field_data_comm_info);
        /// @brief 获取Field数量
        Id GetFieldNum();
        /// @brief 添加Field
        void AddField(std::shared_ptr<Field> field, std::shared_ptr<FieldDataCommInfo> field_data_comm_info);
        /// @brief 移除Field
        void RemoveField(std::shared_ptr<Field> field);
    private:
        /// @brief Field, store all fields
        Array<std::shared_ptr<Field>> m_field;
        /// @brief  field data communication 
        Array<std::shared_ptr<FieldDataCommInfo>> m_field_data_comm_info;
    };
} // namespace zaran
