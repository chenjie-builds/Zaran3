#pragma once
#include"Field.h"
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
        Field* GetField(int index);
        /// @brief 获取Field数量
        int GetFieldSize();
        /// @brief 添加Field
        void AddField(Field* field);
        /// @brief 移除Field
        void RemoveField(Field* field);
    private:
        std::vector<Field*> m_field;

    };
} // namespace zaran
