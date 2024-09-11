#include "GlobalField.h"
namespace zaran
{
    void FieldManager::RemoveField(Field *field)
    {
        m_field.erase(m_field.begin() + field->GetIdx());
        m_field_data_comm_info.erase(m_field_data_comm_info.begin() + field->GetIdx());
    }
    FieldManager::FieldManager()
    {
        
    }
    FieldManager::~FieldManager()
    {
        FreeMemory();
    }

    void FieldManager::FreeMemory()
    {
        for (auto field : m_field)
        {
            delete field;
        }
        m_field.clear();
        for (auto field_data_comm_info : m_field_data_comm_info)
        {
            delete field_data_comm_info;
        }
        m_field_data_comm_info.clear();
    }
    Field *FieldManager::GetField(int idx_field)
    {
        return m_field[idx_field];
    }
    FieldDataCommInfo *FieldManager::GetFieldDataCommInfo(int idx_field)
    {
        return m_field_data_comm_info[idx_field];
    }
    int FieldManager::GetFieldNum()
    {
        return m_field.size();
    }
    void FieldManager::AddField(Field *field, FieldDataCommInfo *field_data_comm_info)
    {
        field->SetIdx(m_field.size());
        m_field.push_back(field);
        m_field_data_comm_info.push_back(field_data_comm_info);
    }
} // namespace zaran