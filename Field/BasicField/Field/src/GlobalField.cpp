#include "GlobalField.h"
namespace zaran
{
    void FieldManager::RemoveField(Field* field)
    {
        for (auto it = m_field.begin(); it != m_field.end(); ++it)
        {
            if (*it == field)
            {
                delete *it;
                m_field.erase(it);
                break;
            }
        }
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
    }
    Field* FieldManager::GetField(int index)
    {
        return m_field[index];
    }
    int FieldManager::GetFieldSize()
    {
        return m_field.size();
    }
    void FieldManager::AddField(Field* field)
    {
        field->SetIdx(m_field.size());
        m_field.push_back(field);
    }
} // namespace zaran