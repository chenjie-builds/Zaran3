#include "Field.h"
namespace zaran
{
    Field::Field(GridBase* grid, FieldType fieldType) :
        m_grid(grid), m_fieldType(fieldType)
    {
        m_fieldData = nullptr;
        m_solver = nullptr;
        m_solver_para = nullptr;
        m_data_manager = nullptr;
        m_fieldData = new FieldData();
    }
    Field::~Field()
    {
        if (m_fieldData != nullptr)
        {
            delete m_fieldData;
        }
        if (m_solver != nullptr)
        {
            delete m_solver;
        }
        if (m_solver_para != nullptr)
        {
            delete m_solver_para;
        }
        if (m_data_manager != nullptr)
        {
            delete m_data_manager;
        }
        if (m_grid != nullptr)
        {
            delete m_grid;
        }
    }
    void Field::Allocate()
    {
        if(m_fieldData != nullptr)
        {
            delete m_fieldData;
        }
        m_fieldData = new FieldData();
    }
}