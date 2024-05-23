#include "Field.h"
namespace zaran
{
    Field::Field(GridBase* grid, FieldType fieldType) :
        m_grid(grid), m_fieldType(fieldType)
    {
    }
    Field::~Field()
    {
        delete m_fieldData;
        delete m_solver;
        delete m_solver_para;
        delete m_grid;
        delete m_dataManager;
    }
    void Field::Allocate()
    {
        m_fieldData = new FieldData();

    }
}