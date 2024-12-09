#include "Field.h"
namespace zaran
{
    Field::Field(std::shared_ptr<GridBase> grid, FieldType fieldType) :
        m_grid(grid), m_fieldType(fieldType)
    {
    }
    Field::~Field()
    {
    }
    void Field::Allocate()
    {
        AllocateFieldData();
        AllocateSolverPara();
        AllocateDataManager();
        AllocateSolver();

    }
    void Field::AllocateFieldData()
    {
		m_data = std::make_shared<FieldData>();
    }
}