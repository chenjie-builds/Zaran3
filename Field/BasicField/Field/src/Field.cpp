#include "Field.h"
namespace zaran
{
    Field::Field()
    {
        m_fieldData = std::make_shared<FieldData>();
    }

    Field::~Field()
    {
    }


}