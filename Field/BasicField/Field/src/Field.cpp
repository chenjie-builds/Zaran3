#include "Field.h"
namespace zaran
{
    Field::Field()
    {
        m_fieldData =new FieldData();
    }

    Field::~Field()
    {
        std::cout<<"FieldÎö¹¹º¯Êý"<<std::endl;
        delete m_fieldData;
        delete m_solver;
        delete m_solverPara;
        delete m_grid;
    }


}