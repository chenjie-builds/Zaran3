#include "FieldDataManager.h"
namespace zaran
{
    DataManager::DataManager( FieldData* fieldData, int equ_num) : m_data(fieldData), m_data_num(equ_num)
    {

    }
    DataManager::~DataManager()
    {
    }
}