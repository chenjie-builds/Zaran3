#include "FieldDataManager.h"
namespace zaran
{
    DataManager::DataManager(std::shared_ptr<FieldData> fieldData, int data_num) : m_data(fieldData), m_data_num(data_num)
    {

    }
    DataManager::~DataManager()
    {
    }
}