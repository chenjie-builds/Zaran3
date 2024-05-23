#pragma once
#include"FieldData.h"
#include "GridBase.h"
namespace zaran
{
    class DataManager
    {
    public:
        DataManager( FieldData* fieldData, int equ_num);
       virtual ~DataManager();
    public:
       virtual void CreateFieldData()=0;
       virtual void RegisterFieldData()=0;
    protected:
        FieldData* m_data;
        int m_equ_num;
    };
}