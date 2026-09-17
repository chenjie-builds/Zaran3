/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file DataManagerNSTwoPhase.cpp
 * \brief Data manager for the two-phase (volume-fraction) NS/Euler field.
 */

#include "DataManagerNSTwoPhase.h"

namespace zaran
{
	DataManagerNSTwoPhase::DataManagerNSTwoPhase(shared_ptr<FieldData> field_data, int ni, int nj, int nk)
		: DataManagerNSStruct(field_data, ni, nj, nk)
	{
	}

	DataManagerNSTwoPhase::~DataManagerNSTwoPhase()
	{
		// m_volume_fraction 指向 FieldData 内部缓冲，不负责释放
	}

	void DataManagerNSTwoPhase::CreateData()
	{
		DataManagerNSStruct::CreateData();
		m_data->AddData("volume_fraction", FieldDataType::real, m_data_num);
	}

	void DataManagerNSTwoPhase::RegisterData()
	{
		DataManagerNSStruct::RegisterData();
		m_data->GetData("volume_fraction", m_volume_fraction);
	}
}
