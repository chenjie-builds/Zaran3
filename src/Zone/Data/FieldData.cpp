#include "FieldData.h"
#include <utility>
#include <iostream>
namespace zaran
{

	FieldData::FieldData()
	{
	}

	FieldData::~FieldData()
	{

	}

	void FieldData::EraseData(const string& name)
	{
		int id = m_name_id[name];
		m_name_id.erase(name);
		m_data_type.erase(m_data_type.begin() + id);
		m_data_num.erase(m_data_num.begin() + id);
		m_data.erase(m_data.begin() + id);
		for (auto& it : m_name_id)
		{
			if (it.second > id)
			{
				it.second--;
			}
		}
	}

	bool FieldData::HasData(const string& name) const
	{
		return m_name_id.find(name) != m_name_id.end();
	}

	void FieldData::AddData(const string& name, FieldDataType type, int num)
	{
#ifdef DEBUG_MODE
		if (HasData(name))
		{
			std::cout << "FieldData::AddData: " << name << " is already exist" << std::endl;
			exit(0);
		}
#endif // DEBUG_MODE
		m_name_id[name] = m_data_type.size();
		m_data_type.push_back(type);
		m_data_num.push_back(num);
		if (type == FieldDataType::integer)
		{
			m_data.push_back(dynamic_array<std::byte>(sizeof(int) * num));
			int* data = reinterpret_cast<int*>(m_data.back().data());
			for (int i = 0; i < num; i++)
			{
				data[i] = 0;
			}
		}
		else if (type == FieldDataType::real)
		{
			m_data.push_back(dynamic_array<std::byte>(sizeof(double) * num));
			double* data = reinterpret_cast<double*>(m_data.back().data());
			for (int i = 0; i < num; i++)
			{
				data[i] = 0.0;
			}
		}
	}

	void FieldData::GetData(const string& name, double*& data)
	{
		int id = m_name_id[name];
#ifdef DEBUG_MODE
		if (HasData(name) == false)
		{
			std::cout << "FieldData::GetData: " << name << " is not exist" << std::endl;
			exit(0);
		}
		if (m_data_type[id] != FieldDataType::real)
		{
			std::cout << "FieldData::GetData: " << name << " is not real type" << std::endl;
			exit(0);
		}
#endif // DEBUG_MODE
		data = reinterpret_cast<double*>(m_data[id].data());
	}
	void FieldData::GetData(const string& name, int*& data)
	{
		int id = m_name_id[name];
#ifdef DEBUG_MODE
		if (HasData(name) == false)
		{
			std::cout << "FieldData::GetData: " << name << " is not exist" << std::endl;
			exit(0);
		}
		if (m_data_type[id] != FieldDataType::integer)
		{
			std::cout << "FieldData::GetData: " << name << " is not integer type" << std::endl;
			exit(0);
		}
#endif // DEBUG_MODE
		data = reinterpret_cast<int*>(m_data[id].data());
	}
	void FieldData::GetDataSize(const string& name, int& num)
	{
		int id = m_name_id[name];
#ifdef DEBUG_MODE
		if (HasData(name) == false)
		{
			std::cout << "FieldData::GetData: " << name << " is not exist" << std::endl;
			exit(0);
		}
#endif // DEBUG_MODE
		num = m_data_num[id];
	}
	double& FieldData::GetData(const string& name, int idx_data)
	{
		int id = m_name_id[name];
		return reinterpret_cast<double*>(m_data[id].data())[idx_data];
	}
}