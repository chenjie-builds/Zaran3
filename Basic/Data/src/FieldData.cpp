#include "FieldData.h"
#include<utility>
#include<iostream>
namespace zaran
{

	FieldData::FieldData()
	{

	}

	FieldData::~FieldData()
	{
		for (int i = 0; i < m_data.size(); i++)
		{
			if (m_data_type[i] == FieldDataType::integer)
			{
				int* data = (int*)m_data[i];
				delete[] data;
			}
			else if (m_data_type[i] == FieldDataType::real)
			{
				double* data = (double*)m_data[i];
				delete[] data;
			}
		}

	}

	void FieldData::EraseData(const string& name)
	{
		int id = m_name_id[name];
		m_name_id.erase(name);
		m_data_type.erase(m_data_type.begin() + id);
		m_data_num.erase(m_data_num.begin() + id);
		if (m_data_type[id] == FieldDataType::integer)
		{
			int* data = (int*)m_data[id];
			delete[] data;
		}
		else if (m_data_type[id] == FieldDataType::real)
		{
			double* data = (double*)m_data[id];
			delete[] data;
		}
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
#ifdef DEBUG
		if (HasData(name))
		{
			cout << "FieldData::AddData: " << name << " is already exist" << endl;
			exit(0);
		}
#endif // DEBUG
		m_name_id[name] = m_data_type.size();
		m_data_type.push_back(type);
		m_data_num.push_back(num);
		if (type == FieldDataType::integer)
		{
			m_data.push_back(new int[num]);
			for (int i = 0;i < num;i++)
			{
				((int*)m_data.back())[i] = 0;
			}
		}
		else if (type == FieldDataType::real)
		{
			m_data.push_back(new double[num]);
			for (int i = 0;i < num;i++)
			{
				((double*)m_data.back())[i] = 0.0;
			}
		}
	}

	void FieldData::GetData(const string& name, double*& data)
	{
		int id = m_name_id[name];
#ifdef DEBUG
		if (HasData(name) == false)
		{
			cout << "FieldData::GetData: " << name << " is not exist" << endl;
			exit(0);
		}
		if (m_data_type[id] != FieldDataType::real)
		{
			cout << "FieldData::GetData: " << name << " is not real type" << endl;
			exit(0);
		}
#endif // DEBUG
		data = (double*)m_data[id];
	}
	void FieldData::GetData(const string& name, int*& data)
	{
		int id = m_name_id[name];
#ifdef DEBUG
		if (HasData(name) == false)
		{
			cout << "FieldData::GetData: " << name << " is not exist" << endl;
			exit(0);
		}
		if (m_data_type[id] != FieldDataType::integer)
		{
			cout << "FieldData::GetData: " << name << " is not integer type" << endl;
			exit(0);
		}
#endif // DEBUG
		data = (int*)m_data[id];
	}
	void FieldData::GetDataSize(const string& name, int& num)
	{
		int id = m_name_id[name];
#ifdef DEBUG
		if (HasData(name) == false)
		{
			cout << "FieldData::GetData: " << name << " is not exist" << endl;
			exit(0);
		}
#endif // DEBUG
		num = m_data_num[id];
	}
}