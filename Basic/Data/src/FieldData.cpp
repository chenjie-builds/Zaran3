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

	}

	void FieldData::AddData(const std::string& name, DArray& data)
	{
		m_data.insert(std::make_pair(name, data));
	}
	void FieldData::EraseData(const std::string& name)
	{
		m_data.erase(name);
	}
	bool FieldData::HasData(const std::string& name) const
	{
		return m_data.find(name) != m_data.end();
	}

	DArray& FieldData::GetData(const std::string& name)
	{
		return m_data.at(name);
	}
}