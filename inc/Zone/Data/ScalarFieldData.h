#pragma once
#include "FieldDataBase.h"
#include "BasicType.h"
template<typename T>
class ScalarField final : public FieldBase 
{
public:
	ScalarField(const std::string& name, int data_num, T initial_value = 0.0)
		: FieldBase(name), m_values(data_num, initial_value) 
	{
	}

	T GetValue(int node_id) const 
	{
		if (node_id < 0 || node_id >= static_cast<int>(m_values.size()))
			throw std::out_of_range("Node ID out of range in ScalarField::GetValue");
		return m_values[node_id];
	}

	void SetValue(int node_id, T value) 
	{
		if (node_id < 0 || node_id >= static_cast<int>(m_values.size()))
			throw std::out_of_range("Node ID out of range in ScalarField::SetValue");
		m_values[node_id] = value;
	}

	std::string Info() const override 
	{
		std::ostringstream oss;
		oss << "ScalarField(name=" << m_name
			<< ", node_count=" << m_values.size() << ")";
		return oss.str();
	}
	count_type GetDataCount() const { return m_values.size(); }

	std::string GetType()
	{
		if constexpr (std::is_same_v<T, float>)
			return "float";
		else if constexpr (std::is_same_v<T, double>)
			return "double";
		else if constexpr (std::is_same_v<T, int>)
			return "int";
		else if constexpr (std::is_same_v<T, unsigned int>)
			return "unsigned int";
		else
			return "unknown";
	}
private:
	//check if the index is valid	
	void CheckIndex(index_type index) const
	{
		if (index < 0 || index >= GetDataCount())
			throw std::out_of_range("Index out of range in VectorField::CheckIndex");
	}
private:
	dynamic_array<T> m_values; // Values of the scalar field
};