/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file VectorField.h
 * \brief VectorField class, used to store vector field data.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include <span.hpp>
#include "FieldDataBase.h"
#include "BasicType.h"
template<typename T>
class VectorField final : public FieldBase
{
public:
	VectorField(const std::string& name, count_type data_num,count_type data_dim, T initial_value = 0.0)
		: FieldBase(name), m_data_dim(data_dim), m_values(data_num* data_dim, initial_value)
	{
	}
	T GetValue(index_type idx, index_type dim) const
	{
		CheckIndex(idx);
		CheckDimension(dim);
		return m_values[idx * m_data_dim + dim];
	}
	void SetValue(index_type idx, index_type dim, T value)
	{
		CheckIndex(idx);
		CheckDimension(dim);
		m_values[idx * m_data_dim + dim] = value;
	}
	std::string Info() const override
	{
		std::ostringstream oss;
		oss << "VectorField(name=" << m_name
			<< ", data count=" << m_data_dim << ", data dimension=" << m_data_dim << ")";
		return oss.str();
	}
	span<T> GetValue(index_type idx) const
	{
		CheckIndex(idx);
		return span<T>(m_values.data() + idx * m_data_dim, m_data_dim);
	}
	T* GetData() { return m_values.data(); }
	const T* GetData() const { return m_values.data(); }
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
	//check if the dimension is valid
	void CheckDimension(index_type dim) const
	{
		if (dim < 0 || dim >= m_data_dim)
			throw std::out_of_range("Dimension out of range in VectorField::CheckDimension");
	}
private:
	dynamic_array<T> m_values;
	count_type m_data_dim;// dimension of data
};