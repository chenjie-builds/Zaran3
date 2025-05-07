/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file FieldData.h
 * \brief FieldData class, used to store field data.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "BasicType.h"
#include <cstring>
namespace zaran
{
	enum class FieldDataType
	{
		integer,
		real,
	};
	class FieldData
	{
	public:
		FieldData();
		~FieldData();
		void EraseData(const string &name);
		bool HasData(const string &name) const;
		void AddData(const string &name, FieldDataType type, int num);
		void GetData(const string &name, double *&data);
		void GetData(const string &name, int *&data);
		double& GetData(const std::string &name, int index);
		void GetDataSize(const string &name, int &num);
		template<typename T>
		T* GetData(const string& name)
		{
			int id = m_name_id[name];
			if (id < 0 || id >= m_data_num[id])
			{
				throw std::runtime_error("FieldData::GetData: invalid data index");
			}
			return reinterpret_cast<T*>(m_data[id].data());
		}
	private:
		/// @brief 用于记录数据名称和数据的对应关系
		map<string, int> m_name_id;
		/// @brief 用于记录数据类型
		dynamic_array<FieldDataType> m_data_type;
		/// @brief 用于记录数据个数
		dynamic_array<int> m_data_num;
		/// @brief 用于记录数据
		/// @note 第一维是数据索引，第二维是数据
		dynamic_array<dynamic_array<std::byte>> m_data;
	};



}
