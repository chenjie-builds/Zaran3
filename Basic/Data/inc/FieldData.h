//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	FieldData.h															||
//*	@brief	场数据类																||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include"BasicType.h"
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
		void EraseData(const string& name);
		bool HasData(const string& name)const;
		void AddData(const string& name, FieldDataType type, int num);
		void GetData(const string& name, double*& data);
		void GetData(const string& name, int*& data);
		void GetDataSize(const string& name, int& num);
	private:
		/// @brief 用于记录数据名称和数据的对应关系
		map<string, int>m_name_id;
		/// @brief 用于记录数据类型
		std::vector<FieldDataType> m_data_type;
		/// @brief 用于记录数据个数
		std::vector<int> m_data_num;
		/// @brief 用于记录数据
		/// @note 第一维是数据索引，第二维是数据
		std::vector<void*> m_data;
	};
}

