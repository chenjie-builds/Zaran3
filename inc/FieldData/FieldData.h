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
#include<vector>
#include<string>
#include <memory>
#include<unordered_map>
namespace zaran
{
	class FieldData
	{
	public:
		FieldData();
		~FieldData();
		void AddData(const std::string& name,  std::vector<double>& data);
		void AddData(const std::string& name, std::vector<double>&& data);
		void EraseData(const std::string& name);
		bool HasData(const std::string& name)const;
		std::vector<double>& GetData(const std::string& name);
	private:
		std::unordered_map<std::string, std::vector<double>> m_data;
	};
}

