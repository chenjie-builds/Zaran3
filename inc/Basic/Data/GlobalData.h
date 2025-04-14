//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	GlobalData.h														||
//*	@brief	Store parameter for global use										||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include"BasicType.h"
#include <span>
namespace zaran
{
	using dataVariant = variant<string, int, double>;
	class GlobalData
	{
	public:
		static GlobalData& Init();
		static bool IsExist(const string& dataName);
		static const int GetInt(const string& dataName);
		static const double GetDouble(const string& dataName);
		static const string GetString(const string& dataName);
		template<typename T>
		static const T GetData(const string& dataName) { return *(std::get_if<T>(&Get(dataName))); }
		static void Update(const string& dataName, const dataVariant& dataValue);
		static void Backup(const string& back_file_name);
	private:
		GlobalData() {};
		~GlobalData() {};
		static const dataVariant& Get(const string& dataName);
		static map<string, dataVariant >dataMap_;
	};
}