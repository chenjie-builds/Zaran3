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
		static void Update(const string& dataName, const dataVariant& dataValue);
	private:
		GlobalData() {};
		~GlobalData() {};
		static const dataVariant& Get(const string& dataName);
		static map<string, dataVariant >dataMap_;
	};
}