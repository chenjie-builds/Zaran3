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
#include <unordered_map>
#include <string>
#include <variant>
using dataVariant = std::variant<std::string, int, double>;
class GlobalData
{
public:
	static GlobalData& Init();
	static bool IsExist(const std::string& dataName);
	static const int GetInt(const std::string& dataName);
	static const double GetDouble(const std::string& dataName);
	static const std::string GetString(const std::string& dataName);
	static void Update(const std::string& dataName, const dataVariant& dataValue);
private:
	GlobalData() {};
	~GlobalData() {};
	static const dataVariant& Get(const std::string& dataName);
	static std::unordered_map<std::string, dataVariant >dataMap_;
};

