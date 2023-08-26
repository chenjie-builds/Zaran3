#include "GlobalData.h"
#include "log.h"

GlobalData& GlobalData::Init()
{
	static GlobalData data;
	return data;
}

bool GlobalData::IsExist(const std::string& dataName)
{
	if (dataMap_.find(dataName) != dataMap_.end())
	{
		return true;
	}
	return false;
}

const dataVariant& GlobalData::Get(const std::string& dataName)
{
	auto data = dataMap_.find(dataName);
	if (data != dataMap_.end())
		return data->second;
	else
	{
		ZaranLog::warn("Not Found:{},in Global Data!", dataName);
		ZaranLog::warn("Exit Now!");
		exit(0);
	}
}

const int GlobalData::GetInt(const std::string& dataName)
{
	return *(std::get_if<int>(&Get(dataName)));
}

const double GlobalData::GetDouble(const std::string& dataName)
{
	return *(std::get_if<double>(&Get(dataName)));
}

const std::string GlobalData::GetString(const std::string& dataName)
{
	return *(std::get_if<std::string>(&Get(dataName)));
}

std::unordered_map<std::string, dataVariant > GlobalData::dataMap_;

void GlobalData::Update(const std::string& varName, const dataVariant& varValue)
{
	dataMap_[varName] = varValue;
}
