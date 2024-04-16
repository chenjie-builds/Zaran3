#include "GlobalData.h"
#include "log.h"
using namespace zaran;
GlobalData& GlobalData::Init()
{
	static GlobalData data;
	return data;
}

bool GlobalData::IsExist(const string& dataName)
{
	if (dataMap_.find(dataName) != dataMap_.end())
	{
		return true;
	}
	return false;
}

const dataVariant& GlobalData::Get(const string& dataName)
{
	auto data = dataMap_.find(dataName);
	if (data != dataMap_.end())
		return data->second;
	else
	{
		Log::warn("Not Found:{},in Global Data!", dataName);
		Log::warn("Exit Now!");
		exit(0);
	}
}

const int GlobalData::GetInt(const string& dataName)
{
	return *(std::get_if<int>(&Get(dataName)));
}

const double GlobalData::GetDouble(const string& dataName)
{
	return *(std::get_if<double>(&Get(dataName)));
}

const string GlobalData::GetString(const string& dataName)
{
	return *(std::get_if<string>(&Get(dataName)));
}

map<string, dataVariant > GlobalData::dataMap_;

void GlobalData::Update(const string& varName, const dataVariant& varValue)
{
	dataMap_[varName] = varValue;
}

void zaran::GlobalData::Backup(const string& back_file_name)
{

}
