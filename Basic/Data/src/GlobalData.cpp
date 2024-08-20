#include "GlobalData.h"
#include "Log.h"
#include <fstream>
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

void zaran::GlobalData::Backup(const string& back_folder)
{
	std::string data_file = "zaran.ini";
	std::ifstream fin(data_file);
	std::ofstream fout(back_folder + "/" + data_file);
	std::string line;
	std::string dataType;
	std::string dataName;
	std::string dataValue;
	while (std::getline(fin, line))
	{
		std::string separator = " =\t#$,;\"";
		if (line.empty())
		{
			fout << line << std::endl;
			continue;
		}
		if (line[0] == '!' || line[0] == '#' || line[0] == '/')
		{
			fout << line << std::endl;
			continue;
		}
		size_t id = line.find_first_not_of(' ');
		line.erase(line.begin(), line.begin() + id);
		id = line.find_first_of(separator);
		dataType = line.substr(0, id);
		line.erase(0, id);
		line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
		if (line.empty())
			continue;
		id = line.find_first_of("=");
		if (id == std::string::npos)
			continue;
		dataName = line.substr(0, id);
		dataValue = line.substr(id + 1);
		if (dataType == "string")
			fout << dataType << "\t\t" << dataName << "\t\t=\t\t" << GetString(dataName) << std::endl;
		else if (dataType == "double")
			fout << dataType << "\t\t" << dataName << "\t\t=\t\t" << GetDouble(dataName) << std::endl;
		else if (dataType == "int")
			fout << dataType << "\t\t" << dataName << "\t\t=\t\t" << GetInt(dataName) << std::endl;
		else
			Log::warn("Unsupported Data Type:{}, Name:{}, Value:{}", dataType, dataName, dataValue);
	}
}
