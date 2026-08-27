#include "GlobalData.h"
#include "Log.h"
#include "ZaranError.h"
#include <fstream>
#include <sstream>

namespace zaran
{

namespace
{
	std::vector<std::string> SplitDotted(const std::string& key)
	{
		std::vector<std::string> parts;
		std::stringstream ss(key);
		std::string part;
		while (std::getline(ss, part, '.'))
			parts.push_back(part);
		return parts;
	}

	// 沿点分键逐级查找节点；找不到返回 nullptr
	toml::node* FindNode(toml::table& root, const std::string& dotted)
	{
		auto parts = SplitDotted(dotted);
		toml::node* cur = &root;
		for (const auto& part : parts)
		{
			toml::table* tbl = cur->as_table();
			if (!tbl)
				return nullptr;
			auto it = tbl->find(part);
			if (it == tbl->end())
				return nullptr;
			cur = &(it->second);
		}
		return cur;
	}

	// 沿点分键逐级确保子表存在（必要时创建），返回前 upto 段路径对应的表
	toml::table* EnsureTable(toml::table& root, const std::vector<std::string>& parts, size_t upto)
	{
		toml::table* cur = &root;
		for (size_t i = 0; i < upto; ++i)
		{
			auto it = cur->find(parts[i]);
			if (it == cur->end())
			{
				cur->insert_or_assign(parts[i], toml::table{});
				it = cur->find(parts[i]);
			}
			cur = it->second.as_table();
			if (!cur)
				return nullptr;
		}
		return cur;
	}
}

GlobalData& GlobalData::Init()
{
	static GlobalData data;
	return data;
}

void GlobalData::Load(const std::string& file_name)
{
	try
	{
		// toml++ 默认开启异常：parse_file 成功返回 toml::table，失败抛 toml::parse_error
		Init().m_table = toml::parse_file(file_name);
	}
	catch (const toml::parse_error& err)
	{
		std::stringstream ss;
		ss << err;
		throw ZaranError("Failed to parse control file: " + file_name + " : " + ss.str());
	}
}

bool GlobalData::IsExist(const std::string& dataName)
{
	return FindNode(Init().m_table, dataName) != nullptr;
}

int GlobalData::GetInt(const std::string& dataName)
{
	toml::node* node = FindNode(Init().m_table, dataName);
	if (!node)
		throw ZaranError("GlobalData: key not found: " + dataName);
	auto v = node->value<int64_t>();
	if (!v)
		throw ZaranError("GlobalData: key is not an integer: " + dataName);
	return static_cast<int>(*v);
}

double GlobalData::GetDouble(const std::string& dataName)
{
	toml::node* node = FindNode(Init().m_table, dataName);
	if (!node)
		throw ZaranError("GlobalData: key not found: " + dataName);
	auto v = node->value<double>();
	if (!v)
		throw ZaranError("GlobalData: key is not a double: " + dataName);
	return *v;
}

bool GlobalData::GetBool(const std::string& dataName)
{
	toml::node* node = FindNode(Init().m_table, dataName);
	if (!node)
		throw ZaranError("GlobalData: key not found: " + dataName);
	auto v = node->value<bool>();
	if (!v)
		throw ZaranError("GlobalData: key is not a boolean: " + dataName);
	return *v;
}

std::string GlobalData::GetString(const std::string& dataName)
{
	toml::node* node = FindNode(Init().m_table, dataName);
	if (!node)
		throw ZaranError("GlobalData: key not found: " + dataName);
	auto v = node->value<std::string>();
	if (!v)
		throw ZaranError("GlobalData: key is not a string: " + dataName);
	return *v;
}

dynamic_array<double> GlobalData::GetDoubleArray(const std::string& dataName)
{
	toml::node* node = FindNode(Init().m_table, dataName);
	if (!node || !node->is_array())
		throw ZaranError("GlobalData: double[] key not found or not an array: " + dataName);
	dynamic_array<double> result;
	const auto& arr = *node->as_array();
	result.reserve(arr.size());
	for (const auto& elem : arr)
	{
		auto v = elem.value<double>();
		if (!v)
			throw ZaranError("GlobalData: array element is not a double: " + dataName);
		result.push_back(*v);
	}
	return result;
}

dynamic_array<int> GlobalData::GetIntArray(const std::string& dataName)
{
	toml::node* node = FindNode(Init().m_table, dataName);
	if (!node || !node->is_array())
		throw ZaranError("GlobalData: int[] key not found or not an array: " + dataName);
	dynamic_array<int> result;
	const auto& arr = *node->as_array();
	result.reserve(arr.size());
	for (const auto& elem : arr)
	{
		auto v = elem.value<int64_t>();
		if (!v)
			throw ZaranError("GlobalData: array element is not an integer: " + dataName);
		result.push_back(static_cast<int>(*v));
	}
	return result;
}

dynamic_array<std::string> GlobalData::GetStringArray(const std::string& dataName)
{
	toml::node* node = FindNode(Init().m_table, dataName);
	if (!node || !node->is_array())
		throw ZaranError("GlobalData: string[] key not found or not an array: " + dataName);
	dynamic_array<std::string> result;
	const auto& arr = *node->as_array();
	result.reserve(arr.size());
	for (const auto& elem : arr)
	{
		auto v = elem.value<std::string>();
		if (!v)
			throw ZaranError("GlobalData: array element is not a string: " + dataName);
		result.push_back(*v);
	}
	return result;
}

void GlobalData::Update(const std::string& dataName, int value)
{
	auto parts = SplitDotted(dataName);
	toml::table* tbl = EnsureTable(Init().m_table, parts, parts.size() - 1);
	if (!tbl)
		throw ZaranError("GlobalData: cannot update key under non-table: " + dataName);
	tbl->insert_or_assign(parts.back(), value);
}

void GlobalData::Update(const std::string& dataName, double value)
{
	auto parts = SplitDotted(dataName);
	toml::table* tbl = EnsureTable(Init().m_table, parts, parts.size() - 1);
	if (!tbl)
		throw ZaranError("GlobalData: cannot update key under non-table: " + dataName);
	tbl->insert_or_assign(parts.back(), value);
}

void GlobalData::Update(const std::string& dataName, bool value)
{
	auto parts = SplitDotted(dataName);
	toml::table* tbl = EnsureTable(Init().m_table, parts, parts.size() - 1);
	if (!tbl)
		throw ZaranError("GlobalData: cannot update key under non-table: " + dataName);
	tbl->insert_or_assign(parts.back(), value);
}

void GlobalData::Update(const std::string& dataName, const std::string& value)
{
	auto parts = SplitDotted(dataName);
	toml::table* tbl = EnsureTable(Init().m_table, parts, parts.size() - 1);
	if (!tbl)
		throw ZaranError("GlobalData: cannot update key under non-table: " + dataName);
	tbl->insert_or_assign(parts.back(), value);
}

void GlobalData::Backup(const std::string& back_folder)
{
	std::ofstream fout(back_folder + "/zaran.toml");
	if (!fout.is_open())
		throw ZaranError("GlobalData: cannot write control file to: " + back_folder);
	toml::toml_formatter formatter{ Init().m_table };
	fout << formatter;
}

} // namespace zaran
