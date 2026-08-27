/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file GlobalData.h
 * \brief Stores parameters for global use (TOML-backed).
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */

#pragma once
#include "BasicType.h"
#include <toml.hpp>

namespace zaran
{
	/// @brief 全局参数单例，后端为 toml::table，键名使用“分节”点分形式（如 task.dimension）
	class GlobalData
	{
	public:
		static GlobalData& Init();
		/// @brief 解析 TOML 控制文件
		static void Load(const std::string& file_name);
		static bool IsExist(const std::string& dataName);

		static int GetInt(const std::string& dataName);
		static double GetDouble(const std::string& dataName);
		static bool GetBool(const std::string& dataName);
		static std::string GetString(const std::string& dataName);

		static dynamic_array<double> GetDoubleArray(const std::string& dataName);
		static dynamic_array<int> GetIntArray(const std::string& dataName);
		static dynamic_array<std::string> GetStringArray(const std::string& dataName);

		static void Update(const std::string& dataName, int value);
		static void Update(const std::string& dataName, double value);
		static void Update(const std::string& dataName, bool value);
		static void Update(const std::string& dataName, const std::string& value);

		/// @brief 将当前全局参数序列化为 TOML 写入备份目录
		static void Backup(const std::string& back_folder);

	private:
		GlobalData() = default;
		~GlobalData() = default;
		toml::table m_table;
	};
}
