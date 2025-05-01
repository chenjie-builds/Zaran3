/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file GlobalData.h
 * \brief Stores parameters for global use.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.	
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */

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