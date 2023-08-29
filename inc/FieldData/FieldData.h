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
#include"BasicType.h"
namespace zaran
{
	class FieldData
	{
	public:
		FieldData();
		~FieldData();
		void AddData(const string& name, DArray& data);
		void EraseData(const string& name);
		bool HasData(const string& name)const;
		DArray& GetData(const string& name);
	private:
		map<string, DArray> m_data;
	};
}

