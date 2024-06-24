//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Visual.h															||
//*	@brief	Á÷³¡ÏÔÊ¾																||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma	once
#include "Field.h"
#include "NSFieldStruct.h"
namespace zaran
{
	class Visual
	{
	public:
		void WriteTecplotBinary(Field* field);
		void WriteTecplotBinary(NSFieldStruct* field);
		void WriteVTK(Field* field);
	};
}