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
#pragma once
#include "Field.h"
#include "NSFieldStruct.h"
#include "NSFieldZaran.h"
namespace zaran
{
	class Visual
	{
	public:
		void WriteTecASCII(FieldManager* field_manager);
		void WriteTecplotASCII(NSFieldStruct* field, std::ostream& os);
		void WriteTecplotASCII(NSFieldZaran* field, std::ostream& os);

		void WriteTecplotBinary(FieldManager* field_manager);
		void WriteTecplotBinary(NSFieldFNFDM* field);
		void WriteTecplotBinary(NSFieldZaran* field);
		void WriteTecplotBinary(NSFieldStruct* field);
	};
}