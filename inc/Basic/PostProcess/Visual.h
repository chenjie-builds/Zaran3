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
		void WriteTecASCII(shared_ptr<FieldManager> field_manager);
		void WriteTecplotASCII(shared_ptr<NSFieldStruct> field, std::ostream& os);
		void WriteTecplotASCII(shared_ptr<NSFieldZaran> field, std::ostream& os);
		void WriteTecplotASCII(shared_ptr<NSFieldFNFDM> field, std::ostream& os);

		void WriteTecplotBinary(shared_ptr<FieldManager> field_manager);
		void WriteTecplotBinary(shared_ptr<NSFieldFNFDM> field);
		void WriteTecplotBinary(shared_ptr<NSFieldZaran> field);
		void WriteTecplotBinary(shared_ptr<NSFieldStruct> field);

		void WriteVtkASCII(shared_ptr<FieldManager> field_manager);
		void WriteVtkASCII(shared_ptr<NSFieldStruct> field, std::ostream& os);
		void WriteVtkBinary(shared_ptr<FieldManager> field_manager);
		void WriteVtkBinary(shared_ptr<NSFieldStruct> field, std::ostream& os);
	};
}