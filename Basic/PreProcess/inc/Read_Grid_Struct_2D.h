//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Read_Grid_Struct_2D.h															||
//*	@brief	读取二维结构网格网格											||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include"GridListFactory.h"
namespace zaran
{
	class GridFactoryStruct2D :public GridListFactory
	{
	public:
		GridFactoryStruct2D() {};
		void Create(Grid*& grid) override;
	};

}

