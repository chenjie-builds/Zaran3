//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	FNFDM2D.h															||
//*	@brief	读取二维维自由节点有限差分网格											||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include"GridListFactory.h"
namespace zaran
{
	class GridListFactoryFNFDM2D :public GridListFactory
	{
	public:
		GridListFactoryFNFDM2D() {};
		void Create(Ptr<GridList>& gridList) override;
	};

}

