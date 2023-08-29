//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	GridListFactory.h													||
//*	@brief	网格数组工厂, 负责生成网格数组											||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "GridList.h"
#include "GlobalData.h"
namespace zaran
{
	class GridListFactory
	{
	public:
		virtual void Create(Ptr<GridList>& gridList);
	private:
		void ReadPlot3D(Ptr<GridList>& gridList);
		void CreateByTest(Ptr<GridList>& gridList);
	};
}