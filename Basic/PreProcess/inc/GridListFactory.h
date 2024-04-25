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
		virtual void Create(Grid*& grid);
	private:
		void ReadPlot3D(Grid*& grid);
		void CreateByTest(Grid*& grid);
		// 生成测试的三维结构网格
		void CreateStructGrid3D(Grid*& grid);
		// 生成测试的二维结构网格
		void CreateStructGrid2D(Grid*& grid);
		// 生成测试的二维自由节点有限差分网格
		void CreateGridFNFDM2D(Grid*& grid);
	};
}