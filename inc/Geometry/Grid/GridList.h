//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	GridList.h															||
//*	@brief	网格数组，存储所有的计算网格											||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "grid.h"
namespace zaran
{
	//网格的集合，所有的网格均存在此处
	class GridList
	{
	public:
		GridList() {}
		GridList(Array<Ptr<Grid>>& gridList);
		~GridList();
		//向网格集合中添加一个网格
		void AddGrid(Ptr<Grid>& grid);
		//检查网格的index是否和网格集合中匹配
		void CheckGrid();
		//返回一个网格
		Ptr<Grid>& GetGrid(const int index);
		//返回网格的个数
		int GetGridNumber() { return gridList_.size(); }
		// 返回网格数组
		Array<Ptr<Grid>>& GetGrid();
		// 重设网格个数
		void Resize(int gridNum);
	private:
		Array<Ptr<Grid>> gridList_;
	};
}