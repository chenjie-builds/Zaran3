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
#include <vector>
#include "grid.h"
namespace zaran
{
	//网格的集合，所有的网格均存在此处
	class GridList
	{
	public:
		GridList() {}
		GridList(std::vector<std::shared_ptr<Grid>>& gridList );
		~GridList();
		//向网格集合中添加一个网格
		void AddGrid(std::shared_ptr<Grid>& grid);
		//检查网格的index是否和网格集合中匹配
		void CheckGrid();
		//返回一个网格
		std::shared_ptr<Grid>& GetGrid(const int index);
		//返回网格的个数
		int GetGridNumber() { return gridList_.size(); }
		// 返回网格数组
		std::vector<std::shared_ptr<Grid>>& GetGrid();
		// 重设网格个数
		void Resize(int gridNum);
	private:
		std::vector<std::shared_ptr<Grid>> gridList_;
	};
}