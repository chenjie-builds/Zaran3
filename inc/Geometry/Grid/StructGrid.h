//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	StructGrid.h														||
//*	@brief	结构网格类，根据结构网格数据结构特点存储信息							||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "grid.h"
#include "NodeTopoInfoStruct.h"
namespace zaran
{
	class StructGrid :public Grid
	{
	public:
		StructGrid()
			:Grid(), ni_(0), nj_(0), nk_(0)
		{}
	public:
		void SetNi(int ni);
		void SetNj(int nj);
		void SetNk(int nk);
		int GetNi();
		int GetNj();
		int GetNk();
		void SetNodeNum(int ni, int nj, int nk);
		void GetNodeNum(int& ni, int& nj, int& nk);
		int GetNodeIndex(int i, int j, int k);
		void GetNodeIndex(int index, int& i, int& j, int& k);
		int GetCellIndex(int i, int j, int k);
		void GetCellIndex(int index, int& i, int& j, int& k);
	public:
		void GetRange(int& iStart, int& iEnd, int& jStart, int& jEnd, int& kStart, int& kEnd);
	private:
		int ni_, nj_, nk_;//三个方向的节点个数
	};
}