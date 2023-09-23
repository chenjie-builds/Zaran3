//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Grid_Struct_2D.h													||
//*	@brief	二维结构网格类，根据结构网格数据结构特点存储信息							||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "grid.h"
#include "NodeTopoInfoStruct.h"
namespace zaran
{

	class Grid_Struct_2D :public Grid
	{
	public:
		Grid_Struct_2D()
			:Grid(), m_ni(0), m_nj(0)
		{}
	public:
		void SetNi(int ni);
		void SetNj(int nj);
		int GetNi();
		int GetNj();
		void SetNodeNum(int ni, int nj);
		void GetNodeNum(int& ni, int& nj);
		int GetNodeIndex(int i, int j);
		void GetNodeIndex(int index, int& i, int& j);
		int GetCellIndex(int i, int j);
		void GetCellIndex(int index, int& i, int& j);
	public:
		void GetRange(int& iStart, int& iEnd, int& jStart, int& jEnd);
	private:
		int m_ni, m_nj;//三个方向的节点个数
	};
}