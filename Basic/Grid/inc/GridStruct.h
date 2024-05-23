//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Grid_Struct_3D.h													||
//*	@brief	三维结构网格类，根据结构网格数据结构特点存储信息						||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "GridBase.h"
#include "NodeStruct.h"
#include "FaceStruct.h"
#include "CellStruct.h"
namespace zaran
{
	class GridStruct :public GridBase
	{
	public:
		GridStruct(const string& name, int index, int dim, GridType type);
		virtual ~GridStruct();
	public:
		int GetNi();
		int GetNj();
		int GetNk();
		void SetNodeNum(int ni, int nj, int nk);
		void GetNodeNum(int& ni, int& nj, int& nk);
	public:
		void GetRange(int& iStart, int& iEnd, int& jStart, int& jEnd, int& kStart, int& kEnd);
	private:
		int m_ni, m_nj, m_nk;//三个方向的节点个数
		NodeStruct* m_node;
		FaceStruct* m_face;
		CellStruct* m_cell;
	};
}