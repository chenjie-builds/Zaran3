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
#include"BoundMapStruct.h"
namespace zaran
{
	class GridStruct :public GridBase
	{
	public:
		GridStruct(const string& name, int index, int dim);
		virtual ~GridStruct();
		void Allocate(int ni, int nj, int nk, int ghost_size);
	public:
		int GetNi();
		int GetNj();
		int GetNk();
		void GetNodeNum(int& ni, int& nj, int& nk);
		int GetTotalNodeNum();
		NodeStruct* GetNode() { return m_node; }
		FaceStruct* GetFace() { return m_face; }
		CellStruct* GetCell() { return m_cell; }
		BoundMapStruct* GetBoundMap() { return m_bound_map; }
	public:
		void GetRange(int& iStart, int& iEnd, int& jStart, int& jEnd, int& kStart, int& kEnd);
	private:
		int m_ghost_size;
		NodeStruct* m_node;
		FaceStruct* m_face;
		CellStruct* m_cell;
		BoundMapStruct* m_bound_map;
	};
}