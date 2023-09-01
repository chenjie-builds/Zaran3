//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	SYSU.h																||
//*	@brief	中山大学网格生成														||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include"GridListFactory.h"
namespace zaran
{
	const int SYSU_INLET = 1;
	const int SYSU_OUTLET = 2;
	const int SYSU_WALL = 3;

	class GridListFactorySYSU :public GridListFactory
	{
	public:
		GridListFactorySYSU();
		void Create(Ptr<GridList>& gridList) override;
	private:
		void ReadFile(Ptr<GridList>& gridList);
	private:
		string m_fileName;
		int m_NodeNum;
		int m_BoundNodeNum;
		DArray m_NodeX, m_NodeY, m_NodeZ;
		map<int, IArray> m_NodeNeiborNodeIndex;
		IArray m_InletNodeIndex, m_OutletNodeIndex, m_WallNodeIndex;
		IArray m_InletNeiborNodeIndex, m_OutletNeiborNodeIndex, m_WallNeiborNodeIndex;
	};
}