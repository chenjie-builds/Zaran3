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
#include<map>
namespace zaran
{
	const int SYSU_INLET = 1;
	const int SYSU_OUTLET = 2;
	const int SYSU_WALL = 3;

	class GridListFactorySYSU :public GridListFactory
	{
	public:
		GridListFactorySYSU();
		void Create(std::shared_ptr<GridList>& gridList) override;
	private:
		void ReadFile();
	private:
		std::string m_fileName;
		int m_NodeNum;
		int m_BoundNodeNum;
		std::vector<double> m_NodeX, m_NodeY, m_NodeZ;
		std::map<int,std::vector<int>> m_NodeNeiborNodeIndex;
		std::vector<int> m_InletNodeIndex, m_OutletNodeIndex, m_WallNodeIndex;
		std::vector<int> m_InletNeiborNodeIndex, m_OutletNeiborNodeIndex, m_WallNeiborNodeIndex;
	};
}