//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	FNFDM3D.h															||
//*	@brief	读取三维自由节点有限差分网格											||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include"GridListFactory.h"
namespace zaran
{
	class GridFactoryFNFDM3D :public GridListFactory
	{
	public:
		GridFactoryFNFDM3D();
		void Create(Grid*& gridList) override;
	private:
		void ReadFile(Grid* gridList);
		void SortNeiborNode(Grid* gridList);
		/// @brief 扩展邻居节点，使用kd树，找到范围内的节点
		/// @param gridList 
		void ExtendNeighborNode(Grid* gridList);
		void ReadCellFile(Grid* gridList);
		void ReadBoundFaceFile(Grid* gridList);
	private:
		string m_node_file_name;
		string m_ele_file_name;
		string m_bnd_file_name;
		int m_NodeNum;
		int m_BoundNodeNum;
		DArray m_NodeX, m_NodeY, m_NodeZ;
		map<int, IArray> m_NodeNeiborNodeIndex;
		IArray m_InletNodeIndex, m_OutletNodeIndex, m_WallNodeIndex;
		IArray m_InletNeiborNodeIndex, m_OutletNeiborNodeIndex, m_WallNeiborNodeIndex;
	};
}