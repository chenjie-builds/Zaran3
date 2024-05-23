//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	InterNodeInfo.h														||
//*	@brief	插值节点信息，用于网格间传递数据										||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include"BasicType.h"
namespace zaran
{

	class InterNode
	{
	public:
		InterNode();
		~InterNode();
		void SetLocalNodeIndex(const IArray& localNodeIndex);
		void SetNeiborGridIndex(const IArray& neighborGridIndex);
		void SetNeighborNodeIndex(const IArray& neighborNodeIndex);
		void SetDataName(const SArray& dataName);
		IArray& GetLocalNodeIndex();
		IArray& GetNeighborGridIndex();
		IArray& GetNeighborNodeIndex();
		SArray& GetDataName();
	private:
		//插值节点在当前节点的位置
		IArray	localNodeIndex_;
		//插值节点对应的节点在其所在的网格上的位置
		IArray neighborNodeIndex_;
		//插值节点对应的邻居网格的索引
		IArray neighborGridIndex_;
		//插值节点需要传递变量名字列表
		SArray dataName_;
	};
}