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
#include <vector>
#include <string>
namespace zaran
{

	class InterNodeInfo
	{
	public:
		InterNodeInfo();
		~InterNodeInfo();
		void SetLocalNodeIndex(const std::vector<int>& localNodeIndex);
		void SetNeiborGridIndex(const std::vector<int>& neighborGridIndex);
		void SetNeighborNodeIndex(const std::vector<int>& neighborNodeIndex);
		void SetDataName(const std::vector<std::string>& dataName);
		std::vector<int>& GetLocalNodeIndex();
		std::vector<int>& GetNeighborGridIndex();
		std::vector<int>& GetNeighborNodeIndex();
		std::vector<std::string>& GetDataName();
	private:
		//插值节点在当前节点的位置
		std::vector<int>	localNodeIndex_;
		//插值节点对应的节点在其所在的网格上的位置
		std::vector<int> neighborNodeIndex_;
		//插值节点对应的邻居网格的索引
		std::vector<int> neighborGridIndex_;
		//插值节点需要传递变量名字列表
		std::vector<std::string> dataName_;
	};
}