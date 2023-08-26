//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	CelTopoInfo.h														||
//*	@brief	Cell Topology setting												||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include <vector>
namespace zaran
{
	// 单元拓扑类，记录单元包含的节点和面元
	class CellTopoInfo
	{
	public:
		CellTopoInfo();
		~CellTopoInfo();
		void SetNode(std::vector<int>& node);
		void SetFace(std::vector<int>& face);
		std::vector<int>& GetFace();
		std::vector<int>& GetNode();
	private:
		// 单元包含的节点
		std::vector<int> node_;
		// 单元包含的面元
		std::vector<int>face_;
	};
}