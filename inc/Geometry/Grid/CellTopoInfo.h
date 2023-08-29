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
#include"BasicType.h"
namespace zaran
{
	// 单元拓扑类，记录单元包含的节点和面元
	class CellTopoInfo
	{
	public:
		CellTopoInfo();
		~CellTopoInfo();
		void SetNode(IArray& node);
		void SetFace(IArray& face);
		IArray& GetFace();
		IArray& GetNode();
	private:
		// 单元包含的节点
		IArray node_;
		// 单元包含的面元
		IArray face_;
	};
}