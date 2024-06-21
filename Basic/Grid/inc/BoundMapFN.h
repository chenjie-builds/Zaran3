//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	BoudaryMap.h														||
//*	@brief	Store all Boundary information										||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "BoundFN.h"
#include"BasicType.h"
namespace zaran
{
	class BoundMapFN
	{
	public:
		void AddBoundary(const string& name, BoundFN& bound);
		int GetBoundaryTypeName() { return m_bound_map.size(); }
		Array<BoundFN>& GetBoundary(const string& name);
		map<string, Array<BoundFN>>& GetBoundaryMap() { return m_bound_map; }
	private:
		map<string, Array<BoundFN>>m_bound_map;
	};
}