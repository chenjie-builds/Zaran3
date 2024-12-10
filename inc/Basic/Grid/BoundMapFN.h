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
		void AddBoundary(const string& name, BoundFN&& bound);
		void CreateBoundary(const string& name);
		int GetBoundaryTypeName() { return m_bound_map.size(); }
		dynamic_array<BoundFN>& GetBoundary(const string& name);
		map<string, dynamic_array<BoundFN>>& GetBoundaryMap() { return m_bound_map; }
	private:
		map<string, dynamic_array<BoundFN>>m_bound_map;
	};
}