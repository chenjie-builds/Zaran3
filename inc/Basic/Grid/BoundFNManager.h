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
	class BoundManagerFN
	{
	public:
		void AddBoundary(const string& name, BoundFN&& bound);
		void AllocateBound(const string& name);
		count_type GetBoundCount() const { return m_bound_map.size(); }
		dynamic_array<BoundFN>& GetBound(const string& name);
		map<string, dynamic_array<BoundFN>>& GetBoundMap() { return m_bound_map; }
	private:
		map<string, dynamic_array<BoundFN>>m_bound_map;
	};
}