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
#include "Boundary.h"
#include"BasicType.h"
namespace zaran
{
	class BoundaryMap
	{
	public:
		void AddBoundary(const string& name, Boundary& bound);
		int GetBoundaryTypeName() { return boundaryMap_.size(); }
		vector<Boundary>& GetBoundary(const string& name);
		map<string, vector<Boundary>>& GetBoundaryMap() { return boundaryMap_; }
	private:
		map<string, vector<Boundary>>boundaryMap_;
	};
}