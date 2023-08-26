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
#include <string>
#include <vector>
#include <unordered_map>
namespace zaran
{
	class BoundaryMap
	{
	public:
		void AddBoundary(const std::string& name, Boundary& bound);
		int GetBoundaryTypeName() { return boundaryMap_.size(); }
		std::vector<Boundary>& GetBoundary(const std::string& name);
		std::unordered_map<std::string, std::vector<Boundary>>& GetBoundaryMap() { return boundaryMap_; }
	private:
		std::unordered_map<std::string, std::vector<Boundary>>boundaryMap_;
	};
}