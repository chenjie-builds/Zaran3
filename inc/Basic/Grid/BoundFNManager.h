/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file BoundaryFN.h
 * \brief Boundary Information for Flexible Node FDM.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "BoundFN.h"
#include "BasicType.h"
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