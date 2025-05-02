/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file BoundStructManager.h
 * \brief Boundary Information Manager for Structured Grid
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "BasicType.h"
#include"BoundStruct.h"
namespace zaran
{
    class BoundStructManager
    {
    public:
        void AddBound(const string& name, BoundStruct& bound);
        void AllocateBound(const string& name);
        bool HasBound(const string& name)const { return m_bound_map.find(name) != m_bound_map.end(); }
        count_type GetBoundCount()const { return m_bound_map.size(); }
        dynamic_array<BoundStruct>& GetBound(const string& name);
        map<string, dynamic_array<BoundStruct>>& GetBoundMap() { return m_bound_map; }
    private:
        map<string, dynamic_array<BoundStruct>>m_bound_map;
    };
}