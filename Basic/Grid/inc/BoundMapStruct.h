#pragma once
#include "BasicType.h"
#include"BoundStruct.h"
namespace zaran
{
    class BoundMapStruct
    {
    public:
        void AddBoundary(const string& name, BoundStruct& bound);
        int GetBoundaryTypeName() { return m_bound_map.size(); }
        Array<BoundStruct>& GetBoundary(const string& name);
        map<string, Array<BoundStruct>>& GetBoundaryMap() { return m_bound_map; }
    private:
        map<string, Array<BoundStruct>>m_bound_map;
    };
}