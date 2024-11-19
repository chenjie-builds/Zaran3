#pragma once
#include "BasicType.h"
#include"BoundStruct.h"
namespace zaran
{
    class BoundManagerStruct
    {
    public:
        void AddBoundary(const string& name, BoundStruct& bound);
        void CreateBoundary(const string& name);
        bool IsBoundaryExist(const string& name) { return m_bound_map.find(name) != m_bound_map.end(); }
        int GetBoundaryTypeName() { return m_bound_map.size(); }
        std::vector<BoundStruct>& GetBoundary(const string& name);
        map<string, std::vector<BoundStruct>>& GetBoundaryMap() { return m_bound_map; }
    private:
        map<string, std::vector<BoundStruct>>m_bound_map;
    };
}