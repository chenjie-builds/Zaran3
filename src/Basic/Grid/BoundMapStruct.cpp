#include "BoundStructManager.h"
#include "Log.h"
namespace zaran
{

    void BoundStructManager::AddBound(const string &name, BoundStruct &bound)
    {
        AllocateBound(name);
        auto &bound_vec = m_bound_map[name];
        for (int iBound = 0; iBound < bound_vec.size(); ++iBound)
        {
            if (bound_vec[iBound] == bound)
            {
                return;
            }
        }
        m_bound_map[name].emplace_back(bound);
    }

    void BoundStructManager::AllocateBound(const string &name)
    {
        auto it = m_bound_map.find(name);
        if (it == m_bound_map.end())
        {
            m_bound_map.insert({name, dynamic_array<BoundStruct>()});
        }
    }

    dynamic_array<BoundStruct> &BoundStructManager::GetBound(const string &name)
    {
        auto it = m_bound_map.find(name);
        if (it == m_bound_map.end())
        {
            Log::warn("Boundary Name:{}, is not found!", name);
        }
        return it->second;
    }
}