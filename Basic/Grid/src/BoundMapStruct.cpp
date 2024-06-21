#include "BoundMapStruct.h"
#include "Log.h"
namespace zaran
{



    void BoundMapStruct::AddBoundary(const string& name, BoundStruct& bound)
    {
        auto it = m_bound_map.find(name);
        if (it == m_bound_map.end())
        {
            m_bound_map.insert({ name,std::vector<BoundStruct>({bound}) });
        }
        else
        {
            auto& bound_vec = it->second;
            for (int iBound = 0; iBound < bound_vec.size(); ++iBound)
            {
                if (bound_vec[iBound] == bound)
                {
                    return;
                }
            }
            bound_vec.emplace_back(bound);
        }
    }

    Array<BoundStruct>& zaran::BoundMapStruct::GetBoundary(const string& name)
    {
        auto it = m_bound_map.find(name);
        if (it == m_bound_map.end())
        {
            Log::info("Boundary Name:{}, is not found!", name);
        }
        return it->second;
    }
}