#include "BoundMapFN.h"
#include "Log.h"
using namespace zaran;
void BoundMapFN::AddBoundary(const std::string& name, BoundFN& bound)
{
	auto it = m_bound_map.find(name);
	if (it == m_bound_map.end())
	{
		m_bound_map.insert({ name,std::vector<BoundFN>({bound}) });
	}
	else
	{
		auto& boundVec = it->second;
		for (int iBound = 0; iBound < boundVec.size(); ++iBound)
		{
			if (boundVec[iBound].GetIdxBound() > bound.GetIdxBound())
			{
				boundVec.insert(boundVec.begin() + iBound, bound);
				return;
			}
			if (boundVec[iBound] == bound)
			{
				return;
			}
		}
		boundVec.emplace_back(bound);
	}
}

std::vector<BoundFN>& BoundMapFN::GetBoundary(const std::string& name)
{
	auto it = m_bound_map.find(name);
	if (it == m_bound_map.end())
	{
		Log::info("Boundary Name:{}, is not found!", name);
	}
	return it->second;

}

