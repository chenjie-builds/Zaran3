#include "BoundFNManager.h"
#include "Log.h"
using namespace zaran;
void BoundManagerFN::AddBoundary(const std::string& name, BoundFN&& bound)
{
	auto it = m_bound_map.find(name);
	if (it == m_bound_map.end())
	{
		m_bound_map.insert({ name,dynamic_array<BoundFN>({bound}) });
	}
	else
	{
		auto& bound_vec = it->second;
		for (index_type iter_bound = 0; iter_bound < bound_vec.size(); ++iter_bound)
		{
			if (bound_vec[iter_bound].GetIdxBound() > bound.GetIdxBound())
			{
				bound_vec.insert(bound_vec.begin() + iter_bound, bound);
				return;
			}
			if (bound_vec[iter_bound] == bound)
			{
				return;
			}
		}
		bound_vec.emplace_back(bound);
	}
}

void zaran::BoundManagerFN::AllocateBound(const string &name)
{
	auto it = m_bound_map.find(name);
	if (it == m_bound_map.end())
	{
		m_bound_map.insert({ name,dynamic_array<BoundFN>() });
	}
}

dynamic_array<BoundFN>& BoundManagerFN::GetBound(const std::string& name)
{
	auto it = m_bound_map.find(name);
	if (it == m_bound_map.end())
	{
		Log::info("Boundary Name:{}, is not found!", name);
	}
	return it->second;

}

