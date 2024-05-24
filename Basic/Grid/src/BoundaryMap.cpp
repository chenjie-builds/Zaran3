#include "BoundaryMap.h"
#include "log.h"
using namespace zaran;
void BoundaryMap::AddBoundary(const std::string& name, Boundary& bound)
{
	auto it = boundaryMap_.find(name);
	if (it == boundaryMap_.end())
	{
		boundaryMap_.insert({ name,std::vector<Boundary>({bound}) });
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

std::vector<Boundary>& BoundaryMap::GetBoundary(const std::string& name)
{
	auto it = boundaryMap_.find(name);
	if (it == boundaryMap_.end())
	{
		Log::info("Boundary Name:{}, is not found!", name);
	}
	return it->second;

}

