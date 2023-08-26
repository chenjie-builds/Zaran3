#include "InterNodeInfo.h"
namespace zaran
{
	void InterNodeInfo::SetLocalNodeIndex(const std::vector<int>& localNodeIndex)
	{
		localNodeIndex_ = localNodeIndex;
	}

	void InterNodeInfo::SetNeiborGridIndex(const std::vector<int>& neighborGridIndex)
	{
		neighborGridIndex_ = neighborGridIndex;
	}

	void InterNodeInfo::SetNeighborNodeIndex(const std::vector<int>& neighborNodeIndex)
	{
		neighborNodeIndex_ = neighborNodeIndex;
	}

	std::vector<int>& InterNodeInfo::GetLocalNodeIndex()
	{
		return localNodeIndex_;
	}

	std::vector<int>& InterNodeInfo::GetNeighborGridIndex()
	{
		return neighborGridIndex_;
	}

	std::vector<int>& InterNodeInfo::GetNeighborNodeIndex()
	{
		return neighborNodeIndex_;
	}

	void InterNodeInfo::SetDataName(const std::vector<std::string>& dataName)
	{
		dataName_ = dataName;
	}

	std::vector<std::string>& InterNodeInfo::GetDataName()
	{
		return dataName_;
	}

	InterNodeInfo::InterNodeInfo()
	{

	}

	InterNodeInfo::~InterNodeInfo()
	{

	}

}