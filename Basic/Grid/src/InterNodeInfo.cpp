#include "InterNodeInfo.h"
using namespace zaran;
void InterNodeInfo::SetLocalNodeIndex(const IArray& localNodeIndex)
{
	localNodeIndex_ = localNodeIndex;
}

void InterNodeInfo::SetNeiborGridIndex(const IArray& neighborGridIndex)
{
	neighborGridIndex_ = neighborGridIndex;
}

void InterNodeInfo::SetNeighborNodeIndex(const IArray& neighborNodeIndex)
{
	neighborNodeIndex_ = neighborNodeIndex;
}

IArray& InterNodeInfo::GetLocalNodeIndex()
{
	return localNodeIndex_;
}

IArray& InterNodeInfo::GetNeighborGridIndex()
{
	return neighborGridIndex_;
}

IArray& InterNodeInfo::GetNeighborNodeIndex()
{
	return neighborNodeIndex_;
}

void InterNodeInfo::SetDataName(const SArray& dataName)
{
	dataName_ = dataName;
}

SArray& InterNodeInfo::GetDataName()
{
	return dataName_;
}

InterNodeInfo::InterNodeInfo()
{

}

InterNodeInfo::~InterNodeInfo()
{

}

