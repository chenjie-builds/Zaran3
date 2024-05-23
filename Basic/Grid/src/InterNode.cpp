#include "InterNode.h"
using namespace zaran;
void InterNode::SetLocalNodeIndex(const IArray& localNodeIndex)
{
	localNodeIndex_ = localNodeIndex;
}

void InterNode::SetNeiborGridIndex(const IArray& neighborGridIndex)
{
	neighborGridIndex_ = neighborGridIndex;
}

void InterNode::SetNeighborNodeIndex(const IArray& neighborNodeIndex)
{
	neighborNodeIndex_ = neighborNodeIndex;
}

IArray& InterNode::GetLocalNodeIndex()
{
	return localNodeIndex_;
}

IArray& InterNode::GetNeighborGridIndex()
{
	return neighborGridIndex_;
}

IArray& InterNode::GetNeighborNodeIndex()
{
	return neighborNodeIndex_;
}

void InterNode::SetDataName(const SArray& dataName)
{
	dataName_ = dataName;
}

SArray& InterNode::GetDataName()
{
	return dataName_;
}

InterNode::InterNode()
{

}

InterNode::~InterNode()
{

}

