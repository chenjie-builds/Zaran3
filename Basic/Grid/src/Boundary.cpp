#include "Boundary.h"
using namespace zaran;
Boundary::Boundary()
{
	nodeIndex_ = -1;
	ghostNodeIndex_ = -1;
	innerNodeIndex_ = -1;
	norm_ = { 0,0,0 };
}

Boundary::Boundary(int nodeIndex, int innerNodeIndex, int ghostNodeIndex, DVector3D& norm)
{
	nodeIndex_ = nodeIndex;
	innerNodeIndex_ = innerNodeIndex;
	ghostNodeIndex_ = ghostNodeIndex;
	norm_ = norm;
}

void Boundary::SetGhostNodeIndex(const int& ghostNodeIndex)
{
	ghostNodeIndex_ = ghostNodeIndex;
}

void Boundary::SetInnerNodeIndex(const int& innnerNodeIndex)
{
	innerNodeIndex_ = innnerNodeIndex;
}

void Boundary::SetNodeIndex(const int& nodeIndex)
{
	nodeIndex_ = nodeIndex;
}
void Boundary::SetNorm(const DVector3D& boundNorm)
{
	norm_ = boundNorm;
}

bool Boundary::operator==(const Boundary& bound)
{
	return nodeIndex_ == bound.nodeIndex_ && ghostNodeIndex_ == bound.ghostNodeIndex_ && innerNodeIndex_ == bound.innerNodeIndex_ && norm_ == bound.norm_;
}

