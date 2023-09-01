#pragma once
#include"grid.h"
#include "log.h"
using namespace zaran;

Grid::Grid()
{
	GridBase();
	boundaryMap_ = std::make_shared<BoundaryMap>();
	interNodeInfo_ = std::make_shared<InterNodeInfo>();	
}

Grid::Grid(const Grid& grid)
{

}

void Grid::SetLevel(const int& level)
{
	level_ = level;
}


void Grid::SetInnerNodeNumber(const int& innerNodeNum)
{
	nInnerNode_ = innerNodeNum;
}

void Grid::SetBoundNodeNumber(const int& boundNodeNum)
{
	nBoundNode_ = boundNodeNum;
}

void Grid::SetTotalNodeNumber(const int& totalNodeNum)
{
	nTotalNode_ = totalNodeNum;
}

void Grid::SetNodeTopo(NodeTopoInfoVec& nodeTopo)
{
	nodeTopo_ = nodeTopo;
}

void Grid::SetFaceTopo(FaceTopoInfoVec& faceTopo)
{
	faceTopo_ = faceTopo;
}

void Grid::SetCellTopoInfo(CellTopoInfoVec& cellTopo)
{
	cellTopo_ = cellTopo;
}

void Grid::SetInterNodeInfo(Ptr<InterNodeInfo> interNodeInfo)
{
	interNodeInfo_ = interNodeInfo;
}


void Grid::SetBoundaryMap(Ptr<BoundaryMap> boundaryMap)
{
	boundaryMap_ = boundaryMap;
}

const int& Grid::GetLevel() const
{
	return level_;
}


const int& Grid::GetTotalNodeNum() const
{
	return nTotalNode_;
}


const int& Grid::GetInnerNodeNum() const
{
	return nInnerNode_;
}

const int& Grid::GetBoundNodeNum() const
{
	return nBoundNode_;
}

zaran::CellTopoInfoVec& Grid::GetCellTopoInfo()
{
	return cellTopo_;
}

zaran::FaceTopoInfoVec& Grid::GetFaceTopoInfo()
{
	return faceTopo_;
}

InterNodeInfoPtr& Grid::GetInterNodeInfo()
{
	return interNodeInfo_;
}

NodeTopoInfoVec& Grid::GetNodeTopoInfo()
{
	return nodeTopo_;
}


BoundaryMapPtr& Grid::GetBoundaryMap()
{
	return boundaryMap_;
}

