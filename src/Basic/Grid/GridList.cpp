#include "GridList.h"
#include "Log.h"
using namespace zaran;
GridList::GridList(dynamic_array<shared_ptr<GridBase>>& gridList)
{
	gridList_ = gridList;
}

GridList::~GridList()
{

}

void GridList::AddGrid(shared_ptr<GridBase> grid)
{
	gridList_.emplace_back(grid);
}

void GridList::CheckGrid()
{
	for (int iGrid = 0; iGrid < gridList_.size(); ++iGrid)
	{
		if (iGrid != gridList_[iGrid]->GetIndex())
		{
			Log::warn("the index of grid are not match to grid set!");
			Log::warn("grid index:{}, grid set index:{}", gridList_[iGrid]->GetIndex(), iGrid);
		}
	}
}

shared_ptr<GridBase>& GridList::GetGrid(const int index)
{
	return gridList_[index];
}

dynamic_array<shared_ptr<zaran::GridBase>>& GridList::GetGrid()
{
	return gridList_;
}

void GridList::Resize(int gridNum)
{
	gridList_.resize(gridNum);
}
