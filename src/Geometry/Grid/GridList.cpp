#include "GridList.h"
#include "log.h"
namespace zaran
{


	GridList::GridList(std::vector<std::shared_ptr<Grid>>& gridList)
	{
		gridList_ = gridList;
	}

	GridList::~GridList()
	{

	}

	void GridList::AddGrid(std::shared_ptr<Grid>& grid)
	{
		gridList_.emplace_back(grid);
	}

	void GridList::CheckGrid()
	{
		for (int iGrid = 0; iGrid < gridList_.size(); ++iGrid)
		{
			if (iGrid != gridList_[iGrid]->GetIndex())
			{
				ZaranLog::warn("the index of grid are not match to grid set!");
				ZaranLog::warn("grid index:{}, grid set index:{}", gridList_[iGrid]->GetIndex(), iGrid);
			}
		}
	}

	std::shared_ptr<Grid>& GridList::GetGrid(const int index)
	{
		return gridList_[index];
	}

	std::vector<std::shared_ptr<zaran::Grid>>& GridList::GetGrid()
	{
		return gridList_;
	}

	void GridList::Resize(int gridNum)
	{
		gridList_.resize(gridNum);
	}

}