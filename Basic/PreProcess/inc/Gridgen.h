#pragma once
#include "GridListFactory.h"
namespace zaran
{
	class Gridgen :public GridListFactory
	{
	public:
		void Create(Ptr<GridList>& gridList) override;
	private:
		void ReadGridgenCoord();
		void SetGridFileName(std::string& fileName);
	private:
		string	gridFileName_;
		int nBlock_;
		IArray ni_, nj_, nk_;
		Array<DArray> x_, y_, z_;

	};

}