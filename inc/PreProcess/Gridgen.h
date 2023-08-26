#pragma once
#include "GridListFactory.h"
namespace zaran
{
	class Gridgen:public GridListFactory
	{
	public:
		void Create(std::shared_ptr<GridList>& gridList) override;
	private:
		void ReadGridgenCoord();
		void SetGridFileName(std::string& fileName);
	private:
		std::string	gridFileName_;
		int nBlock_;
		std::vector<int> ni_, nj_, nk_;
		std::vector<std::vector<double>> x_, y_, z_;

	};

}