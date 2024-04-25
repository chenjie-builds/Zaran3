#pragma once
#include "GridListFactory.h"
namespace zaran
{
	class Gridgen :public GridListFactory
	{
	public:
		void Create(Grid*& grid) override;
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