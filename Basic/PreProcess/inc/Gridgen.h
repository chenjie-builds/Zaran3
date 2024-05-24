#pragma once
#include "GridBuilder.h"
namespace zaran
{
	class Gridgen :public GridCreater
	{
	public:
		Gridgen();
		~Gridgen();
		GridBase* CreateGrid() override;
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