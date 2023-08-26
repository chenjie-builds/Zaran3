#include "Gridgen.h"
#include <fstream>
namespace zaran
{


	void Gridgen::ReadGridgenCoord()
	{

	}

	void Gridgen::SetGridFileName(std::string& fileName)
	{
		gridFileName_ = fileName;
	}

	void Gridgen::Create(std::shared_ptr<GridList>& gridList)
	{

		std::fstream fin;
		fin.open(gridFileName_, std::ios_base::in);
		fin >> nBlock_;
		if (!gridList)
		{
			gridList = std::make_shared<GridList>();
		}
		gridList->Resize(nBlock_);
		for (size_t i = 0; i < gridList->GetGridNumber(); i++)
		{
			auto& currentGridPtr = gridList->GetGrid(i);
			currentGridPtr = std::make_shared<Grid>();
		}
	}

}