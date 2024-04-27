#include"Read_Grid_Struct_2D.h"
#include"Log.h"
#include"Grid_Struct_2D.h"
#include<fstream>
void zaran::GridFactoryStruct2D::Create(Grid*& grid)
{
	if (grid != nullptr)
	{
		delete[] grid;
	}
	Grid_Struct_2D* gridPtr = new Grid_Struct_2D();
	gridPtr->SetDimension(Dimension::two);
	gridPtr->SetName("Struct_Grid_2D");
	gridPtr->SetType(GridType::Structured);
	grid = dynamic_cast<Grid*>(gridPtr);
	std::string fileName = "grid.dat";
	std::ifstream fin;
	fin.open(fileName, std::ios_base::in);
	if (!fin.is_open())
	{
		Log::error("Grid File:{} is NOT exist! Please Check!", fileName);
		system("pause");
	}
	int ghostNum;
	int ni, nj;
	fin >> ni >> nj;
	ni += ghostNum;
	nj += ghostNum;
	gridPtr->SetTotalNodeNum(ni * nj);
	gridPtr->SetNodeNum(ni, nj);
	NodeTopo* nodeTopo = gridPtr->GetNodeTopo();
	auto& nodeCoord = nodeTopo->GetCoordinate();
	nodeCoord.resize(gridPtr->GetTotalNodeNum());
	for (int j = ghostNum; j < nj + ghostNum; j++)
	{
		for (int i = ghostNum; i < ni + ghostNum; i++)
		{
			fin >> nodeCoord[i + j * (ni + ghostNum * 2)][0] >> nodeCoord[i + j * (ni + ghostNum * 2)][1];
		}
	}


}
