#include"Read_Grid_Struct_2D.h"
#include"Log.h"
#include"Grid_Struct_2D.h"
#include<fstream>
void zaran::GridListFactoryStruct2D::Create(Ptr<GridList>& gridList)
{
	if (!gridList)
	{
		gridList = std::make_shared<GridList>();
	}
	Ptr <Grid_Struct_2D> gridPtr = std::make_shared<Grid_Struct_2D>();
	gridPtr->SetDimension(Dimension::two);
	gridPtr->SetName("Struct_Grid_2D");
	gridList->AddGrid(gridPtr);
	std::string fileName = "grid.dat";
	std::ifstream fin;
	fin.open(fileName, std::ios_base::in);
	if (!fin.is_open())
	{
		Log::error("Grid File:{} is NOT exist! Please Check!", fileName);
		system("pause");
	}
	int ghostNum;
	int ni,nj;
	fin >> ni >> nj;
	ni += ghostNum;
	nj += ghostNum;
	gridPtr->SetTotalNodeNum(ni*nj);
	gridPtr->SetNodeNum(ni,nj);
	auto& nodeTopo = gridPtr->GetNodeTopo();
	auto& nodeCoord = nodeTopo->GetCoordinate();
	nodeCoord.resize(gridPtr->GetTotalNodeNum());
	for (int j = ghostNum; j < nj + ghostNum; j++)
	{
		for (int i = ghostNum; i < ni + ghostNum ; i++)
		{
			fin >> nodeCoord[i + j * (ni + ghostNum * 2)][0] >> nodeCoord[i + j * (ni + ghostNum * 2)][1];
		}
	}


}
