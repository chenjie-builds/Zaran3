#include"GridGeneratorBlockExplosion.h"
#include<fstream>
namespace zaran
{
	GridGeneratorBlockExplosion::GridGeneratorBlockExplosion()
	{
	}
	GridGeneratorBlockExplosion::~GridGeneratorBlockExplosion()
	{

	}
	void GridGeneratorBlockExplosion::CreateGrid(shared_ptr<GridBlock> grid, GridBlockInfo& grid_info)
	{
		GridGeneratorBlock::CreateGrid(grid, grid_info);
		SetIBlankByBuilding();
		ExportTecplot();
	}

	void GridGeneratorBlockExplosion::SetIBlankByBuilding()
	{

	}

	void GridGeneratorBlockExplosion::ExportTecplot()
	{
		auto grid = GetGrid();
		int ni = grid->GetNi();
		int nj = grid->GetNj();
		int nk = grid->GetNk();
		auto node = grid->GetNode();
		IdProxyStruct& idx_proxy = grid->GetIdxProxy();
		std::ofstream fout("block_explosion.dat");
		fout << "TITLE = \"block_explosion\"" << std::endl;
		fout << "VARIABLES = \"X\", \"Y\", \"Z\", \"IBLANK\"" << std::endl;
		fout << "ZONE I = " << ni << ", J = " << nj << ", K = " << nk << ", F = POINT" << std::endl;
		for (int k = 0; k < nk; ++k)
		{
			for (int j = 0; j < nj; ++j)
			{
				for (int i = 0; i < ni; ++i)
				{
					auto id = idx_proxy(i, j, k);
					auto coord = node->GetCoord(id);
					fout << coord[0] << " " << coord[1] << " " << coord[2] << " " << int(grid->GetIBlank(i, j, k)) << std::endl;
				}
			}
		}
	}


}