#include "FieldGeneratorBuildingExplosion.h"
#include <omp.h>
#include <fstream>
#include "Log.h"
namespace zaran
{

	FieldGeneratorBuildingExplosion::FieldGeneratorBuildingExplosion(GridType grid_type, FieldSolverType solver_type, Dimension dim)
		:FieldGenerator(GridType::Block, FieldSolverType::NS_Struct, Dimension::three)
	{

	}

	std::shared_ptr<FieldManager> FieldGeneratorBuildingExplosion::Create()
	{
		CreateGrid();
		TagGrid();
		Log::info("Tag Grid!");
		OutTest();
		CreateField();
		return m_field_manager;
	}



	void FieldGeneratorBuildingExplosion::CreateGrid()
	{

		GridBlockInfo grid_info;
		grid_info.ni = 101;
		grid_info.nj = 101;
		grid_info.nk = 101;
		grid_info.bound_box.x_min = -0;
		grid_info.bound_box.x_max = 50;
		grid_info.bound_box.y_min = -0;
		grid_info.bound_box.y_max = 50;
		grid_info.bound_box.z_min = -0;
		grid_info.bound_box.z_max = 50;
		grid_info.bound_type_i_minus = GlobalData::GetString("boundTypeIMinus");
		grid_info.bound_type_i_plus = GlobalData::GetString("boundTypeIPlus");
		grid_info.bound_type_j_minus = GlobalData::GetString("boundTypeJMinus");
		grid_info.bound_type_j_plus = GlobalData::GetString("boundTypeJPlus");
		grid_info.bound_type_k_minus = GlobalData::GetString("boundTypeKMinus");
		grid_info.bound_type_k_plus = GlobalData::GetString("boundTypeKPlus");
		GridGeneratorBlock* grid_factory = new GridGeneratorBlock();
		grid_factory->CreateGrid(m_grid, grid_info);

	}

	void FieldGeneratorBuildingExplosion::TagGrid()
	{
		auto grid = std::static_pointer_cast<GridBlock>(m_grid);
		auto idx_proxy = grid->GetIdxProxy();
		auto node = grid->GetNode();
		int ni = grid->GetNi();
		int nj = grid->GetNj();
		int nk = grid->GetNk();
		for (int k = 0; k < nk; k++)
		{
			for (int j = 0; j < nj; j++)
			{
				for (int i = 0; i < ni; i++)
				{
					auto id = idx_proxy->GetIdx(i, j, k);
					auto iblank = grid->GetIBlank(i, j, k);
					iblank = IBlank::Trans;
					auto coord = node->GetCoord(id);
					if (abs(coord[2] - 0.5) < SMALL_NUMBER || abs(coord[2] - 10.5) < SMALL_NUMBER || abs(coord[2] - 20.5) < SMALL_NUMBER || abs(coord[2] - 30.5) < SMALL_NUMBER || abs(coord[2] - 40.5) < SMALL_NUMBER || abs(coord[2] - 50.5) < SMALL_NUMBER || abs(coord[2] - 60.5) < SMALL_NUMBER || abs(coord[2] - 70.5) < SMALL_NUMBER || abs(coord[2] - 80.5) < SMALL_NUMBER || abs(coord[2] - 90.5) < SMALL_NUMBER)
					{
						if (coord[0] < 40 || coord[0] > 101 || coord[1] < 40 || coord[1] > 121)
							iblank = IBlank::Fluid;
						else
							iblank = IBlank::Solid;
					}
					else if (abs(coord[1] - 40.5) < SMALL_NUMBER || abs(coord[1] - 75.5) < SMALL_NUMBER || abs(coord[1] - 85.5) < SMALL_NUMBER || abs(coord[1] - 120.5) < SMALL_NUMBER)
					{
						if (coord[0] < 40 || coord[0] > 101 || coord[2] < 0 || coord[2] > 91)
							iblank = IBlank::Fluid;
						else
							iblank = IBlank::Solid;
					}
					else if (abs(coord[0] - 40.5) < SMALL_NUMBER || abs(coord[0] - 100.5) < SMALL_NUMBER)
					{
						if (coord[1] < 40 || coord[1] > 121 || coord[2] < 0 || coord[2] > 91)
							iblank = IBlank::Fluid;
						else
							iblank = IBlank::Solid;
					}
					else if (abs(coord[0] - 70.5) < SMALL_NUMBER)
					{
						if (coord[1] < 40 || coord[1] > 121 || (coord[1] > 75 && coord[1] < 86) || coord[2] < 0 || coord[2] > 91)
							iblank = IBlank::Fluid;
						else
							iblank = IBlank::Solid;
					}
					else
						iblank = IBlank::Fluid;
					if (abs(coord[0] - 100.5) < SMALL_NUMBER)
					{
						if (coord[1] < 60.5 && coord[1] > 55.5 || coord[1] > 100.5 && coord[1] < 105.5)
						{
							if (coord[2] > 2.5 && coord[2] < 7.5 ||
								coord[2] > 12.5 && coord[2] < 17.5 ||
								coord[2] > 22.5 && coord[2] < 27.5 ||
								coord[2] > 32.5 && coord[2] < 37.5 ||
								coord[2] > 42.5 && coord[2] < 47.5 ||
								coord[2] > 52.5 && coord[2] < 57.5 ||
								coord[2] > 62.5 && coord[2] < 67.5 ||
								coord[2] > 72.5 && coord[2] < 77.5 ||
								coord[2] > 82.5 && coord[2] < 87.5)
								iblank = IBlank::Fluid;
						}
					}

					grid->SetIBlank(i, j, k, iblank);
				}
			}
		}
	}
	void FieldGeneratorBuildingExplosion::OutTest()
	{
		//out tecplot data format
		std::ofstream out("test.dat");
		auto grid = std::static_pointer_cast<GridBlock>(m_grid);
		auto idx_proxy = grid->GetIdxProxy();
		auto node = grid->GetNode();
		int ni = grid->GetNi();
		int nj = grid->GetNj();
		int nk = grid->GetNk();
		out << "TITLE = \"test\"" << std::endl;
		out << "VARIABLES = \"X\", \"Y\", \"Z\", \"IBlank\"" << std::endl;
		out << "ZONE I=" << ni << ", J=" << nj << ", K=" << nk << ", F=POINT" << std::endl;
		for (int k = 0; k < nk; k++)
		{
			for (int j = 0; j < nj; j++)
			{
				for (int i = 0; i < ni; i++)
				{
					auto id = idx_proxy->GetIdx(i, j, k);
					auto coord = node->GetCoord(id);
					out << coord[0] << " " << coord[1] << " " << coord[2] << " " << int(grid->GetIBlank(i, j, k)) << std::endl;
				}
			}
		}
		out.close();

	}
	void FieldGeneratorBuildingExplosion::CreateField()
	{
		m_field_manager = std::make_shared<FieldManager>();
		std::shared_ptr<NSFieldBlockExplosion> field = std::make_shared<NSFieldBlockExplosion>(m_grid);
		m_field_manager->AddField(field, nullptr);
		field->Allocate();
	}

}