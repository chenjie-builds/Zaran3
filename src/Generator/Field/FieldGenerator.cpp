#include "FieldGenerator.h"
#include "GridGeneratorFlexibleSYSU.h"
#include "GridGeneratorFlexible.h"
#include "GridGeneratorBlock.h"
#include "GridGenerator.h"
#include "GridGeneratorStructGridgen.h"
#include "Log.h"
#include "NSFieldFN.h"
#include "NSFieldStruct.h"
#include "NSFieldBlockExplosion.h"
#include "NSFieldZaran.h"
#include "PolyData.h"
#include "ReadSTL.h"

namespace zaran
{
	shared_ptr<FieldManager> FieldGenerator::Create()
	{

		if (m_grid_type == GridType::Zaran)
		{
			return CreateFieldZaran();
		}
		shared_ptr<FieldManager>field_manager = make_shared<FieldManager>();
		shared_ptr < GridGenerator> grid_factory;
		if (m_grid_type == GridType::Flexible)
		{
			std::string work_dir = GlobalData::GetString("work_dir");
			std::string node_file_name = work_dir + "/node.dat";
			std::string ele_file_name = work_dir + "/cell.dat";
			std::string bnd_file_name = work_dir + "/bound.dat";
			//注意SYSU的网格文件格式不再使用
			//grid_factory = make_shared<GridBuilderSYSU_FN>(node_file_name, ele_file_name, bnd_file_name);
			grid_factory = make_shared<GridGeneratorFlexible>(node_file_name, ele_file_name, bnd_file_name);
		}
		else if (m_grid_type == GridType::Structured)
		{
			std::string work_dir = GlobalData::GetString("work_dir");
			grid_factory = make_shared<GridBuilderStructGridgen>();
		}
		else
		{
			Log::warn("Unsupported Dimension! Please Check!");
			system("pause");
		}
		dynamic_array<shared_ptr<GridBase>> grid_list;
		grid_factory->CreateGrid(grid_list);
		dynamic_array<shared_ptr<Field>> field_list;
		dynamic_array<shared_ptr<FieldDataCommInfo>> field_data_comm_info;
		field_list.resize(grid_list.size());
		field_data_comm_info.resize(grid_list.size());
		int recv_per_node = 5;
		dynamic_array<std::string> recv_data_name = { "primitive_0", "primitive_1",
												   "primitive_2", "primitive_3",
												   "primitive_4" };
		dynamic_array<index_type> recv_node_idx_local, recv_field_idx_global,
			recv_node_idx_global;

		for (index_type iField = 0; iField < grid_list.size(); iField++)
		{
			if (m_solver_type == FieldSolverType::NS_FNFDM)
			{
				auto grid = std::dynamic_pointer_cast<GridFN>(grid_list[iField]);
				field_list[iField] = make_shared < NSFieldFNFDM>(grid);
				field_list[iField]->SetIdx(iField);
			}
			else if (m_solver_type == FieldSolverType::NS_Struct)
			{
				auto grid_struct = std::dynamic_pointer_cast<GridStruct>(grid_list[iField]);
				field_list[iField] = make_shared < NSFieldStruct>(grid_struct);
				field_list[iField]->SetIdx(iField);
				int ni = grid_struct->GetNi();
				int nj = grid_struct->GetNj();
				int nk = grid_struct->GetNk();
				auto bound_map = grid_struct->GetBoundMap();
				if (bound_map->HasBound("connection")) {
					auto& connect_bound = bound_map->GetBound("connection");
					IdProxyStruct* src_idx_proxy = new IdProxyStruct(ni, nj, nk);
					recv_node_idx_local.resize(connect_bound.size());
					recv_field_idx_global.resize(connect_bound.size());
					recv_node_idx_global.resize(connect_bound.size());
					for (int i = 0; i < connect_bound.size(); i++)
					{
						index_type idx_i, idx_j, idx_k;
						connect_bound[i].GetIdx(idx_i, idx_j, idx_k);
						int idx = src_idx_proxy->GetIdx(idx_i, idx_j, idx_k);
						recv_node_idx_local[i] = idx;
						recv_field_idx_global[i] = connect_bound[i].GetTargetBlock();
						auto recv_grid_struct = std::dynamic_pointer_cast<GridStruct>(grid_list[recv_field_idx_global[i]]);
						int recv_ni = recv_grid_struct->GetNi();
						int recv_nj = recv_grid_struct->GetNj();
						int recv_nk = recv_grid_struct->GetNk();
						IdProxyStruct* tgt_idx_proxy =
							new IdProxyStruct(recv_ni, recv_nj, recv_nk);
						connect_bound[i].GetIdxTgt(idx_i, idx_j, idx_k);
						recv_node_idx_global[i] = tgt_idx_proxy->GetIdx(idx_i, idx_j, idx_k);
						delete tgt_idx_proxy;
					}
					delete src_idx_proxy;
				}
			}
			else
			{
				Log::warn("Unsupported Solver Type! Please Check!");
				system("pause");
			}
			field_list[iField]->Allocate();
			field_data_comm_info[iField] = make_shared<FieldDataCommInfo>(
				recv_node_idx_local.size(), recv_data_name, recv_node_idx_local.data(),
				recv_field_idx_global.data(), recv_node_idx_global.data());

			field_manager->AddField(field_list[iField], field_data_comm_info[iField]);
		}
		return field_manager;
	}

	void FieldGenerator::CreateGrid() {}

	void FieldGenerator::CreateField() {}

	void FieldGenerator::CreateSolver() {}

	shared_ptr<FieldManager> FieldGenerator::CreateFieldZaran()
	{
		// Read model file
		STLReader stl_reader;
		string mode_file = GlobalData::GetString("modelFileName");
		std::string work_dir = GlobalData::GetString("work_dir");
		mode_file = work_dir + "/" + mode_file;
		stl_reader.ReadSTLFile(mode_file.c_str());
		// Create model manager
		shared_ptr < PolyDataModel > poly_data_model = make_shared <PolyDataModel>();
		poly_data_model->SetPolyData(stl_reader.GetPolyData(), 1e-6);
		shared_ptr < ModelManager > model_manager = make_shared <ModelManager>();
		model_manager->AddModel(poly_data_model);
		// Create field manager
		shared_ptr<FieldManager> field_manager = make_shared<FieldManager>();
		shared_ptr<NSFieldZaran> zaran_field = make_shared<NSFieldZaran>();
		zaran_field->SetModelManager(model_manager);
		// Create master grid
		GridBlockInfo grid_info;
		grid_info.ni = GlobalData::GetInt("nx");
		grid_info.nj = GlobalData::GetInt("ny");
		grid_info.nk = GlobalData::GetInt("nz");
		grid_info.bound_box.x_min = GlobalData::GetDouble("xMin");
		grid_info.bound_box.x_max = GlobalData::GetDouble("xMax");
		grid_info.bound_box.y_min = GlobalData::GetDouble("yMin");
		grid_info.bound_box.y_max = GlobalData::GetDouble("yMax");
		grid_info.bound_box.z_min = GlobalData::GetDouble("zMin");
		grid_info.bound_box.z_max = GlobalData::GetDouble("zMax");
		grid_info.bound_type_i_minus = GlobalData::GetString("boundTypeIMinus");
		grid_info.bound_type_i_plus = GlobalData::GetString("boundTypeIPlus");
		grid_info.bound_type_j_minus = GlobalData::GetString("boundTypeJMinus");
		grid_info.bound_type_j_plus = GlobalData::GetString("boundTypeJPlus");
		grid_info.bound_type_k_minus = GlobalData::GetString("boundTypeKMinus");
		grid_info.bound_type_k_plus = GlobalData::GetString("boundTypeKPlus");
		GridGeneratorBlock* grid_factory = new GridGeneratorBlock();
		shared_ptr<GridBlock> grid;
		grid_factory->CreateGrid(grid, grid_info);
		zaran_field->SetGrid(grid);
		field_manager->AddField(zaran_field, nullptr);
		// Create slave grid
		zaran_field->CreateSlaveField(field_manager);
		zaran_field->Allocate();

		return field_manager;
	}
	shared_ptr<FieldManager> FieldGenerator::CreateFieldExplosion()
	{
		// Create master grid
		GridBlockInfo grid_info;
		grid_info.ni = 100;
		grid_info.nj = 100;
		grid_info.nk = 100;
		grid_info.bound_box.x_min = -100;
		grid_info.bound_box.x_max = 100;
		grid_info.bound_box.y_min = -100;
		grid_info.bound_box.y_max = 100;
		grid_info.bound_box.z_min = -100;
		grid_info.bound_box.z_max = 100;
		grid_info.bound_type_i_minus = GlobalData::GetString("boundTypeIMinus");
		grid_info.bound_type_i_plus = GlobalData::GetString("boundTypeIPlus");
		grid_info.bound_type_j_minus = GlobalData::GetString("boundTypeJMinus");
		grid_info.bound_type_j_plus = GlobalData::GetString("boundTypeJPlus");
		grid_info.bound_type_k_minus = GlobalData::GetString("boundTypeKMinus");
		grid_info.bound_type_k_plus = GlobalData::GetString("boundTypeKPlus");
		GridGeneratorBlock* grid_factory = new GridGeneratorBlock();
		shared_ptr<GridBlock> grid;
		grid_factory->CreateGrid(grid, grid_info);
		shared_ptr<FieldManager> field_manager = make_shared <FieldManager>();
		shared_ptr<NSFieldBlockExplosion> field = make_shared <NSFieldBlockExplosion>(grid);
		field_manager->AddField(field, nullptr);
		// Create slave grid
		field->Allocate();

		return field_manager;
	}
} // namespace zaran