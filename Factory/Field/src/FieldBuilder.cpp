#include "FieldBuilder.h"
#include "GridFactory.h"
#include "FNGridFactory.h"
#include "GridStructFactoryGridgen.h"
#include "GridBlockFactory.h"
#include "NSFieldFN.h"
#include "NSFieldStruct.h"
#include "Log.h"
#include "NSFieldZaran.h"
#include "ReadSTL.h"
#include "PolyData.h"
namespace zaran
{
    FieldManager *FieldBuilder::Create()
    {
        return CreateFieldZaran();
        if(m_grid_type==GridType::Zaran)
        {
            return CreateFieldZaran();
        }
        FieldManager *field_manager = new FieldManager();
        GridFactory *grid_factory;
        if (m_grid_type == GridType::Flexible)
        {
            grid_factory = new GridFNFactorySYSU();
        }
        else if (m_grid_type == GridType::Structured)
        {
            grid_factory = new GridStructFactoryGridgen();
        }
        else
        {
            Log::warn("Unsupported Dimension! Please Check!");
            system("pause");
        }
        GridBase **grid;
        int grid_num;
        grid_factory->CreateGrid(grid, grid_num);
        Field **field_list;
        FieldDataCommInfo **field_data_comm_info;
        field_list = new Field *[grid_num];
        field_data_comm_info = new FieldDataCommInfo *[grid_num];
        int recv_per_node = 5;
        std::vector<std::string> recv_data_name = {"density", "x_velocity", "y_velocity", "z_velocity", "pressure"};
        std::vector<int> recv_node_idx_local, recv_field_idx_global, recv_node_idx_global;

        for (int iField = 0; iField < grid_num; iField++)
        {
            if (m_solver_type == FieldSolverType::NS_FNFDM)
            {
                field_list[iField] = new NSFieldFNFDM(grid[iField]);
            }
            else if (m_solver_type == FieldSolverType::NS_Struct)
            {
                field_list[iField] = new NSFieldStruct(grid[iField]);
                auto grid_struct = dynamic_cast<GridStruct *>(grid[iField]);
                auto bound_map = grid_struct->GetBoundMap();
                auto &connect_bound = bound_map->GetBoundary("connection");
                StructIdxProxy *src_idx_proxy = new StructIdxProxy(grid_struct);
                recv_node_idx_local.resize(connect_bound.size());
                recv_field_idx_global.resize(connect_bound.size());
                recv_node_idx_global.resize(connect_bound.size());
                for (int i = 0; i < connect_bound.size(); i++)
                {
                    int idx_i, idx_j, idx_k;
                    connect_bound[i].GetIdx(idx_i, idx_j, idx_k);
                    int idx = src_idx_proxy->GetIdx(idx_i, idx_j, idx_k);
                    recv_node_idx_local[i] = idx;
                    recv_field_idx_global[i] = connect_bound[i].GetTargetBlock();
                    StructIdxProxy *tgt_idx_proxy = new StructIdxProxy(dynamic_cast<GridStruct *>(grid[recv_field_idx_global[i]]));
                    connect_bound[i].GetIdxTgt(idx_i, idx_j, idx_k);
                    recv_node_idx_global[i] = tgt_idx_proxy->GetIdx(idx_i, idx_j, idx_k);
                    delete tgt_idx_proxy;
                }
                delete src_idx_proxy;
            }
            else
            {
                Log::warn("Unsupported Solver Type! Please Check!");
                system("pause");
            }
            field_list[iField]->Allocate();
            field_data_comm_info[iField] = new FieldDataCommInfo(recv_per_node, recv_node_idx_local.size(), recv_data_name, recv_node_idx_local.data(), recv_field_idx_global.data(), recv_node_idx_global.data());
            field_manager->AddField(field_list[iField], field_data_comm_info[iField]);
        }
        delete grid_factory;
        return field_manager;
    }

    void FieldBuilder::CreateGrid()
    {
    }

    void FieldBuilder::CreateField()
    {
    }

    void FieldBuilder::CreateSolver()
    {
    }

    FieldManager *FieldBuilder::CreateFieldZaran()
    {
        // Read model file
        STLReader stl_reader;
        string mode_file = GlobalData::GetString("modelFileName");
        stl_reader.ReadSTLFile(mode_file.c_str());
        // Create model manager
        PolyDataModel *poly_data_model = new PolyDataModel();
        poly_data_model->SetPolyData(stl_reader.GetPolyData(), 0.0001);
        ModelManager *model_manager = new ModelManager();
        model_manager->AddModel(poly_data_model);
        // Create field manager
        FieldManager *field_manager = new FieldManager();
        NSFieldZaran *zaran_field = new NSFieldZaran();
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
        grid_info.bound_type_i_minus=GlobalData::GetString("boundTypeIMinus");
        grid_info.bound_type_i_plus=GlobalData::GetString("boundTypeIPlus");
        grid_info.bound_type_j_minus=GlobalData::GetString("boundTypeJMinus");
        grid_info.bound_type_j_plus=GlobalData::GetString("boundTypeJPlus");
        grid_info.bound_type_k_minus=GlobalData::GetString("boundTypeKMinus");
        grid_info.bound_type_k_plus=GlobalData::GetString("boundTypeKPlus");
        GridBlockFactory *grid_factory = new GridBlockFactory();
        GridBlock *grid;
        grid_factory->CreateGrid(grid, grid_info);
        zaran_field->SetGrid(grid);
        field_manager->AddField(zaran_field, nullptr);
        // Create slave grid
        zaran_field->CreateSlaveField(field_manager);
        zaran_field->Allocate();

        return field_manager;
    }
}