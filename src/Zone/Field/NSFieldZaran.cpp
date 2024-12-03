#include "NSFieldZaran.h"
#include "GlobalData.h"
#include "Log.h"
#include <fstream>
#include <omp.h>
#include "GridFNFactoryZaran.h"
#include "NSSolverBlock.h"
namespace zaran
{
    NSFieldZaran::NSFieldZaran() : FieldNS(nullptr, FieldType::NS_Zaran)
    {
        m_slave_field = nullptr;
        m_model_manager = nullptr;
        m_idx_proxy = nullptr;
    }

    NSFieldZaran::~NSFieldZaran()
    {
        if (m_idx_proxy != nullptr)
            delete m_idx_proxy;
        if (m_model_manager != nullptr)
            delete m_model_manager;
    }

    void NSFieldZaran::Allocate()
    {
        FieldNS::Allocate();
        AllocateIdxProxy();
    }

    void NSFieldZaran::SetModelManager(ModelManager *model_manager)
    {
        m_model_manager = model_manager;
    }

    void NSFieldZaran::AllocateSolver()
    {
        if (m_solver != nullptr)
        {
            delete m_solver;
        }
        auto para = GetSolverPara();
        m_solver = new NSSolverBlock(GetIdx(), "NS_Block", GetSolverPara(), GetGrid(), GetDataManager());
    }

    void NSFieldZaran::AllocateDataManager()
    {
        if (m_data_manager != nullptr)
        {
            delete m_data_manager;
        }
        int ni = GetGrid()->GetNi();
        int nj = GetGrid()->GetNj();
        int nk = GetGrid()->GetNk();
        m_data_manager = new DataManagerNSStruct(GetFieldData(), ni, nj, nk);
        m_data_manager->CreateData();
        m_data_manager->RegisterData();
    }

    void NSFieldZaran::AllocateIdxProxy()
    {
        if (m_idx_proxy != nullptr)
        {
            delete m_idx_proxy;
        }
        int ni = GetGrid()->GetNi();
        int nj = GetGrid()->GetNj();
        int nk = GetGrid()->GetNk();

        m_idx_proxy = new IdxStruct(ni, nj, nk);
    }

    void NSFieldZaran::AllocateSolverPara()
    {
        if (m_solver_para != nullptr)
        {
            delete m_solver_para;
        }
        m_solver_para = new FlowSolverParamStruct();
        GetSolverPara()->Init();
    }

    GridBlock *NSFieldZaran::GetGrid()
    {
        return static_cast<GridBlock *>(Field::GetGrid());
    }
    FlowSolverParamStruct *NSFieldZaran::GetSolverPara()
    {
        return static_cast<FlowSolverParamStruct *>(Field::GetSolverPara());
    }
    DataManagerNSStruct *NSFieldZaran::GetDataManager()
    {
        return static_cast<DataManagerNSStruct *>(Field::GetDataManager());
    }
    void NSFieldZaran::CalcResidual()
    {
        auto grid = GetGrid();
        auto node = grid->GetNode();
        auto para = GetSolverPara();
        int equ_num = para->GetEqNum();
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        int total_node_num = (ie - is + 1) * (je - js + 1) * (ke - ks + 1);
        double norm_inf = -LARGE_NUMBER;
        double norm_l2 = 0;
        int norm_inf_node = -1;
        double norm_inf_coord[3];
        for (int iEqu = 0; iEqu < equ_num; iEqu++)
        {
            norm_inf = -LARGE_NUMBER;
            norm_l2 = 0;
            auto res = GetDataManager()->GetResidual(iEqu);
            for (int k = ks; k <= ke; k++)
            {
                for (int j = js; j <= je; j++)
                {
                    for (int i = is; i <= ie; i++)
                    {
                        int idx = m_idx_proxy->GetIdx(i, j, k);
                        if (abs(res[idx]) > norm_inf)
                        {
                            norm_inf = abs(res[idx]);
                            norm_inf_node = idx;
                            for (int iDim = 0; iDim < grid->GetDim(); iDim++)
                            {
                                norm_inf_coord[iDim] = node->GetCoord(i, j, k)[iDim];
                            }
                        }
                        norm_l2 += res[idx] * res[idx];
                    }
                }
            }
            norm_l2 = sqrt(norm_l2 / total_node_num);
            m_res_info->SetInfNorm(iEqu, norm_inf);
            m_res_info->SetL2Norm(iEqu, norm_l2);
            m_res_info->SetInfNormCoord(iEqu, norm_inf_coord);
            m_res_info->SetInfNormIdx(iEqu, norm_inf_node);
        }
    }
    void NSFieldZaran::CreateSlaveField(FieldManager *field_manager)
    {
        int idx = field_manager->GetFieldNum();
        auto block = GetGrid();
        int dim = block->GetDim();
        GridFNFactoryZaran grid_factory;
        GridFN *grid = new GridFN("ZaranSlaveGrid", idx, dim);
        grid_factory.CreateGrid(block, grid, GetModelManager());
        std::vector<std::string> fn_recv_data_name, block_recv_data_name;
        fn_recv_data_name = {"primitive_0", "primitive_1", "primitive_2", "primitive_3", "primitive_4"};
        block_recv_data_name = {"primitive_0", "primitive_1", "primitive_2", "primitive_3", "primitive_4"};
        int fn_recv_num, block_recv_num;
        std::vector<int> fn_recv_node_idx_src, block_recv_node_idx_src;
        std::vector<int> fn_recv_field_idx_tgt, block_recv_field_idx_tgt;
        std::vector<int> fn_recv_node_idx_tgt, block_recv_node_idx_tgt;
        auto &fn_grid_info = grid_factory.GetFNGridInfo();
        auto &ref_node = grid_factory.GetRefNode();
        fn_recv_num = ref_node.size();
        fn_recv_node_idx_src.resize(fn_recv_num);
        fn_recv_field_idx_tgt.resize(fn_recv_num);
        fn_recv_node_idx_tgt.resize(fn_recv_num);
        int i = 0;
        for (auto &ref : ref_node)
        {
            fn_recv_node_idx_src[i] = ref.idx_local_layer;
            fn_recv_field_idx_tgt[i] = this->GetIdx();
            fn_recv_node_idx_tgt[i] = ref.idx_block;
            i++;
        }

        auto &trans_node = grid_factory.GetTransNode();
        block_recv_num = trans_node.size();
        block_recv_node_idx_src.resize(block_recv_num);
        block_recv_field_idx_tgt.resize(block_recv_num);
        block_recv_node_idx_tgt.resize(block_recv_num);
        i = 0;
        for (auto &trans : trans_node)
        {
            block_recv_node_idx_src[i] = trans.idx_block;
            block_recv_field_idx_tgt[i] = idx;
            block_recv_node_idx_tgt[i] = fn_grid_info.node[1][trans.idx_local_layer].idx;
            i++;
        }

        FieldDataCommInfo *fn_comm_info = new FieldDataCommInfo(fn_recv_num, fn_recv_data_name, fn_recv_node_idx_src.data(), fn_recv_field_idx_tgt.data(), fn_recv_node_idx_tgt.data());
        FieldDataCommInfo *block_comm_info = new FieldDataCommInfo(block_recv_num, block_recv_data_name, block_recv_node_idx_src.data(), block_recv_field_idx_tgt.data(), block_recv_node_idx_tgt.data());
        field_manager->SetFieldDataCommInfo(this->GetIdx(), block_comm_info);
        m_slave_field = new NSFieldFNFDM(grid);
        field_manager->AddField(m_slave_field, fn_comm_info);
        m_slave_field->Allocate();
    }
    //     void NSFieldZaran::CreateMasterGrid()
    //     {
    //         Log::info("AllocateMasterGrid");
    //         AllocateMasterGrid();
    //         Log::info("SetNodeCoord");
    //         SetNodeCoord();
    //         Log::info("SetBoundary");
    //         // SetBoundary();
    //         Log::info("TagCells");
    //         TagCells();
    //         Log::info("TagNodes");
    //         TagNodes();
    //         Log::info("WriteNodeTag");
    //         m_slave_grid.node.resize(11);
    //         WriteNodeTag();
    //         BuildRefNode();
    //         BuildTransNode();
    //         BuildProjectNode();
    //         BuildMiddleNode();
    //         ProjectToModel();
    //         BuildTransFace();
    //         WriteProjectNode();
    //         WriteModelSurface();
    //         WriteTransFace();
    //         WriteSlaveGrid();
    //         BuildProjectNodeNeighbor();
    //         BuildTransNodeNeighbor();
    //         exit(0);
    //     }
    //     void NSFieldZaran::AllocateMasterGrid()
    //     {
    //         int nx = GlobalData::GetInt("nx");
    //         int ny = GlobalData::GetInt("ny");
    //         int nz = GlobalData::GetInt("nz");
    //         int ghost_size = 3;
    //         GridStruct *grid = new GridStruct("ZaranMasterGrid", 1, 3);
    //         auto node = grid->GetNode();
    //         m_idx_proxy = new IdxStruct(grid);
    //         grid->Allocate(nx, ny, nz, ghost_size);
    //         SetGrid(grid);
    //         m_cell_type = new PhysicalType[(nx + 2 * ghost_size) * (ny + 2 * ghost_size) * (nz + 2 * ghost_size)];
    //         m_node_type = new PhysicalType[(nx + 2 * ghost_size) * (ny + 2 * ghost_size) * (nz + 2 * ghost_size)];
    //     }
    //     void NSFieldZaran::SetNodeCoord()
    //     {
    //         m_box.x_min = GlobalData::GetDouble("xMin");
    //         m_box.x_max = GlobalData::GetDouble("xMax");
    //         m_box.y_min = GlobalData::GetDouble("yMin");
    //         m_box.y_max = GlobalData::GetDouble("yMax");
    //         m_box.z_min = GlobalData::GetDouble("zMin");
    //         m_box.z_max = GlobalData::GetDouble("zMax");
    //         int nx = GlobalData::GetInt("nx");
    //         int ny = GlobalData::GetInt("ny");
    //         int nz = GlobalData::GetInt("nz");
    //         int ghost_size = 3;
    //         m_dx = (m_box.x_max - m_box.x_min) / nx;
    //         m_dy = (m_box.y_max - m_box.y_min) / ny;
    //         m_dz = (m_box.z_max - m_box.z_min) / nz;
    //         auto grid = GetGrid();
    //         auto node = grid->GetNode();
    //         double coord[3];
    //         for (int i = 0; i < nx + 2 * ghost_size; i++)
    //         {
    //             for (int j = 0; j < ny + 2 * ghost_size; j++)
    //             {
    //                 for (int k = 0; k < nz + 2 * ghost_size; k++)
    //                 {
    //                     coord[0] = m_box.x_min + (i - ghost_size) * m_dx;
    //                     coord[1] = m_box.y_min + (j - ghost_size) * m_dy;
    //                     coord[2] = m_box.z_min + (k - ghost_size) * m_dz;
    //                     node->SetCoord(i, j, k, coord);
    //                 }
    //             }
    //         }
    //     }
    //     void NSFieldZaran::SetBoundary()
    //     {
    //         int i_bnd, j_bnd, k_bnd;
    //         int bnd_dir[3];
    //         double norm[3];
    //         auto grid = GetGrid();
    //         auto bnd_manager = grid->GetBoundMap();
    //         int is, ie, js, je, ks, ke;
    //         grid->GetRange(is, ie, js, je, ks, ke);
    //         for (int k = ks; k < ke; ++k)
    //         {
    //             for (int j = js; j < je; ++j)
    //             {
    //                 i_bnd = is, j_bnd = j, k_bnd = k;
    //                 bnd_dir[0] = -1, bnd_dir[1] = 0, bnd_dir[2] = 0;
    //                 norm[0] = -1, norm[1] = 0, norm[2] = 0;
    //                 BoundStruct bound(i_bnd, j_bnd, k_bnd, bnd_dir, norm);
    //                 bnd_manager->AddBoundary("inlet", bound);
    //                 i_bnd = ie, j_bnd = j, k_bnd = k;
    //                 bnd_dir[0] = 1, bnd_dir[1] = 0, bnd_dir[2] = 0;
    //                 norm[0] = 1, norm[1] = 0, norm[2] = 0;
    //                 BoundStruct bound1(i_bnd, j_bnd, k_bnd, bnd_dir, norm);
    //                 bnd_manager->AddBoundary("outlet", bound1);
    //             }
    //         }
    //         for (int i = is; i < ie; ++i)
    //         {
    //             for (int j = js; j < je; ++j)
    //             {
    //                 k_bnd = ks, i_bnd = i, j_bnd = j;
    //                 bnd_dir[0] = 0, bnd_dir[1] = 0, bnd_dir[2] = -1;
    //                 norm[0] = 0, norm[1] = 0, norm[2] = -1;
    //                 BoundStruct bound2(i_bnd, j_bnd, k_bnd, bnd_dir, norm);
    //                 bnd_manager->AddBoundary("outlet", bound2);
    //                 k_bnd = ke, i_bnd = i, j_bnd = j;
    //                 bnd_dir[0] = 0, bnd_dir[1] = 0, bnd_dir[2] = 1;
    //                 norm[0] = 0, norm[1] = 0, norm[2] = 1;
    //                 BoundStruct bound3(i_bnd, j_bnd, k_bnd, bnd_dir, norm);
    //                 bnd_manager->AddBoundary("outlet", bound3);
    //             }
    //         }
    //         for (int i = is; i < ie; ++i)
    //         {
    //             for (int k = ks; k < ke; ++k)
    //             {
    //                 j_bnd = js, i_bnd = i, k_bnd = k;
    //                 bnd_dir[0] = 0, bnd_dir[1] = -1, bnd_dir[2] = 0;
    //                 norm[0] = 0, norm[1] = -1, norm[2] = 0;
    //                 BoundStruct bound4(i_bnd, j_bnd, k_bnd, bnd_dir, norm);
    //                 bnd_manager->AddBoundary("outlet", bound4);
    //                 j_bnd = je, i_bnd = i, k_bnd = k;
    //                 bnd_dir[0] = 0, bnd_dir[1] = 1, bnd_dir[2] = 0;
    //                 norm[0] = 0, norm[1] = 1, norm[2] = 0;
    //                 BoundStruct bound5(i_bnd, j_bnd, k_bnd, bnd_dir, norm);
    //                 bnd_manager->AddBoundary("outlet", bound5);
    //             }
    //         }
    //     }
    //     void NSFieldZaran::TagCells()
    //     {
    //         auto grid = GetGrid();
    //         int ni, nj, nk;
    //         ni = grid->GetNi();
    //         nj = grid->GetNj();
    //         nk = grid->GetNk();
    //         int ghost_size = grid->GetGhostLevel();
    //         Log::info("tag cells");
    //         // tag the solid cells by distance
    //         double tol = 0.5 * sqrt(m_dx * m_dx + m_dy * m_dy + m_dz * m_dz);
    //         Log::info("tag cells");
    //         for (int i = 0; i < ni - 1; i++)
    //         {
    //             for (int j = 0; j < nj - 1; j++)
    //             {
    //                 for (int k = 0; k < nk - 1; k++)
    //                 {
    //                     int idx = m_idx_proxy->GetIdx(i, j, k);
    //                     m_cell_type[idx] = PhysicalType::Solid;
    //                 }
    //             }
    //         }
    //         Log::info("tag cells");
    //         auto &box = m_model_manager->GetBox();
    //         Log::info("Model box: x_min={}, x_max={}", box.x_min, box.x_max);
    //         Log::info("Model box: y_min={}, y_max={}", box.y_min, box.y_max);
    //         Log::info("Model box: z_min={}, z_max={}", box.z_min, box.z_max);

    // #pragma omp parallel for
    //         for (int iCell = 0; iCell < ni * nj * nk; iCell++)
    //         {
    //             int i, j, k;
    //             m_idx_proxy->GetIdxStruct(iCell, i, j, k);
    //             double cell_center[3];
    //             cell_center[0] = m_box.x_min + (i - ghost_size + 0.5) * m_dx;
    //             cell_center[1] = m_box.y_min + (j - ghost_size + 0.5) * m_dy;
    //             cell_center[2] = m_box.z_min + (k - ghost_size + 0.5) * m_dz;
    //             if (cell_center[0] < box.x_min - tol || cell_center[0] > box.x_max + tol || cell_center[1] < box.y_min - tol || cell_center[1] > box.y_max + tol || cell_center[2] < box.z_min - tol || cell_center[2] > box.z_max + tol)
    //             {
    //                 m_cell_type[iCell] = PhysicalType::Fluid;
    //             }
    //             else
    //             {
    //                 double dist = m_model_manager->GetClosestDistance(cell_center);
    //                 if (dist < tol)
    //                 {
    //                     m_cell_type[iCell] = PhysicalType::Solid;
    //                 }
    //                 else
    //                 {
    //                     m_cell_type[iCell] = PhysicalType::Unset;
    //                 }
    //             }
    //         }
    //         Log::info("tag cells");
    //         // tag the fluid cells from the fist cell
    //         m_cell_type[m_idx_proxy->GetIdx(0, 0, 0)] = PhysicalType::Fluid;
    //         for (int i = 0; i < ni - 2; i++)
    //         {
    //             for (int j = 0; j < nj - 2; j++)
    //             {
    //                 for (int k = 0; k < nk - 2; k++)
    //                 {
    //                     int idx = m_idx_proxy->GetIdx(i, j, k);
    //                     if (m_cell_type[idx] == PhysicalType::Fluid)
    //                     {
    //                         if (m_cell_type[m_idx_proxy->GetIdx(i + 1, j, k)] == PhysicalType::Unset)
    //                         {
    //                             m_cell_type[m_idx_proxy->GetIdx(i + 1, j, k)] = PhysicalType::Fluid;
    //                         }
    //                         if (m_cell_type[m_idx_proxy->GetIdx(i, j + 1, k)] == PhysicalType::Unset)
    //                         {
    //                             m_cell_type[m_idx_proxy->GetIdx(i, j + 1, k)] = PhysicalType::Fluid;
    //                         }
    //                         if (m_cell_type[m_idx_proxy->GetIdx(i, j, k + 1)] == PhysicalType::Unset)
    //                         {
    //                             m_cell_type[m_idx_proxy->GetIdx(i, j, k + 1)] = PhysicalType::Fluid;
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //         Log::info("tag cells");

    //         // tag the solid cells after the fluid cells
    //         for (int i = 0; i < ni - 1; i++)
    //         {
    //             for (int j = 0; j < nj - 1; j++)
    //             {
    //                 for (int k = 0; k < nk - 1; k++)
    //                 {
    //                     int idx = m_idx_proxy->GetIdx(i, j, k);
    //                     if (m_cell_type[idx] == PhysicalType::Unset)
    //                     {
    //                         m_cell_type[idx] = PhysicalType::Solid;
    //                     }
    //                 }
    //             }
    //         }
    //         Log::info("tag cells");

    //         // retag the bad fluid cells as solid
    //         // int new_solid_num = 1;
    //         // while (new_solid_num > 0)
    //         // {
    //         //     new_solid_num = 0;
    //         //     for (int i = 0; i < ni - 1; i++)
    //         //     {
    //         //         for (int j = 0; j < nj - 1; j++)
    //         //         {
    //         //             for (int k = 0; k < -1; k++)
    //         //             {
    //         //                 int idx = m_idx_proxy->GetIdx(i, j, k);
    //         //                 if (m_cell_type[idx] != PhysicalType::Fluid)
    //         //                     continue;
    //         //                 if (m_cell_type[m_idx_proxy->GetIdx(i + 1, j, k)] == PhysicalType::Solid && m_cell_type[m_idx_proxy->GetIdx(i - 1, j, k)] == PhysicalType::Solid)
    //         //                 {
    //         //                     m_cell_type[idx] = PhysicalType::Solid;
    //         //                     new_solid_num++;
    //         //                 }
    //         //                 else if (m_cell_type[m_idx_proxy->GetIdx(i, j + 1, k)] == PhysicalType::Solid && m_cell_type[m_idx_proxy->GetIdx(i, j - 1, k)] == PhysicalType::Solid)
    //         //                 {
    //         //                     m_cell_type[idx] = PhysicalType::Solid;
    //         //                     new_solid_num++;
    //         //                 }
    //         //                 else if (m_cell_type[m_idx_proxy->GetIdx(i, j, k + 1)] == PhysicalType::Solid && m_cell_type[m_idx_proxy->GetIdx(i, j, k - 1)] == PhysicalType::Solid)
    //         //                 {
    //         //                     m_cell_type[idx] = PhysicalType::Solid;
    //         //                     new_solid_num++;
    //         //                 }
    //         //             }
    //         //         }
    //         //     }
    //         // }
    //         Log::info("tag cells");

    //         // tag the fluid-solid cells
    //         for (int k = 1; k < nk - 1; k++)
    //         {
    //             for (int j = 1; j < nj - 1; j++)
    //             {
    //                 for (int i = 1; i < ni - 1; i++)
    //                 {
    //                     int idx = m_idx_proxy->GetIdx(i, j, k);
    //                     if (m_cell_type[idx] != PhysicalType::Fluid)
    //                         continue;
    //                     if (m_cell_type[m_idx_proxy->GetIdx(i + 1, j, k)] == PhysicalType::Solid)
    //                     {
    //                         m_cell_type[idx] = PhysicalType::FluidSolid;
    //                     }
    //                     if (m_cell_type[m_idx_proxy->GetIdx(i - 1, j, k)] == PhysicalType::Solid)
    //                     {
    //                         m_cell_type[idx] = PhysicalType::FluidSolid;
    //                     }
    //                     if (m_cell_type[m_idx_proxy->GetIdx(i, j + 1, k)] == PhysicalType::Solid)
    //                     {
    //                         m_cell_type[idx] = PhysicalType::FluidSolid;
    //                     }
    //                     if (m_cell_type[m_idx_proxy->GetIdx(i, j - 1, k)] == PhysicalType::Solid)
    //                     {
    //                         m_cell_type[idx] = PhysicalType::FluidSolid;
    //                     }
    //                     if (m_cell_type[m_idx_proxy->GetIdx(i, j, k + 1)] == PhysicalType::Solid)
    //                     {
    //                         m_cell_type[idx] = PhysicalType::FluidSolid;
    //                     }
    //                     if (m_cell_type[m_idx_proxy->GetIdx(i, j, k - 1)] == PhysicalType::Solid)
    //                     {
    //                         m_cell_type[idx] = PhysicalType::FluidSolid;
    //                     }
    //                 }
    //             }
    //         }
    //         Log::info("tag cells");
    //     }
    //     void NSFieldZaran::TagNodes()
    //     {
    //         auto grid = GetGrid();
    //         auto node = grid->GetNode();
    //         int ni, nj, nk;
    //         ni = grid->GetNi();
    //         nj = grid->GetNj();
    //         nk = grid->GetNk();
    //         auto IsValidNode = [&](int i, int j, int k) -> bool
    //         {
    //             return i >= 0 && i < ni && j >= 0 && j < nj && k >= 0 && k < nk;
    //         };
    //         auto IsValidCell = [&](int i, int j, int k) -> bool
    //         {
    //             return i >= 0 && i < ni - 1 && j >= 0 && j < nj - 1 && k >= 0 && k < nk - 1;
    //         };
    //         int is, ie, js, je, ks, ke;
    //         grid->GetRange(is, ie, js, je, ks, ke);
    //         // initialize the node type
    //         Log::info("initialize the node type");
    //         for (int idx = 0; idx < ni * nj * nk; idx++)
    //         {
    //             m_node_type[idx] = PhysicalType::Unset;
    //         }
    //         Log::info("initialize the node type");
    //         std::set<int> trans_node_set;
    //         std::set<TransFace> trans_face_set;
    //         // find the fluid-solid nodes
    //         for (int i = 0; i < ni; i++)
    //         {
    //             for (int j = 0; j < nj; j++)
    //             {
    //                 for (int k = 0; k < nk; k++)
    //                 {
    //                     int idx = m_idx_proxy->GetIdx(i, j, k);
    //                     if (IsValidCell(i - 1, j, k))
    //                     {
    //                         int idx_left = m_idx_proxy->GetIdx(i - 1, j, k);
    //                         if (m_cell_type[idx_left] == PhysicalType::FluidSolid && m_cell_type[idx] == PhysicalType::Fluid || m_cell_type[idx_left] == PhysicalType::Fluid && m_cell_type[idx] == PhysicalType::FluidSolid)
    //                         {
    //                             trans_face_set.insert(TransFace{{m_idx_proxy->GetIdx(i, j, k), m_idx_proxy->GetIdx(i, j + 1, k), m_idx_proxy->GetIdx(i, j + 1, k + 1), m_idx_proxy->GetIdx(i, j, k + 1)}});
    //                             m_node_type[idx] = PhysicalType::FluidSolid;
    //                             m_node_type[m_idx_proxy->GetIdx(i, j + 1, k)] = PhysicalType::FluidSolid;
    //                             m_node_type[m_idx_proxy->GetIdx(i, j, k + 1)] = PhysicalType::FluidSolid;
    //                             m_node_type[m_idx_proxy->GetIdx(i, j + 1, k + 1)] = PhysicalType::FluidSolid;
    //                         }
    //                     }
    //                     if (IsValidCell(i, j - 1, k))
    //                     {
    //                         int idx_down = m_idx_proxy->GetIdx(i, j - 1, k);
    //                         if (m_cell_type[idx_down] == PhysicalType::FluidSolid && m_cell_type[idx] == PhysicalType::Fluid || m_cell_type[idx_down] == PhysicalType::Fluid && m_cell_type[idx] == PhysicalType::FluidSolid)
    //                         {
    //                             trans_face_set.insert(TransFace{{m_idx_proxy->GetIdx(i, j, k), m_idx_proxy->GetIdx(i + 1, j, k), m_idx_proxy->GetIdx(i + 1, j, k + 1), m_idx_proxy->GetIdx(i, j, k + 1)}});
    //                             m_node_type[idx] = PhysicalType::FluidSolid;
    //                             m_node_type[m_idx_proxy->GetIdx(i + 1, j, k)] = PhysicalType::FluidSolid;
    //                             m_node_type[m_idx_proxy->GetIdx(i, j, k + 1)] = PhysicalType::FluidSolid;
    //                             m_node_type[m_idx_proxy->GetIdx(i + 1, j, k + 1)] = PhysicalType::FluidSolid;
    //                         }
    //                     }
    //                     if (IsValidCell(i, j, k - 1))
    //                     {
    //                         int idx_back = m_idx_proxy->GetIdx(i, j, k - 1);
    //                         if (m_cell_type[idx_back] == PhysicalType::FluidSolid && m_cell_type[idx] == PhysicalType::Fluid || m_cell_type[idx_back] == PhysicalType::Fluid && m_cell_type[idx] == PhysicalType::FluidSolid)
    //                         {
    //                             trans_face_set.insert(TransFace{{m_idx_proxy->GetIdx(i, j, k), m_idx_proxy->GetIdx(i + 1, j, k), m_idx_proxy->GetIdx(i + 1, j + 1, k), m_idx_proxy->GetIdx(i, j + 1, k)}});
    //                             m_node_type[idx] = PhysicalType::FluidSolid;
    //                             m_node_type[m_idx_proxy->GetIdx(i + 1, j, k)] = PhysicalType::FluidSolid;
    //                             m_node_type[m_idx_proxy->GetIdx(i, j + 1, k)] = PhysicalType::FluidSolid;
    //                             m_node_type[m_idx_proxy->GetIdx(i + 1, j + 1, k)] = PhysicalType::FluidSolid;
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //         Log::info("initialize the node type");
    //         m_trans_face.resize(trans_face_set.size());
    //         int idx = 0;
    //         for (auto &trans : trans_face_set)
    //         {
    //             m_trans_face[idx++].idx_master = trans.idx_master;
    //         }

    //         m_node_type[m_idx_proxy->GetIdx(0, 0, 0)] = PhysicalType::Fluid;
    //         // find the fluid nodes
    //         Log::info("find the fluid nodes");
    //         for (int k = 0; k < nk; k++)
    //         {
    //             for (int j = 0; j < nj; j++)
    //             {
    //                 for (int i = 0; i < ni; i++)
    //                 {
    //                     if (!IsValidNode(i, j, k))
    //                     {
    //                         Log::info("Invalid node: i={}, j={}, k={}, idx={}, ni={}, nj={}, nk={}", i, j, k, m_idx_proxy->GetIdx(i, j, k), ni, nj, nk);
    //                         exit(0);
    //                     }
    //                     int idx = m_idx_proxy->GetIdx(i, j, k);
    //                     if (m_node_type[idx] != PhysicalType::Fluid)
    //                         continue;
    //                     if (IsValidNode(i + 1, j, k))
    //                     {
    //                         if (m_node_type[m_idx_proxy->GetIdx(i + 1, j, k)] == PhysicalType::Unset)
    //                             m_node_type[m_idx_proxy->GetIdx(i + 1, j, k)] = PhysicalType::Fluid;
    //                     }
    //                     if (IsValidNode(i, j + 1, k))
    //                     {
    //                         if (m_node_type[m_idx_proxy->GetIdx(i, j + 1, k)] == PhysicalType::Unset)
    //                             m_node_type[m_idx_proxy->GetIdx(i, j + 1, k)] = PhysicalType::Fluid;
    //                     }
    //                     if (IsValidNode(i, j, k + 1))
    //                     {
    //                         if (m_node_type[m_idx_proxy->GetIdx(i, j, k + 1)] == PhysicalType::Unset)
    //                             m_node_type[m_idx_proxy->GetIdx(i, j, k + 1)] = PhysicalType::Fluid;
    //                     }
    //                 }
    //             }
    //         }
    //         // find the solid nodes
    //         Log::info("find the solid nodes");
    //         for (int iNode = 0; iNode < ni * nj * nk; iNode++)
    //         {
    //             if (m_node_type[iNode] == PhysicalType::Unset)
    //             {
    //                 m_node_type[iNode] = PhysicalType::Solid;
    //             }
    //         }
    //     }
    //     void NSFieldZaran::WriteNodeTag()
    //     {
    //         std::ofstream out("node_tag.dat");
    //         auto grid = GetGrid();
    //         auto node = grid->GetNode();
    //         int is, ie, js, je, ks, ke;
    //         grid->GetRange(is, ie, js, je, ks, ke);
    //         out << "variables = x, y, z, type" << std::endl;
    //         out << "Zone T=test\n I=" << ie - is << " J=" << je - js << " K=" << ke - ks << " F=POINT" << std::endl;
    //         for (int k = ks; k < ke; k++)
    //         {
    //             for (int j = js; j < je; j++)
    //             {
    //                 for (int i = is; i < ie; i++)
    //                 {
    //                     int idx = m_idx_proxy->GetIdx(i, j, k);
    //                     out << node->GetCoord(i, j, k)[0] << " " << node->GetCoord(i, j, k)[1] << " " << node->GetCoord(i, j, k)[2] << " " << (int)m_node_type[idx] << std::endl;
    //                 }
    //             }
    //         }
    //     }
    //     void NSFieldZaran::WriteProjectNode()
    //     {
    //         // std::ofstream out("project_node.dat");
    //         // out << "variables = x, y, z" << std::endl;
    //         // out << "Zone T=test" << std::endl;
    //         // for (auto &patch : m_bound_patch)
    //         // {
    //         //     out << patch.coord[0] << " " << patch.coord[1] << " " << patch.coord[2] << std::endl;
    //         // }
    //     }
    //     void NSFieldZaran::WriteModelSurface()
    //     {
    //         // std::ofstream out("model_surface.dat");
    //         // out << "variables = x, y, z" << std::endl;
    //         // out << "Zone T=test" << std::endl;
    //         // out << " N=" << m_bound_patch.size() << " E=" << m_trans_face.size() << " ZONETYPE=FEQuadrilateral" << std::endl;
    //         // out << "DATAPACKING=POINT" << std::endl;
    //         // for (auto &patch : m_bound_patch)
    //         // {
    //         //     out << patch.coord[0] << " " << patch.coord[1] << " " << patch.coord[2] << std::endl;
    //         // }
    //         // for (int iFace = 0; iFace < m_trans_face.size(); iFace++)
    //         // {
    //         //     for (int iNode = 0; iNode < m_trans_face[iFace].idx_slave.size(); iNode++)
    //         //     {
    //         //         out << m_trans_face[iFace].idx_slave[iNode] + 1 << " ";
    //         //     }
    //         //     out << std::endl;
    //         // }
    //         // out.close();
    //     }
    //     void NSFieldZaran::WriteTransFace()
    //     {
    //         auto grid = GetGrid();
    //         auto node = grid->GetNode();
    //         int ni, nj, nk;
    //         ni = grid->GetNi();
    //         nj = grid->GetNj();
    //         nk = grid->GetNk();
    //         int node_num = ni * nj * nk;
    //         std::ofstream out("trans_face.dat");
    //         out << "variables = x, y, z" << std::endl;
    //         out << "Zone T=test" << std::endl;
    //         out << " N=" << node_num << " E=" << m_trans_face.size() << " ZONETYPE=FEQuadrilateral" << std::endl;
    //         out << "DATAPACKING=POINT" << std::endl;
    //         for (int k = 0; k < nk; k++)
    //         {
    //             for (int j = 0; j < nj; j++)
    //             {
    //                 for (int i = 0; i < ni; i++)
    //                 {
    //                     out << node->GetCoord(i, j, k)[0] << " " << node->GetCoord(i, j, k)[1] << " " << node->GetCoord(i, j, k)[2] << std::endl;
    //                 }
    //             }
    //         }
    //         for (int iFace = 0; iFace < m_trans_face.size(); iFace++)
    //         {
    //             for (int iNode = 0; iNode < m_trans_face[iFace].idx_master.size(); iNode++)
    //             {
    //                 out << m_trans_face[iFace].idx_master[iNode] + 1 << " ";
    //             }
    //             out << std::endl;
    //         }
    //         out.close();
    //     }
    //     void NSFieldZaran::WriteSlaveGrid()
    //     {
    //         std::ofstream out("slave_grid.dat");
    //         out << "variables = x, y, z" << std::endl;
    //         out << "Zone T=test" << std::endl;
    //         int node_num = 0;
    //         for (int iLayer = 0; iLayer < m_slave_grid.node.size(); iLayer++)
    //         {
    //             node_num += m_slave_grid.node[iLayer].size();
    //         }
    //         int cell_num = 0;
    //         cell_num = m_trans_face.size() * (m_slave_grid.node.size() - 2);
    //         out << " N=" << node_num << " E=" << cell_num << " F=FEPOINT, ET=Brick" << std::endl;
    //         for (int iLayer = 0; iLayer < m_slave_grid.node.size(); iLayer++)
    //         {
    //             for (int iNode = 0; iNode < m_slave_grid.node[iLayer].size(); iNode++)
    //             {
    //                 out << m_slave_grid.node[iLayer][iNode].coord[0] << " " << m_slave_grid.node[iLayer][iNode].coord[1] << " " << m_slave_grid.node[iLayer][iNode].coord[2] << std::endl;
    //             }
    //         }
    //         for (int iFace = 0; iFace < m_trans_face.size(); iFace++)
    //         {
    //             auto trans_face = m_trans_face[iFace].idx_master;
    //             for (int iNode = 0; iNode < trans_face.size(); iNode++)
    //             {
    //                 for (auto &slave : m_slave_connect_to_master)
    //                 {
    //                     if (slave.idx_master == trans_face[iNode])
    //                     {
    //                         trans_face[iNode] = slave.idx_slave_layer;
    //                         break;
    //                     }
    //                 }
    //             }

    //             int start_idx = 0;
    //             int next_idx = 0;
    //             for (int iLayer = 1; iLayer < m_slave_grid.node.size() - 1; iLayer++)
    //             {
    //                 start_idx += m_slave_grid.node[iLayer - 1].size();
    //                 next_idx = start_idx + m_slave_grid.node[iLayer].size();
    //                 for (int iNode = 0; iNode < trans_face.size(); iNode++)
    //                 {
    //                     out << start_idx + trans_face[iNode] + 1 << " ";
    //                 }
    //                 for (int iNode = 0; iNode < trans_face.size(); iNode++)
    //                 {
    //                     out << next_idx + trans_face[iNode] + 1 << " ";
    //                 }
    //                 out << std::endl;
    //             }
    //         }
    //         out.close();
    //     }
    //     void NSFieldZaran::BuildRefNode()
    //     {
    //         auto grid = GetGrid();
    //         auto node = grid->GetNode();
    //         int is, ie, js, je, ks, ke;
    //         grid->GetRange(is, ie, js, je, ks, ke);
    //         std::set<int> ref_node_set;
    //         for (int k = ks; k <= ke; k++)
    //         {
    //             for (int j = js; j <= je; j++)
    //             {
    //                 for (int i = is; i <= ie; i++)
    //                 {
    //                     int idx = m_idx_proxy->GetIdx(i, j, k);
    //                     if (m_node_type[idx] == PhysicalType::FluidSolid)
    //                     {
    //                         if (m_node_type[m_idx_proxy->GetIdx(i + 1, j, k)] == PhysicalType::Fluid)
    //                         {
    //                             ref_node_set.insert(m_idx_proxy->GetIdx(i + 1, j, k));
    //                         }
    //                         if (m_node_type[m_idx_proxy->GetIdx(i - 1, j, k)] == PhysicalType::Fluid)
    //                         {
    //                             ref_node_set.insert(m_idx_proxy->GetIdx(i - 1, j, k));
    //                         }
    //                         if (m_node_type[m_idx_proxy->GetIdx(i, j + 1, k)] == PhysicalType::Fluid)
    //                         {
    //                             ref_node_set.insert(m_idx_proxy->GetIdx(i, j + 1, k));
    //                         }
    //                         if (m_node_type[m_idx_proxy->GetIdx(i, j - 1, k)] == PhysicalType::Fluid)
    //                         {
    //                             ref_node_set.insert(m_idx_proxy->GetIdx(i, j - 1, k));
    //                         }
    //                         if (m_node_type[m_idx_proxy->GetIdx(i, j, k + 1)] == PhysicalType::Fluid)
    //                         {
    //                             ref_node_set.insert(m_idx_proxy->GetIdx(i, j, k + 1));
    //                         }
    //                         if (m_node_type[m_idx_proxy->GetIdx(i, j, k - 1)] == PhysicalType::Fluid)
    //                         {
    //                             ref_node_set.insert(m_idx_proxy->GetIdx(i, j, k - 1));
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //         m_slave_grid.node[0].resize(ref_node_set.size());
    //         int idx = 0;
    //         for (auto &ref_node : ref_node_set)
    //         {
    //             int i, j, k;
    //             m_idx_proxy->GetIdxStruct(ref_node, i, j, k);
    //             m_slave_grid.node[0][idx].coord[0] = node->GetCoord(i, j, k)[0];
    //             m_slave_grid.node[0][idx].coord[1] = node->GetCoord(i, j, k)[1];
    //             m_slave_grid.node[0][idx].coord[2] = node->GetCoord(i, j, k)[2];
    //             m_slave_grid.node[0][idx].idx = idx;
    //             m_slave_connect_to_master.insert(SlaveConnectToMaster{ref_node, 0, idx});
    //             idx++;
    //         }
    //     }
    //     void NSFieldZaran::BuildTransNode()
    //     {
    //         auto grid = GetGrid();
    //         auto node = grid->GetNode();
    //         int is, ie, js, je, ks, ke;
    //         grid->GetRange(is, ie, js, je, ks, ke);
    //         std::set<int> trans_node_set;
    //         for (int k = ks; k <= ke; k++)
    //         {
    //             for (int j = js; j <= je; j++)
    //             {
    //                 for (int i = is; i <= ie; i++)
    //                 {
    //                     int idx = m_idx_proxy->GetIdx(i, j, k);
    //                     if (m_node_type[idx] == PhysicalType::FluidSolid)
    //                     {
    //                         trans_node_set.insert(idx);
    //                     }
    //                 }
    //             }
    //         }
    //         m_slave_grid.node[1].resize(trans_node_set.size());
    //         int idx = 0;
    //         for (auto &trans : trans_node_set)
    //         {
    //             int i, j, k;
    //             m_idx_proxy->GetIdxStruct(trans, i, j, k);
    //             m_slave_grid.node[1][idx].coord[0] = node->GetCoord(i, j, k)[0];
    //             m_slave_grid.node[1][idx].coord[1] = node->GetCoord(i, j, k)[1];
    //             m_slave_grid.node[1][idx].coord[2] = node->GetCoord(i, j, k)[2];
    //             m_slave_grid.node[1][idx].idx = idx + m_slave_grid.node[0].size();
    //             m_slave_connect_to_master.insert(SlaveConnectToMaster{trans, 1, idx});
    //             idx++;
    //         }
    //     }
    //     void NSFieldZaran::BuildProjectNode()
    //     {
    //         m_slave_grid.node[10].resize(m_slave_grid.node[1].size());
    //         for (int iNode = 0; iNode < m_slave_grid.node[1].size(); iNode++)
    //         {
    //             int i, j, k;
    //             auto trans_coord = m_slave_grid.node[1][iNode].coord;
    //             double proj_coord[3];
    //             m_model_manager->GetClosestPoint(trans_coord, proj_coord);
    //             m_slave_grid.node[10][iNode].coord[0] = proj_coord[0];
    //             m_slave_grid.node[10][iNode].coord[1] = proj_coord[1];
    //             m_slave_grid.node[10][iNode].coord[2] = proj_coord[2];
    //             m_slave_grid.node[10][iNode].idx = iNode + m_slave_grid.node[0].size() + m_slave_grid.node[1].size() * 9;
    //         }
    //     }
    //     void NSFieldZaran::BuildMiddleNode()
    //     {
    //         for (int iLayer = 2; iLayer < 10; iLayer++)
    //         {
    //             m_slave_grid.node[iLayer].resize(m_slave_grid.node[1].size());
    //         }
    //         for (int iNode = 0; iNode < m_slave_grid.node[1].size(); iNode++)
    //         {
    //             int i, j, k;
    //             auto trans_coord = m_slave_grid.node[1][iNode].coord;
    //             auto proj_coord = m_slave_grid.node[10][iNode].coord;
    //             for (int iLayer = 2; iLayer < 10; iLayer++)
    //             {
    //                 m_slave_grid.node[iLayer][iNode].coord[0] = (iLayer - 1) / 9.0 * (proj_coord[0] - trans_coord[0]) + trans_coord[0];
    //                 m_slave_grid.node[iLayer][iNode].coord[1] = (iLayer - 1) / 9.0 * (proj_coord[1] - trans_coord[1]) + trans_coord[1];
    //                 m_slave_grid.node[iLayer][iNode].coord[2] = (iLayer - 1) / 9.0 * (proj_coord[2] - trans_coord[2]) + trans_coord[2];
    //                 m_slave_grid.node[iLayer][iNode].idx = iNode + m_slave_grid.node[0].size() + m_slave_grid.node[1].size() * (iLayer - 1);
    //             }
    //         }
    //     }
    //     void NSFieldZaran::BuildProjectNodeNeighbor()
    //     {
    //         // direct neighbor
    //         std::vector<std::set<int>> node_neighbor_origin;
    //         node_neighbor_origin.resize(m_slave_grid.node[1].size());
    //         for (int iFace = 0; iFace < m_trans_face.size(); iFace++)
    //         {
    //             auto trans_face = m_trans_face[iFace].idx_master;
    //             for (int iNode = 0; iNode < trans_face.size(); iNode++)
    //             {
    //                 for (auto &slave : m_slave_connect_to_master)
    //                 {
    //                     if (slave.idx_master == trans_face[iNode])
    //                     {
    //                         trans_face[iNode] = slave.idx_slave_layer;
    //                         break;
    //                     }
    //                 }
    //             }
    //             int prior_node, next_node;
    //             for (int iNode = 0; iNode < trans_face.size(); iNode++)
    //             {
    //                 int current_node = trans_face[iNode];
    //                 prior_node = trans_face[(iNode + trans_face.size() - 1) % trans_face.size()];
    //                 next_node = trans_face[(iNode + 1) % trans_face.size()];
    //                 node_neighbor_origin[current_node].insert(prior_node);
    //                 node_neighbor_origin[current_node].insert(next_node);
    //             }
    //         }
    //         // extend neighbor, if the neighbor node is less than 3, extend the neighbor node
    //         std::vector<std::set<int>> node_neighbor_extend;
    //         node_neighbor_extend.resize(m_slave_grid.node[1].size());
    //         for (int iNode = 0; iNode < node_neighbor_origin.size(); iNode++)
    //         {
    //             node_neighbor_extend[iNode] = node_neighbor_origin[iNode];

    //             if (node_neighbor_extend[iNode].size() <= 3)
    //             {
    //                 for (auto &neighbor : node_neighbor_origin[iNode])
    //                 {
    //                     for (auto &neighbor_neighbor : node_neighbor_origin[neighbor])
    //                     {
    //                         node_neighbor_extend[iNode].insert(neighbor_neighbor);
    //                     }
    //                 }
    //             }
    //         }
    //         // build the neighbor node
    //         for (int iNode = 0; iNode < m_slave_grid.node[1].size(); iNode++)
    //         {
    //             for (int iLayer = 2; iLayer < 10; iLayer++)
    //             {
    //                 m_slave_grid.node[iLayer][iNode].neighbor_node.resize(node_neighbor_extend[iNode].size() + 2);
    //                 m_slave_grid.node[iLayer][iNode].neighbor_node[0] = m_slave_grid.node[iLayer - 1][iNode].idx;
    //                 m_slave_grid.node[iLayer][iNode].neighbor_node[1] = m_slave_grid.node[iLayer + 1][iNode].idx;
    //                 for (auto &neighbor : node_neighbor_extend[iNode])
    //                 {
    //                     m_slave_grid.node[iLayer][iNode].neighbor_node.push_back(m_slave_grid.node[iLayer][neighbor].idx);
    //                 }
    //             }
    //             m_slave_grid.node[10][iNode].neighbor_node.resize(node_neighbor_extend[iNode].size() + 2);
    //             m_slave_grid.node[10][iNode].neighbor_node[0] = m_slave_grid.node[9][iNode].idx;
    //             m_slave_grid.node[10][iNode].neighbor_node[1] = -1;
    //             for (auto &neighbor : node_neighbor_extend[iNode])
    //             {
    //                 m_slave_grid.node[10][iNode].neighbor_node.push_back(m_slave_grid.node[10][neighbor].idx);
    //             }
    //         }
    //     }
    //     void NSFieldZaran::BuildTransNodeNeighbor()
    //     {
    //         for (int iNode = 0; iNode < m_slave_grid.node[1].size(); iNode++)
    //         {
    //             m_slave_grid.node[1][iNode].neighbor_node.resize(6);
    //         }
    //         for (auto &nodes : m_slave_connect_to_master)
    //         {
    //             if (m_node_type[nodes.idx_master] == PhysicalType::FluidSolid)
    //             {
    //                 int i, j, k;
    //                 m_idx_proxy->GetIdxStruct(nodes.idx_master, i, j, k);
    //                 m_slave_grid.node[1][nodes.idx_slave_layer].neighbor_node[0] = m_idx_proxy->GetIdx(i - 1, j, k);
    //                 m_slave_grid.node[1][nodes.idx_slave_layer].neighbor_node[1] = m_idx_proxy->GetIdx(i + 1, j, k);
    //                 m_slave_grid.node[1][nodes.idx_slave_layer].neighbor_node[2] = m_idx_proxy->GetIdx(i, j - 1, k);
    //                 m_slave_grid.node[1][nodes.idx_slave_layer].neighbor_node[3] = m_idx_proxy->GetIdx(i, j + 1, k);
    //                 m_slave_grid.node[1][nodes.idx_slave_layer].neighbor_node[4] = m_idx_proxy->GetIdx(i, j, k - 1);
    //                 m_slave_grid.node[1][nodes.idx_slave_layer].neighbor_node[5] = m_idx_proxy->GetIdx(i, j, k + 1);
    //             }
    //         }
    //         for (int iNode = 0; iNode < m_slave_grid.node[1].size(); iNode++)
    //         {
    //             int next_layer_node = m_slave_grid.node[2][iNode].idx;
    //             for (int iNeighbor = 0; iNeighbor < 6; iNeighbor++)
    //             {
    //                 int idx_master = m_slave_grid.node[1][iNode].neighbor_node[iNeighbor];
    //                 if (m_node_type[idx_master] == PhysicalType::FluidSolid || m_node_type[idx_master] == PhysicalType::Fluid)
    //                 {
    //                     for (auto &nodes : m_slave_connect_to_master)
    //                     {
    //                         if (nodes.idx_master == idx_master)
    //                         {
    //                             m_slave_grid.node[1][iNode].neighbor_node[iNeighbor] = m_slave_grid.node[nodes.idx_slave_layer][nodes.idx_slave_proj].idx;
    //                             break;
    //                         }
    //                     }
    //                 }
    //                 else if (m_node_type[idx_master] == PhysicalType::Solid)
    //                 {
    //                     m_slave_grid.node[1][iNode].neighbor_node[iNeighbor] = next_layer_node;
    //                 }
    //                 else
    //                 {
    //                     Log::error("Invalid node type");
    //                     exit(0);
    //                 }
    //             }
    //         }
    //     }
    //     void NSFieldZaran::BuildCell()
    //     {
    //         m_slave_grid.cell.resize(m_trans_face.size() * 8);
    //         for (int iFace = 0; iFace < m_trans_face.size(); iFace++)
    //         {
    //             auto trans_face = m_trans_face[iFace].idx_master;
    //             for (int iNode = 0; iNode < trans_face.size(); iNode++)
    //             {
    //                 for (auto &slave : m_slave_connect_to_master)
    //                 {
    //                     if (slave.idx_master == trans_face[iNode])
    //                     {
    //                         trans_face[iNode] = slave.idx_slave_layer;
    //                         break;
    //                     }
    //                 }
    //             }
    //             for (int iLevel = 1; iLevel < 10; iLevel++)
    //             {
    //                 int iCell = iFace * 8 + iLevel - 1;
    //                 m_slave_grid.cell[iCell].resize(8);
    //                 for (int iNode = 0; iNode < trans_face.size(); iNode++)
    //                 {
    //                     m_slave_grid.cell[iCell][iNode] = m_slave_grid.node[iLevel][trans_face[iNode]].idx;
    //                 }
    //                 for (int iNode = 0; iNode < trans_face.size(); iNode++)
    //                 {
    //                     m_slave_grid.cell[iCell][iNode + 4] = m_slave_grid.node[iLevel + 1][trans_face[iNode]].idx;
    //                 }
    //             }
    //         }
    //     }
    //     void NSFieldZaran::ProjectToModel()
    //     {
    //         // auto grid = GetGrid();
    //         // auto node = grid->GetNode();
    //         // for (int iNode = 0; iNode < m_bound_patch.size(); iNode++)
    //         // {
    //         //     int i, j, k;
    //         //     m_idx_proxy->GetIdxStruct(m_bound_patch[iNode].idx_master, i, j, k);
    //         //     double coord[3];
    //         //     coord[0] = node->GetCoord(i, j, k)[0];
    //         //     coord[1] = node->GetCoord(i, j, k)[1];
    //         //     coord[2] = node->GetCoord(i, j, k)[2];
    //         //     m_model_manager->GetClosestPoint(coord, m_bound_patch[iNode].coord);
    //         // }
    //     }
    //     void NSFieldZaran::BuildTransFace()
    //     {
    //         // for (int iFace = 0; iFace < m_trans_face.size(); iFace++)
    //         // {
    //         //     auto &idx_master = m_trans_face[iFace].idx_master;
    //         //     auto &idx_slave = m_trans_face[iFace].idx_slave;
    //         //     idx_slave.resize(idx_master.size());
    //         //     for (int iNode = 0; iNode < idx_master.size(); iNode++)
    //         //     {
    //         //         for (int iPatch = 0; iPatch < m_bound_patch.size(); iPatch++)
    //         //         {
    //         //             if (m_bound_patch[iPatch].idx_master == idx_master[iNode])
    //         //             {
    //         //                 idx_slave[iNode] = iPatch;
    //         //                 break;
    //         //             }
    //         //         }
    //         //     }
    //         // }
    //     }
} // namespace zaran