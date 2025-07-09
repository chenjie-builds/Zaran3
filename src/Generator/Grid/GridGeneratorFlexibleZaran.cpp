#include "GridGeneratorFlexibleZaran.h"
#include "Log.h"
#include "MathBasic.h"
#include <cmath>
#include <omp.h>
namespace zaran
{
void GridFNFactoryZaran::CreateGrid(const shared_ptr<GridBlock> &block, const shared_ptr<GridFN> &grid,
                                    const shared_ptr<ModelManager> &model_manager)
{
    m_model_manager = model_manager;
    m_block_grid = block;
    m_fn_grid = grid;
    count_type ni = m_block_grid->GetNi();
    count_type nj = m_block_grid->GetNj();
    count_type nk = m_block_grid->GetNk();
    m_idx_proxy = make_shared<IdProxyStruct>(ni, nj, nk);
    m_layer_num = GlobalData::GetInt("projection_layer");
    TagCells();
    TagNodes();
    Log::info("CheckTransNode...");
    while (!CheckTransNode())
    {
        ReTagBlockGrid();
    }
    Log::info("CheckTransNode success");
    CheckTransFace();
    SetNodeTag();
    Log::info("BuildFNNodeCoord");
    BuildFNNodeCoord();
    Log::info("SetFNGridNode");
    SetFNGridNode();
    WriteProjectNode();
    WriteModelSurface();
    WriteTransFace();
    Log::info("BuildCell");
    BuildFNCell();
    Log::info("SetFNGridCell");
    SetFNGridCell();
    WriteSlaveGrid();
    Log::info("BuildNodeNeighbor");
    BuildNodeNeighbor();
    Log::info("SetFNGridNodeNeighbor");
    SetFNGridNodeNeighbor();
    Log::info("SetFNGridBoundary");
    SetFNGridBoundary();
    Log::info("SetFNGridBoundaryFace");
    SetFNGridBoundaryFace();
}

void GridFNFactoryZaran::TagBlockGrid()
{
    TagCells();
    TagNodes();
}

void GridFNFactoryZaran::TagCells()
{
    const auto grid = GetBlockGrid();
    const count_type ni = grid->GetNi();
    const count_type nj = grid->GetNj();
    const count_type nk = grid->GetNk();
    m_cell_type.resize(ni * nj * nk);
    for (int iCell = 0; iCell < ni * nj * nk; iCell++)
    {
        m_cell_type[iCell] = PhysicalType::Unset;
    }
    const double dx = grid->GetDx();
    const double dy = grid->GetDy();
    const double dz = grid->GetDz();
    const double tol_factor = GlobalData::GetDouble("tol_factor");
    double tol = tol_factor * sqrt(dx * dx + dy * dy + dz * dz);
    const auto model_manager = GetModelManager();
    auto &box = grid->GetBoundBox();
    auto &model_box = model_manager->GetBox();
    IdProxyStruct &idx_proxy = *m_idx_proxy;
    ProcessCell(0, ni, 0, nj, 0, nk);
    // #pragma omp parallel for
    //         for (int iCell = 0; iCell < ni * nj * nk; iCell++)
    //         {
    //             // if(iCell % 10==0)
    //             // {
    //             //     Log::info("thread_id:{}, iCell={}", thread_id, iCell);
    //             // }
    //             int i, j, k;
    //             idx_proxy->GetIdxStruct(iCell, i, j, k);
    //             double cell_center[3];
    //             cell_center[0] = box.x_min + (i - ghost_size + 0.5) * dx;
    //             cell_center[1] = box.y_min + (j - ghost_size + 0.5) * dy;
    //             cell_center[2] = box.z_min + (k - ghost_size + 0.5) * dz;
    //             if (cell_center[0] < model_box.x_min - tol || cell_center[0] >
    //             model_box.x_max + tol || cell_center[1] < model_box.y_min - tol
    //             || cell_center[1] > model_box.y_max + tol || cell_center[2] <
    //             model_box.z_min - tol || cell_center[2] > model_box.z_max +
    //             tol)
    //             {
    //                 m_cell_type[iCell] = PhysicalType::Fluid;
    //             }
    //             else
    //             {
    //                 double dist =
    //                 model_manager->GetClosestDistance(cell_center); if (dist <
    //                 tol)
    //                 {
    //                     m_cell_type[iCell] = PhysicalType::Solid;
    //                 }
    //                 else
    //                 {
    //                     m_cell_type[iCell] = PhysicalType::Unset;
    //                 }
    //             }
    //         }

    auto IsValidCell = [&](int i, int j, int k) -> bool {
        return i >= 0 && i < ni && j >= 0 && j < nj && k >= 0 && k < nk;
    };
    m_cell_type[idx_proxy(0, 0, 0)] = PhysicalType::Fluid;
    int new_fluid_cell = 1;
    while (new_fluid_cell != 0)
    {
        new_fluid_cell = 0;
        for (int k = 0; k < nk; k++)
        {
            for (int j = 0; j < nj; j++)
            {
                for (int i = 0; i < ni; i++)
                {
                    int idx = idx_proxy(i, j, k);
                    if (m_cell_type[idx] == PhysicalType::Fluid)
                    {
                        if (IsValidCell(i + 1, j, k) && m_cell_type[idx_proxy(i + 1, j, k)] == PhysicalType::Unset)
                        {
                            m_cell_type[idx_proxy(i + 1, j, k)] = PhysicalType::Fluid;
                            new_fluid_cell++;
                        }
                        if (IsValidCell(i - 1, j, k) && m_cell_type[idx_proxy(i - 1, j, k)] == PhysicalType::Unset)
                        {
                            m_cell_type[idx_proxy(i - 1, j, k)] = PhysicalType::Fluid;
                            new_fluid_cell++;
                        }
                        if (IsValidCell(i, j + 1, k) && m_cell_type[idx_proxy(i, j + 1, k)] == PhysicalType::Unset)
                        {
                            m_cell_type[idx_proxy(i, j + 1, k)] = PhysicalType::Fluid;
                            new_fluid_cell++;
                        }
                        if (IsValidCell(i, j - 1, k) && m_cell_type[idx_proxy(i, j - 1, k)] == PhysicalType::Unset)
                        {
                            m_cell_type[idx_proxy(i, j - 1, k)] = PhysicalType::Fluid;
                            new_fluid_cell++;
                        }
                        if (IsValidCell(i, j, k + 1) && m_cell_type[idx_proxy(i, j, k + 1)] == PhysicalType::Unset)
                        {
                            m_cell_type[idx_proxy(i, j, k + 1)] = PhysicalType::Fluid;
                            new_fluid_cell++;
                        }
                        if (IsValidCell(i, j, k - 1) && m_cell_type[idx_proxy(i, j, k - 1)] == PhysicalType::Unset)
                        {
                            m_cell_type[idx_proxy(i, j, k - 1)] = PhysicalType::Fluid;
                            new_fluid_cell++;
                        }
                    }
                }
            }
        }
    }

    // tag the solid cells after the fluid cells
    for (int i = 0; i < ni - 1; i++)
    {
        for (int j = 0; j < nj - 1; j++)
        {
            for (int k = 0; k < nk - 1; k++)
            {
                int idx = idx_proxy(i, j, k);
                if (m_cell_type[idx] == PhysicalType::Unset)
                {
                    m_cell_type[idx] = PhysicalType::Solid;
                }
            }
        }
    }
    int new_solid_num = 1;
    while (new_solid_num > 0)
    {
        new_solid_num = 0;
        for (int i = 0; i < ni - 1; i++)
        {
            for (int j = 0; j < nj - 1; j++)
            {
                for (int k = 0; k < -1; k++)
                {
                    int idx = idx_proxy(i, j, k);
                    if (m_cell_type[idx] != PhysicalType::Fluid)
                        continue;
                    if (m_cell_type[idx_proxy(i + 1, j, k)] == PhysicalType::Solid &&
                        m_cell_type[idx_proxy(i - 1, j, k)] == PhysicalType::Solid)
                    {
                        m_cell_type[idx] = PhysicalType::Solid;
                        new_solid_num++;
                    }
                    else if (m_cell_type[idx_proxy(i, j + 1, k)] == PhysicalType::Solid &&
                             m_cell_type[idx_proxy(i, j - 1, k)] == PhysicalType::Solid)
                    {
                        m_cell_type[idx] = PhysicalType::Solid;
                        new_solid_num++;
                    }
                    else if (grid->GetDim() == 3 && m_cell_type[idx_proxy(i, j, k + 1)] == PhysicalType::Solid &&
                             m_cell_type[idx_proxy(i, j, k - 1)] == PhysicalType::Solid)
                    {
                        m_cell_type[idx] = PhysicalType::Solid;
                        new_solid_num++;
                    }
                }
            }
        }
        if (new_solid_num > 0)
            Log::info("new_solid_num={}", new_solid_num);
    }

    // tag the fluid-solid cells
    for (int k = 1; k < nk - 1; k++)
    {
        for (int j = 1; j < nj - 1; j++)
        {
            for (int i = 1; i < ni - 1; i++)
            {
                int idx = idx_proxy(i, j, k);
                if (m_cell_type[idx] != PhysicalType::Fluid)
                    continue;
                if (m_cell_type[idx_proxy(i + 1, j, k)] == PhysicalType::Solid ||
                    m_cell_type[idx_proxy(i + 1, j, k)] == PhysicalType::Unset)
                {
                    m_cell_type[idx] = PhysicalType::FluidSolid;
                }
                if (m_cell_type[idx_proxy(i - 1, j, k)] == PhysicalType::Solid ||
                    m_cell_type[idx_proxy(i - 1, j, k)] == PhysicalType::Unset)
                {
                    m_cell_type[idx] = PhysicalType::FluidSolid;
                }
                if (m_cell_type[idx_proxy(i, j + 1, k)] == PhysicalType::Solid ||
                    m_cell_type[idx_proxy(i, j + 1, k)] == PhysicalType::Unset)
                {
                    m_cell_type[idx] = PhysicalType::FluidSolid;
                }
                if (m_cell_type[idx_proxy(i, j - 1, k)] == PhysicalType::Solid ||
                    m_cell_type[idx_proxy(i, j - 1, k)] == PhysicalType::Unset)
                {
                    m_cell_type[idx] = PhysicalType::FluidSolid;
                }
                if (grid->GetDim() == 3)
                {

                    if (m_cell_type[idx_proxy(i, j, k + 1)] == PhysicalType::Solid ||
                        m_cell_type[idx_proxy(i, j, k + 1)] == PhysicalType::Unset)
                    {
                        m_cell_type[idx] = PhysicalType::FluidSolid;
                    }
                    if (m_cell_type[idx_proxy(i, j, k - 1)] == PhysicalType::Solid ||
                        m_cell_type[idx_proxy(i, j, k - 1)] == PhysicalType::Unset)
                    {
                        m_cell_type[idx] = PhysicalType::FluidSolid;
                    }
                }
            }
        }
    }
    WriteCellTag();
}
void GridFNFactoryZaran::TagNodes()
{
    auto grid = GetBlockGrid();
    if (grid->GetDim() == 2)
    {
        TagNodes2D();
    }
    else if (grid->GetDim() == 3)
    {
        TagNodes3D();
    }
    else
    {
        Log::error("GridFNFactoryZaran::TagNodes: Unsupported grid dimension: {}", grid->GetDim());
        throw std::runtime_error("Unsupported grid dimension");
    }
}
void GridFNFactoryZaran::TagNodes3D()
{
    auto grid = GetBlockGrid();
    auto node = grid->GetNode();
    IdProxyStruct &idx_proxy = *m_idx_proxy;
    int ni, nj, nk;
    ni = grid->GetNi();
    nj = grid->GetNj();
    nk = grid->GetNk();
    auto IsValidNode = [&](int i, int j, int k) -> bool {
        return i >= 0 && i < ni && j >= 0 && j < nj && k >= 0 && k < nk;
    };
    auto IsValidCell = [&](int i, int j, int k) -> bool {
        return i >= 0 && i < ni - 1 && j >= 0 && j < nj - 1 && k >= 0 && k < nk - 1;
    };
    index_type is, ie, js, je, ks, ke;
    grid->GetRange(is, ie, js, je, ks, ke);
    m_node_type.resize(ni * nj * nk);
    for (int idx = 0; idx < ni * nj * nk; idx++)
    {
        index_type i, j, k;
        idx_proxy.GetIdxStruct(idx, i, j, k);
        if (i < is || i >= ie || j < js || j >= je || k < ks || k >= ke)
        {
            m_node_type[idx] = PhysicalType::Fluid;
        }
        else
            m_node_type[idx] = PhysicalType::Unset;
    }
    std::set<int> trans_node_set;
    std::set<TransFace> trans_face_set;
    // find the fluid-solid nodes
    for (int k = ks; k <= ke; k++)
    {
        for (int j = js; j <= je; j++)
        {
            for (int i = is; i <= ie; i++)
            {
                auto idx = idx_proxy(i, j, k);
                if (IsValidCell(i - 1, j, k))
                {
                    int idx_left = idx_proxy(i - 1, j, k);
                    if (m_cell_type[idx_left] == PhysicalType::FluidSolid && m_cell_type[idx] == PhysicalType::Fluid ||
                        m_cell_type[idx_left] == PhysicalType::Fluid && m_cell_type[idx] == PhysicalType::FluidSolid)
                    {
                        trans_face_set.insert(TransFace{{idx_proxy(i, j, k), idx_proxy(i, j + 1, k),
                                                         idx_proxy(i, j + 1, k + 1), idx_proxy(i, j, k + 1)}});
                        m_node_type[idx_proxy(i, j, k)] = PhysicalType::FluidSolid;
                        m_node_type[idx_proxy(i, j + 1, k)] = PhysicalType::FluidSolid;
                        m_node_type[idx_proxy(i, j, k + 1)] = PhysicalType::FluidSolid;
                        m_node_type[idx_proxy(i, j + 1, k + 1)] = PhysicalType::FluidSolid;
                    }
                }
                if (IsValidCell(i, j - 1, k))
                {
                    int idx_down = idx_proxy(i, j - 1, k);
                    if (m_cell_type[idx_down] == PhysicalType::FluidSolid && m_cell_type[idx] == PhysicalType::Fluid ||
                        m_cell_type[idx_down] == PhysicalType::Fluid && m_cell_type[idx] == PhysicalType::FluidSolid)
                    {
                        trans_face_set.insert(TransFace{{idx_proxy(i, j, k), idx_proxy(i + 1, j, k),
                                                         idx_proxy(i + 1, j, k + 1), idx_proxy(i, j, k + 1)}});
                        m_node_type[idx_proxy(i, j, k)] = PhysicalType::FluidSolid;
                        m_node_type[idx_proxy(i + 1, j, k)] = PhysicalType::FluidSolid;
                        m_node_type[idx_proxy(i, j, k + 1)] = PhysicalType::FluidSolid;
                        m_node_type[idx_proxy(i + 1, j, k + 1)] = PhysicalType::FluidSolid;
                    }
                }
                if (IsValidCell(i, j, k - 1))
                {
                    int idx_back = idx_proxy(i, j, k - 1);
                    if (m_cell_type[idx_back] == PhysicalType::FluidSolid && m_cell_type[idx] == PhysicalType::Fluid ||
                        m_cell_type[idx_back] == PhysicalType::Fluid && m_cell_type[idx] == PhysicalType::FluidSolid)
                    {
                        trans_face_set.insert(TransFace{{idx_proxy(i, j, k), idx_proxy(i + 1, j, k),
                                                         idx_proxy(i + 1, j + 1, k), idx_proxy(i, j + 1, k)}});
                        m_node_type[idx_proxy(i, j, k)] = PhysicalType::FluidSolid;
                        m_node_type[idx_proxy(i + 1, j, k)] = PhysicalType::FluidSolid;
                        m_node_type[idx_proxy(i, j + 1, k)] = PhysicalType::FluidSolid;
                        m_node_type[idx_proxy(i + 1, j + 1, k)] = PhysicalType::FluidSolid;
                    }
                }
            }
        }
    }
    m_trans_face.clear();
    m_trans_face.resize(trans_face_set.size());
    int idx = 0;
    for (auto &trans : trans_face_set)
    {
        m_trans_face[idx++].idx_block = trans.idx_block;
    }

    m_node_type[idx_proxy(0, 0, 0)] = PhysicalType::Fluid;
    // find the fluid nodes
    for (int k = 0; k < nk; k++)
    {
        for (int j = 0; j < nj; j++)
        {
            for (int i = 0; i < ni; i++)
            {
                int idx = idx_proxy(i, j, k);
                if (m_node_type[idx] != PhysicalType::Fluid)
                    continue;
                if (IsValidNode(i + 1, j, k))
                {
                    if (m_node_type[idx_proxy(i + 1, j, k)] == PhysicalType::Unset)
                        m_node_type[idx_proxy(i + 1, j, k)] = PhysicalType::Fluid;
                }
                if (IsValidNode(i, j + 1, k))
                {
                    if (m_node_type[idx_proxy(i, j + 1, k)] == PhysicalType::Unset)
                        m_node_type[idx_proxy(i, j + 1, k)] = PhysicalType::Fluid;
                }
                if (IsValidNode(i, j, k + 1))
                {
                    if (m_node_type[idx_proxy(i, j, k + 1)] == PhysicalType::Unset)
                        m_node_type[idx_proxy(i, j, k + 1)] = PhysicalType::Fluid;
                }
            }
        }
    }
    for (int iNode = 0; iNode < ni * nj * nk; iNode++)
    {
        if (m_node_type[iNode] == PhysicalType::Unset)
        {
            m_node_type[iNode] = PhysicalType::Solid;
        }
    }
}
void GridFNFactoryZaran::TagNodes2D()
{
    auto grid = GetBlockGrid();
    auto node = grid->GetNode();
    IdProxyStruct &idx_proxy = *m_idx_proxy;
    int ni, nj, nk;
    ni = grid->GetNi();
    nj = grid->GetNj();
    nk = grid->GetNk();
    auto IsValidNode = [&](int i, int j) -> bool { return i >= 0 && i < ni && j >= 0 && j < nj; };
    auto IsValidCell = [&](int i, int j) -> bool { return i >= 0 && i < ni - 1 && j >= 0 && j < nj - 1; };
    index_type is, ie, js, je, ks, ke;
    grid->GetRange(is, ie, js, je, ks, ke);
    m_node_type.resize(ni * nj * nk);
    for (int idx = 0; idx < ni * nj * nk; idx++)
    {
        index_type i, j, k;
        idx_proxy.GetIdxStruct(idx, i, j, k);

        m_node_type[idx] = PhysicalType::Unset;
    }
    std::set<int> trans_node_set;
    std::set<TransFace> trans_face_set;
    // find the fluid-solid nodes
    for (int k = ks; k <= ke; k++)
    {
        for (int j = js; j <= je; j++)
        {
            for (int i = is; i <= ie; i++)
            {
                auto idx = idx_proxy(i, j, k);
                if (IsValidCell(i - 1, j))
                {
                    int idx_left = idx_proxy(i - 1, j, k);
                    if (m_cell_type[idx_left] == PhysicalType::FluidSolid && m_cell_type[idx] == PhysicalType::Fluid ||
                        m_cell_type[idx_left] == PhysicalType::Fluid && m_cell_type[idx] == PhysicalType::FluidSolid)
                    {
                        trans_face_set.insert(TransFace{{idx_proxy(i, j, k), idx_proxy(i, j + 1, k)}});
                        m_node_type[idx_proxy(i, j, k)] = PhysicalType::FluidSolid;
                        m_node_type[idx_proxy(i, j + 1, k)] = PhysicalType::FluidSolid;
                    }
                }
                if (IsValidCell(i, j - 1))
                {
                    int idx_down = idx_proxy(i, j - 1, k);
                    if (m_cell_type[idx_down] == PhysicalType::FluidSolid && m_cell_type[idx] == PhysicalType::Fluid ||
                        m_cell_type[idx_down] == PhysicalType::Fluid && m_cell_type[idx] == PhysicalType::FluidSolid)
                    {
                        trans_face_set.insert(TransFace{{idx_proxy(i, j, k), idx_proxy(i + 1, j, k)}});
                        m_node_type[idx_proxy(i, j, k)] = PhysicalType::FluidSolid;
                        m_node_type[idx_proxy(i + 1, j, k)] = PhysicalType::FluidSolid;
                    }
                }
            }
        }
    }
    m_trans_face.clear();
    m_trans_face.resize(trans_face_set.size());
    int idx = 0;
    for (auto &trans : trans_face_set)
    {
        m_trans_face[idx++].idx_block = trans.idx_block;
    }

    m_node_type[idx_proxy(0, 0, ks)] = PhysicalType::Fluid;
    // find the fluid nodes
    while (true)
    {
        int new_fluid_node = 0;
        for (int k = ks; k <= ke; k++)
        {
            for (int j = 0; j < nj; j++)
            {
                for (int i = 0; i < ni; i++)
                {
                    int idx = idx_proxy(i, j, k);
                    if (m_node_type[idx] != PhysicalType::Fluid)
                        continue;
                    if (IsValidNode(i + 1, j))
                    {
                        if (m_node_type[idx_proxy(i + 1, j, k)] == PhysicalType::Unset)
                        {
                            m_node_type[idx_proxy(i + 1, j, k)] = PhysicalType::Fluid;
                            new_fluid_node++;
                        }
                    }
                    if (IsValidNode(i - 1, j))
                    {
                        if (m_node_type[idx_proxy(i - 1, j, k)] == PhysicalType::Unset)
                        {
                            m_node_type[idx_proxy(i - 1, j, k)] = PhysicalType::Fluid;
                            new_fluid_node++;
                        }
                    }
                    if (IsValidNode(i, j + 1))
                    {
                        if (m_node_type[idx_proxy(i, j + 1, k)] == PhysicalType::Unset)
                        {
                            m_node_type[idx_proxy(i, j + 1, k)] = PhysicalType::Fluid;
                            new_fluid_node++;
                        }
                    }
                    if (IsValidNode(i, j - 1))
                    {
                        if (m_node_type[idx_proxy(i, j - 1, k)] == PhysicalType::Unset)
                        {
                            m_node_type[idx_proxy(i, j - 1, k)] = PhysicalType::Fluid;
                            new_fluid_node++;
                        }
                    }
                }
            }
        }
        if (new_fluid_node == 0)
            break;
    }
    for (int k = 0; k < nk; k++)
    {
        for (int j = 0; j < nj; j++)
        {
            for (int i = 0; i < ni; i++)
            {
                m_node_type[idx_proxy(i, j, k)] = m_node_type[idx_proxy(i, j, ks)];
            }
        }
    }
    for (int iNode = 0; iNode < ni * nj * nk; iNode++)
    {
        if (m_node_type[iNode] == PhysicalType::Unset)
        {
            m_node_type[iNode] = PhysicalType::Solid;
        }
    }
}
void GridFNFactoryZaran::ReTagBlockGrid()
{
    ReTagCells();
    TagNodes();
}

void GridFNFactoryZaran::ProcessCell(const index_type start_i, const index_type end_i, const index_type start_j,
                                     const index_type end_j, const index_type start_k, const index_type end_k)
{
    // Log::info("start_i={}, end_i={}, start_j={}, end_j={}, start_k={},
    // end_k={}", start_i, end_i, start_j, end_j, start_k, end_k);
    const auto grid = GetBlockGrid();
    IdProxyStruct &idx_proxy = grid->GetIdxProxy();
    auto &grid_box = grid->GetBoundBox();
    const auto model_manager = GetModelManager();
    const auto &model_box = model_manager->GetBox();
    const index_type mid_i = (start_i + end_i) / 2;
    const index_type mid_j = (start_j + end_j) / 2;
    const index_type mid_k = (start_k + end_k) / 2;
    double x_min = grid_box.x_min + (int(start_i) - grid->GetGhostLevel() + 0.5) * grid->GetDx();
    double x_max = grid_box.x_min + (int(end_i) - 1 - grid->GetGhostLevel() + 0.5) * grid->GetDx();
    double y_min = grid_box.y_min + (int(start_j) - grid->GetGhostLevel() + 0.5) * grid->GetDy();
    double y_max = grid_box.y_min + (int(end_j) - 1 - grid->GetGhostLevel() + 0.5) * grid->GetDy();
    double z_min = grid_box.z_min + (int(start_k) - grid->GetGhostLevel() + 0.5) * grid->GetDz();
    double z_max = grid_box.z_min + (int(end_k) - 1 - grid->GetGhostLevel() + 0.5) * grid->GetDz();
    if (grid->GetDim() == 2)
    {
        z_min = grid_box.z_min;
        z_max = grid_box.z_max;
    }
    double cell_center[3];
    cell_center[0] = 0.5 * (x_min + x_max);
    cell_center[1] = 0.5 * (y_min + y_max);
    cell_center[2] = 0.5 * (z_min + z_max);
    double tol_x = x_max - x_min;
    double tol_y = y_max - y_min;
    double tol_factor = GlobalData::GetDouble("tol_factor");
    double tol, tol1;
    if (grid->GetDim() == 2)
    {
        tol = tol_factor * sqrt(tol_x * tol_x + tol_y * tol_y);
        tol1 = tol_factor * sqrt(grid->GetDx() * grid->GetDx() + grid->GetDy() * grid->GetDy());
    }
    else
    {
        double tol_z = z_max - z_min;
        tol = tol_factor * sqrt(tol_x * tol_x + tol_y * tol_y + tol_z * tol_z);
        tol1 = tol_factor *
               sqrt(grid->GetDx() * grid->GetDx() + grid->GetDy() * grid->GetDy() + grid->GetDz() * grid->GetDz());
    }
    tol += tol1;
    PhysicalType cell_type = PhysicalType::Unset;
    double dist = 0;
    if (x_max < model_box.x_min - tol1 || x_min > model_box.x_max + tol1 || y_max < model_box.y_min - tol1 ||
        y_min > model_box.y_max + tol1)
    {
        if (grid->GetDim() == 2)
        {
            cell_type = PhysicalType::Unset;
        }
        else if (z_max < model_box.z_min - tol1 || z_min > model_box.z_max + tol1)
        {
            cell_type = PhysicalType::Unset;
        }
    }
    else
    {
        dist = model_manager->GetClosestDistance(cell_center);
        if (dist < tol)
        {
            cell_type = PhysicalType::Solid;
        }
        else
        {
            cell_type = PhysicalType::Unset;
        }
    }
    if (cell_type != PhysicalType::Solid)
    {
        return;
    }
    if (start_i == mid_i)
    {
        if (start_j == mid_j)
        {
            if (start_k == mid_k)
            {
                int idx = idx_proxy(start_i, start_j, start_k);
                if (cell_type == PhysicalType::Solid)
                {
                    m_cell_type[idx] = PhysicalType::Solid;
                }
                return;
            }
            else
            {
                ProcessCell(start_i, end_i, start_j, end_j, start_k, mid_k);
                ProcessCell(start_i, end_i, start_j, end_j, mid_k, end_k);
            }
        }
        else
        {
            if (start_k == mid_k)
            {
                ProcessCell(start_i, end_i, start_j, mid_j, start_k, end_k);
                ProcessCell(start_i, end_i, mid_j, end_j, start_k, end_k);
            }
            else
            {
                ProcessCell(start_i, end_i, start_j, mid_j, start_k, mid_k);
                ProcessCell(start_i, end_i, mid_j, end_j, start_k, mid_k);
                ProcessCell(start_i, end_i, start_j, mid_j, mid_k, end_k);
                ProcessCell(start_i, end_i, mid_j, end_j, mid_k, end_k);
            }
        }
    }
    else
    {
        if (start_j == mid_j)
        {
            if (start_k == mid_k)
            {
                ProcessCell(start_i, mid_i, start_j, end_j, start_k, end_k);
                ProcessCell(mid_i, end_i, start_j, end_j, start_k, end_k);
            }
            else
            {
                ProcessCell(start_i, mid_i, start_j, end_j, start_k, mid_k);
                ProcessCell(mid_i, end_i, start_j, end_j, start_k, mid_k);
                ProcessCell(start_i, mid_i, start_j, end_j, mid_k, end_k);
                ProcessCell(mid_i, end_i, start_j, end_j, mid_k, end_k);
            }
        }
        else
        {
            if (start_k == mid_k)
            {
                ProcessCell(start_i, mid_i, start_j, mid_j, start_k, end_k);
                ProcessCell(mid_i, end_i, start_j, mid_j, start_k, end_k);
                ProcessCell(start_i, mid_i, mid_j, end_j, start_k, end_k);
                ProcessCell(mid_i, end_i, mid_j, end_j, start_k, end_k);
            }
            else
            {
                ProcessCell(start_i, mid_i, start_j, mid_j, start_k, mid_k);
                ProcessCell(mid_i, end_i, start_j, mid_j, start_k, mid_k);
                ProcessCell(start_i, mid_i, mid_j, end_j, start_k, mid_k);
                ProcessCell(mid_i, end_i, mid_j, end_j, start_k, mid_k);
                ProcessCell(start_i, mid_i, start_j, mid_j, mid_k, end_k);
                ProcessCell(mid_i, end_i, start_j, mid_j, mid_k, end_k);
                ProcessCell(start_i, mid_i, mid_j, end_j, mid_k, end_k);
                ProcessCell(mid_i, end_i, mid_j, end_j, mid_k, end_k);
            }
        }
    }
}

void GridFNFactoryZaran::ReTagCells()
{
    // 去掉所有的FluidSolid标记
    for (int iCell = 0; iCell < m_cell_type.size(); iCell++)
    {
        if (m_cell_type[iCell] == PhysicalType::FluidSolid)
        {
            m_cell_type[iCell] = PhysicalType::Fluid;
        }
    }

    auto grid = GetBlockGrid();
    IdProxyStruct &idx_proxy = *m_idx_proxy;
    int ni, nj, nk;
    ni = grid->GetNi();
    nj = grid->GetNj();
    nk = grid->GetNk();
    // 检查是否有Solid单元被标记为Fluid
    int new_solid_num = 1;
    while (new_solid_num > 0)
    {
        new_solid_num = 0;
        for (int i = 0; i < ni - 1; i++)
        {
            for (int j = 0; j < nj - 1; j++)
            {
                for (int k = 0; k < -1; k++)
                {
                    int idx = idx_proxy(i, j, k);
                    if (m_cell_type[idx] != PhysicalType::Fluid)
                        continue;
                    if (m_cell_type[idx_proxy(i + 1, j, k)] == PhysicalType::Solid &&
                        m_cell_type[idx_proxy(i - 1, j, k)] == PhysicalType::Solid)
                    {
                        m_cell_type[idx] = PhysicalType::Solid;
                        new_solid_num++;
                    }
                    else if (m_cell_type[idx_proxy(i, j + 1, k)] == PhysicalType::Solid &&
                             m_cell_type[idx_proxy(i, j - 1, k)] == PhysicalType::Solid)
                    {
                        m_cell_type[idx] = PhysicalType::Solid;
                        new_solid_num++;
                    }
                    else if (m_cell_type[idx_proxy(i, j, k + 1)] == PhysicalType::Solid &&
                             m_cell_type[idx_proxy(i, j, k - 1)] == PhysicalType::Solid)
                    {
                        m_cell_type[idx] = PhysicalType::Solid;
                        new_solid_num++;
                    }
                }
            }
        }
        if (new_solid_num > 0)
            Log::info("new_solid_num={}", new_solid_num);
    }

    // retag the fluid-solid cells
    for (int k = 1; k < nk - 1; k++)
    {
        for (int j = 1; j < nj - 1; j++)
        {
            for (int i = 1; i < ni - 1; i++)
            {
                int idx = idx_proxy(i, j, k);
                if (m_cell_type[idx] != PhysicalType::Fluid)
                    continue;
                if (m_cell_type[idx_proxy(i + 1, j, k)] == PhysicalType::Solid)
                {
                    m_cell_type[idx] = PhysicalType::FluidSolid;
                }
                if (m_cell_type[idx_proxy(i - 1, j, k)] == PhysicalType::Solid)
                {
                    m_cell_type[idx] = PhysicalType::FluidSolid;
                }
                if (m_cell_type[idx_proxy(i, j + 1, k)] == PhysicalType::Solid)
                {
                    m_cell_type[idx] = PhysicalType::FluidSolid;
                }
                if (m_cell_type[idx_proxy(i, j - 1, k)] == PhysicalType::Solid)
                {
                    m_cell_type[idx] = PhysicalType::FluidSolid;
                }
                if (m_cell_type[idx_proxy(i, j, k + 1)] == PhysicalType::Solid)
                {
                    m_cell_type[idx] = PhysicalType::FluidSolid;
                }
                if (m_cell_type[idx_proxy(i, j, k - 1)] == PhysicalType::Solid)
                {
                    m_cell_type[idx] = PhysicalType::FluidSolid;
                }
            }
        }
    }
}

void GridFNFactoryZaran::SetNodeTag()
{
    auto grid = GetBlockGrid();
    auto node = grid->GetNode();
    IdProxyStruct &idx_proxy = *m_idx_proxy;
    index_type is, ie, js, je, ks, ke;
    grid->GetRange(is, ie, js, je, ks, ke);
    for (int k = ks; k <= ke; k++)
    {
        for (int j = js; j <= je; j++)
        {
            for (int i = is; i <= ie; i++)
            {
                int idx = idx_proxy(i, j, k);
                if (m_node_type[idx] == PhysicalType::FluidSolid)
                {
                    grid->SetIBlank(i, j, k, IBlank::Trans);
                }
                else if (m_node_type[idx] == PhysicalType::Fluid)
                {
                    grid->SetIBlank(i, j, k, IBlank::Fluid);
                }
                else if (m_node_type[idx] == PhysicalType::Solid)
                {
                    grid->SetIBlank(i, j, k, IBlank::Solid);
                }
                else
                {
                    grid->SetIBlank(i, j, k, IBlank::Unset);
                }
            }
        }
    }
}

void GridFNFactoryZaran::WriteNodeTag()
{
    std::ofstream out("node_tag.dat");
    auto grid = GetBlockGrid();
    auto node = grid->GetNode();
    index_type is, ie, js, je, ks, ke;
    grid->GetRange(is, ie, js, je, ks, ke);
    out << "variables = x, y, z, type" << std::endl;
    out << "Zone T=test\n I=" << ie - is + 1 << " J=" << je - js + 1 << " K=" << ke - ks + 1 << " F=POINT" << std::endl;
    for (int k = ks; k <= ke; k++)
    {
        for (int j = js; j <= je; j++)
        {
            for (int i = is; i <= ie; i++)
            {
                int idx = grid->GetIdxProxy()(i, j, k);
                out << node->GetCoord(i, j, k)[0] << " " << node->GetCoord(i, j, k)[1] << " "
                    << node->GetCoord(i, j, k)[2] << " " << (int)m_node_type[idx] << std::endl;
            }
        }
    }
}

void GridFNFactoryZaran::WriteProjectNode()
{
    std::ofstream out("project_node.dat");
    out << "variables = x, y, z" << std::endl;
    out << "Zone T=test" << std::endl;
    for (auto &node : m_fn_info.node[m_layer_num])
    {
        out << node.coord[0] << " " << node.coord[1] << " " << node.coord[2] << std::endl;
    }
}

void GridFNFactoryZaran::WriteModelSurface()
{
    std::ofstream out("model_surface.dat");
    out << "variables = x, y, z" << std::endl;
    out << "Zone T=test" << std::endl;
    out << " N=" << m_fn_info.node[m_layer_num].size() << " E=" << m_trans_face.size() << " ZONETYPE=FEQuadrilateral"
        << std::endl;
    out << "DATAPACKING=POINT" << std::endl;
    for (auto &coord : m_fn_info.node[m_layer_num])
    {
        out << coord.coord[0] << " " << coord.coord[1] << " " << coord.coord[2] << std::endl;
    }
    for (int iFace = 0; iFace < m_trans_face.size(); iFace++)
    {
        // auto face_node_idx = m_trans_face[iFace].idx_master;
        // for (int iNode = 0; iNode < face_node_idx.size(); iNode++)
        // {
        //     for (auto &slave : m_trans_node)
        //     {
        //         if (slave.idx_block == face_node_idx[iNode])
        //         {
        //             face_node_idx[iNode] = slave.idx_local_layer;
        //             break;
        //         }
        //     }
        // }
        auto &face_node_idx = m_trans_face[iFace].idx_slave;

        int prior_node, next_node;
        for (int iNode = 0; iNode < face_node_idx.size(); iNode++)
        {
            out << face_node_idx[iNode] + 1 << " ";
        }
        out << std::endl;
    }
    out.close();
}

void GridFNFactoryZaran::WriteTransFace()
{
    std::ofstream out("trans_face.dat");
    out << "variables = x, y, z" << std::endl;
    out << "Zone T=test" << std::endl;
    auto grid = GetBlockGrid();
    if (grid->GetDim() == 2)
    {
        out << " N=" << m_fn_info.node[1].size() << " E=" << m_trans_face.size() << " ZONETYPE=FELineSeg" << std::endl;
    }
    else
        out << " N=" << m_fn_info.node[1].size() << " E=" << m_trans_face.size() << " ZONETYPE=FEQuadrilateral"
            << std::endl;
    out << "DATAPACKING=POINT" << std::endl;
    for (auto &coord : m_fn_info.node[1])
    {
        out << coord.coord[0] << " " << coord.coord[1] << " " << coord.coord[2] << std::endl;
    }
    for (int iFace = 0; iFace < m_trans_face.size(); iFace++)
    {
        // auto face_node_idx = m_trans_face[iFace].idx_master;
        // dynamic_array<int> idx_layer;
        // idx_layer.resize(face_node_idx.size());
        // for (int iNode = 0; iNode < face_node_idx.size(); iNode++)
        // {
        //     for (auto &slave : m_trans_node)
        //     {
        //         if (slave.idx_block == face_node_idx[iNode])
        //         {
        //             face_node_idx[iNode] = slave.idx_local_layer;
        //             idx_layer[iNode] = slave.idx_local_layer;
        //             break;
        //         }
        //     }
        // }
        auto &face_node_idx = m_trans_face[iFace].idx_slave;

        for (int iNode = 0; iNode < face_node_idx.size(); iNode++)
        {
            out << face_node_idx[iNode] + 1 << " ";
        }
        out << std::endl;
    }
    out.close();
}

void GridFNFactoryZaran::WriteSlaveGrid()
{
    std::ofstream out("slave_grid.dat");
    out << "variables = x, y, z,idx" << std::endl;
    out << "Zone T=test" << std::endl;
    int node_num = 0;
    for (int iLayer = 0; iLayer < m_fn_info.node.size(); iLayer++)
    {
        node_num += m_fn_info.node[iLayer].size();
    }
    int cell_num = 0;
    cell_num = m_fn_info.cell.size();
    auto grid = GetBlockGrid();
    if (grid->GetDim() == 2)
        out << " N=" << node_num << " E=" << cell_num << " F=FEPOINT, ET=QUADRILATERAL" << std::endl;
    else
        out << " N=" << node_num << " E=" << cell_num << " F=FEPOINT, ET=Brick" << std::endl;
    for (int iLayer = 0; iLayer < m_fn_info.node.size(); iLayer++)
    {
        for (int iNode = 0; iNode < m_fn_info.node[iLayer].size(); iNode++)
        {
            out << m_fn_info.node[iLayer][iNode].coord[0] << " " << m_fn_info.node[iLayer][iNode].coord[1] << " "
                << m_fn_info.node[iLayer][iNode].coord[2] << " " << m_fn_info.node[iLayer][iNode].idx << std::endl;
        }
    }
    for (int iCell = 0; iCell < m_fn_info.cell.size(); iCell++)
    {
        for (int iNode = 0; iNode < m_fn_info.cell[iCell].size(); iNode++)
        {
            out << m_fn_info.cell[iCell][iNode] + 1 << " ";
        }
        out << std::endl;
    }
    out.close();
}

void GridFNFactoryZaran::WriteCellTag()
{
    std::ofstream out("cell_tag.dat");
    auto grid = GetBlockGrid();
    auto node = grid->GetNode();
    index_type is, ie, js, je, ks, ke;
    grid->GetRange(is, ie, js, je, ks, ke);
    out << "variables = x, y, z, type" << std::endl;
    for (int k = ks; k <= ke; k++)
    {
        for (int j = js; j <= je - 1; j++)
        {
            for (int i = is; i <= ie - 1; i++)
            {
                auto idx = grid->GetIdxProxy()(i, j, k);
                out << node->GetCoord(i, j, k)[0] + grid->GetDx() << " " << node->GetCoord(i, j, k)[1] + grid->GetDy()
                    << " " << node->GetCoord(i, j, k)[2] << " " << (int)m_cell_type[idx] << std::endl;
            }
        }
    }
}

void GridFNFactoryZaran::BuildFNGridInfo()
{
    BuildFNNodeCoord();
    WriteProjectNode();
    WriteModelSurface();
    WriteTransFace();
    BuildFNCell();
    WriteSlaveGrid();
    BuildNodeNeighbor();
}

void GridFNFactoryZaran::BuildFNNodeCoord()
{

    m_fn_info.node.resize(m_layer_num + 1);
    BuildRefNode();
    BuildTransNode();
    BuildWallNode();
    OptimizeWallNode();
    BuildProjectNode();
}

void GridFNFactoryZaran::SetFNGrid()
{
    SetFNGridNode();
    SetFNGridNodeNeighbor();
    SetFNGridCell();
    SetFNGridBoundary();
    SetFNGridBoundaryFace();
}

void GridFNFactoryZaran::SetFNGridNode()
{
    auto grid = GetFNGrid();
    int node_num = 0;
    for (int iLayer = 0; iLayer < m_fn_info.node.size(); iLayer++)
    {
        node_num += m_fn_info.node[iLayer].size();
    }
    NodeFN &node = grid->GetNode();
    node.SetNodeNum(node_num);
    for (int iLayer = 0; iLayer < m_fn_info.node.size(); iLayer++)
    {
        for (int iNode = 0; iNode < m_fn_info.node[iLayer].size(); iNode++)
        {
            int idx = m_fn_info.node[iLayer][iNode].idx;
            auto coord = m_fn_info.node[iLayer][iNode].coord;
            node.SetCoord(idx, coord);
            if (iLayer == 0)
            {
                node.SetType(idx, NodeType::ghost);
            }
            else if (iLayer == m_layer_num)
            {
                node.SetType(idx, NodeType::wall);
            }
            else
            {
                node.SetType(idx, NodeType::inner);
            }
        }
    }
    grid->InitNode();
}

void GridFNFactoryZaran::SetFNGridNodeNeighbor()
{
    auto grid = GetFNGrid();
    NodeFN &node = grid->GetNode();
    index_type node_num = node.GetCount();
    dynamic_array<index_type> neighbor_node_num;
    neighbor_node_num.resize(node_num);
    dynamic_array<dynamic_array<index_type>> neighbor_node(node_num);
    for (index_type iLayer = 0; iLayer < m_fn_info.node.size(); iLayer++)
    {
        for (index_type iNode = 0; iNode < m_fn_info.node[iLayer].size(); iNode++)
        {
            index_type idx = m_fn_info.node[iLayer][iNode].idx;
            neighbor_node[idx] = m_fn_info.node[iLayer][iNode].neighbor_node;
        }
    }
    node.SetNeighborNode(neighbor_node);
}

void GridFNFactoryZaran::SetFNGridCell()
{
    auto grid = GetFNGrid();
    index_type cell_num = m_fn_info.cell.size();
    CellFN &cell = grid->GetCell();
    cell.SetCellNum(cell_num);
    cell.SetNode(m_fn_info.cell);
}

void GridFNFactoryZaran::SetFNGridBoundary()
{
    auto grid = GetFNGrid();
    auto &bound = grid->GetBoundaryMap();
    bound.AllocateBound("wall");
    auto &wall = bound.GetBound("wall");
    index_type node_num = m_fn_info.node[m_layer_num].size();
    wall.resize(node_num);
    for (index_type iNode = 0; iNode < node_num; iNode++)
    {
        index_type bound_idx = m_fn_info.node[m_layer_num][iNode].idx;
        index_type ref_idx = m_fn_info.node[m_layer_num - 1][iNode].idx;
        double normal[3];
        auto bound_coord = m_fn_info.node[m_layer_num][iNode].coord;
        auto ref_coord = m_fn_info.node[m_layer_num - 1][iNode].coord;
        normal[0] = bound_coord[0] - ref_coord[0];
        normal[1] = bound_coord[1] - ref_coord[1];
        normal[2] = bound_coord[2] - ref_coord[2];
        // normal[0] = -bound_coord[0];
        // normal[1] = -bound_coord[1];
        // normal[2] = -bound_coord[2];
        double len = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
        normal[0] /= len;
        normal[1] /= len;
        normal[2] /= len;
        wall[iNode].SetIdxBound(bound_idx);
        wall[iNode].SetIdxRef(ref_idx);
        wall[iNode].SetNorm(normal);
    }
}

void GridFNFactoryZaran::SetFNGridBoundaryFace()
{
    auto grid = GetFNGrid();
    index_type face_num = m_trans_face.size();
    dynamic_array<index_type> face_node_num;
    face_node_num.resize(face_num);
    for (index_type iFace = 0; iFace < face_num; iFace++)
    {
        face_node_num[iFace] = m_trans_face[iFace].idx_block.size();
    }
    FaceFN &face = grid->GetFace();
    face.Allocate(face_num, face_node_num.data());
    double normal[3];
    double area = 0;
    for (index_type iFace = 0; iFace < face_num; iFace++)
    {
        // auto face_node_idx = m_trans_face[iFace].idx_master;
        // for (int iNode = 0; iNode < face_node_idx.size(); iNode++)
        // {
        //     for (auto &slave : m_trans_node)
        //     {
        //         if (slave.idx_block == face_node_idx[iNode])
        //         {
        //             face_node_idx[iNode] = slave.idx_local_layer;
        //             break;
        //         }
        //     }
        // }
        auto &face_node_idx = m_trans_face[iFace].idx_slave;

        for (index_type iNode = 0; iNode < face_node_idx.size(); iNode++)
        {
            face_node_idx[iNode] = m_fn_info.node[m_layer_num][face_node_idx[iNode]].idx;
        }
        face.SetFace2Node(iFace, face_node_idx.data(), face_node_idx.size());
        face.SetNormal(iFace, normal);
        face.SetArea(iFace, area);
    }
}

void GridFNFactoryZaran::BuildRefNode()
{
    auto grid = GetBlockGrid();
    auto node = grid->GetNode();
    IdProxyStruct &idx_proxy = grid->GetIdxProxy();
    index_type is, ie, js, je, ks, ke;
    grid->GetRange(is, ie, js, je, ks, ke);
    std::set<index_type> ref_node_idx_set;
    for (int k = ks; k <= ke; k++)
    {
        for (int j = js; j <= je; j++)
        {
            for (int i = is; i <= ie; i++)
            {
                int idx = idx_proxy(i, j, k);
                if (m_node_type[idx] == PhysicalType::FluidSolid)
                {
                    if (m_node_type[idx_proxy(i + 1, j, k)] == PhysicalType::Fluid)
                    {
                        ref_node_idx_set.insert(idx_proxy(i + 1, j, k));
                    }
                    if (m_node_type[idx_proxy(i - 1, j, k)] == PhysicalType::Fluid)
                    {
                        ref_node_idx_set.insert(idx_proxy(i - 1, j, k));
                    }
                    if (m_node_type[idx_proxy(i, j + 1, k)] == PhysicalType::Fluid)
                    {
                        ref_node_idx_set.insert(idx_proxy(i, j + 1, k));
                    }
                    if (m_node_type[idx_proxy(i, j - 1, k)] == PhysicalType::Fluid)
                    {
                        ref_node_idx_set.insert(idx_proxy(i, j - 1, k));
                    }
                    if (m_node_type[idx_proxy(i, j, k + 1)] == PhysicalType::Fluid)
                    {
                        ref_node_idx_set.insert(idx_proxy(i, j, k + 1));
                    }
                    if (m_node_type[idx_proxy(i, j, k - 1)] == PhysicalType::Fluid)
                    {
                        ref_node_idx_set.insert(idx_proxy(i, j, k - 1));
                    }
                }
            }
        }
    }
    m_fn_info.node[0].resize(ref_node_idx_set.size());
    index_type idx = 0;
    for (auto &ref_node_idx : ref_node_idx_set)
    {
        index_type i, j, k;
        idx_proxy.GetIdxStruct(ref_node_idx, i, j, k);
        m_fn_info.node[0][idx].coord[0] = node->GetCoord(i, j, k)[0];
        m_fn_info.node[0][idx].coord[1] = node->GetCoord(i, j, k)[1];
        m_fn_info.node[0][idx].coord[2] = node->GetCoord(i, j, k)[2];
        m_fn_info.node[0][idx].idx = idx;
        m_ref_node.insert(ConnectInfo{ref_node_idx, 0, idx});
        idx++;
    }
}
void GridFNFactoryZaran::BuildTransNode()
{
    auto grid = GetBlockGrid();
    auto node = grid->GetNode();
    IdProxyStruct &idx_proxy = grid->GetIdxProxy();
    // int is, ie, js, je, ks, ke;
    // grid->GetRange(is, ie, js, je, ks, ke);
    std::set<index_type> trans_node_idx_set;
    // for (int k = ks; k <= ke; k++)
    // {
    //     for (int j = js; j <= je; j++)
    //     {
    //         for (int i = is; i <= ie; i++)
    //         {
    //             int idx = m_idx_proxy->GetIdx(i, j, k);
    //             if (m_node_type[idx] == PhysicalType::FluidSolid)
    //             {
    //                 trans_node_idx_set.insert(idx);
    //             }
    //         }
    //     }
    // }

    for (int iFace = 0; iFace < m_trans_face.size(); iFace++)
    {
        auto &face_node_idx = m_trans_face[iFace].idx_block;
        for (int iNode = 0; iNode < face_node_idx.size(); iNode++)
        {
            trans_node_idx_set.insert(face_node_idx[iNode]);
        }
    }

    m_fn_info.node[1].resize(trans_node_idx_set.size());
    index_type idx = 0;
    for (auto &trans_idx : trans_node_idx_set)
    {
        index_type i, j, k;
        idx_proxy.GetIdxStruct(trans_idx, i, j, k);
        m_fn_info.node[1][idx].coord[0] = node->GetCoord(i, j, k)[0];
        m_fn_info.node[1][idx].coord[1] = node->GetCoord(i, j, k)[1];
        m_fn_info.node[1][idx].coord[2] = node->GetCoord(i, j, k)[2];
        m_fn_info.node[1][idx].idx = idx + m_fn_info.node[0].size();
        m_trans_node.insert(ConnectInfo{trans_idx, 1, idx});
        idx++;
    }
    for (int iFace = 0; iFace < m_trans_face.size(); iFace++)
    {
        auto &idx_master = m_trans_face[iFace].idx_block;
        auto &idx_slave = m_trans_face[iFace].idx_slave;
        idx_slave.resize(idx_master.size());
        bool find = false;
        for (int iNode = 0; iNode < idx_slave.size(); iNode++)
        {
            for (auto &slave : m_trans_node)
            {
                if (slave.idx_block == idx_master[iNode])
                {
                    idx_slave[iNode] = slave.idx_local_layer;
                    find = true;
                    break;
                }
            }
            if (!find)
            {
                Log::error("Can not find trans node:{}", idx_master[iNode]);
            }
        }
    }
}
bool GridFNFactoryZaran::CheckTransNode()
{
    auto grid = GetBlockGrid();
    if (grid->GetDim() == 2)
    {
        return CheckTransNode2D();
    }
    else if (grid->GetDim() == 3)
    {
        return CheckTransNode3D();
    }
    else
    {
        Log::error("Grid dimension is not supported: {}", grid->GetDim());
        return false;
    }
}
bool GridFNFactoryZaran::CheckTransNode2D()
{
    auto grid = GetBlockGrid();
    index_type is, ie, js, je, ks, ke;
    grid->GetRange(is, ie, js, je, ks, ke);
    IdProxyStruct &idx_proxy = *m_idx_proxy;
    int total_error_num = 0;
    int error_num = 0;
    for (index_type k = ks; k <= ke; k++)
    {
        for (index_type j = js; j <= je; j++)
        {
            for (index_type i = is; i <= ie; i++)
            {
                index_type idx = idx_proxy(i, j, k);
                if (m_node_type[idx] != PhysicalType::FluidSolid)
                    continue;
                auto type_ip = m_node_type[idx_proxy(i + 1, j, k)];
                auto type_im = m_node_type[idx_proxy(i - 1, j, k)];
                auto type_jp = m_node_type[idx_proxy(i, j + 1, k)];
                auto type_jm = m_node_type[idx_proxy(i, j - 1, k)];
                int error_num = 0;
                if (type_ip == PhysicalType::Solid && type_im == PhysicalType::Solid)
                {
                    error_num++;
                }
                if (type_jp == PhysicalType::Solid && type_jm == PhysicalType::Solid)
                {
                    error_num++;
                }
                if (error_num > 0)
                {
                    m_cell_type[idx_proxy(i, j, k)] = PhysicalType::Solid;
                    m_cell_type[idx_proxy(i, j - 1, k)] = PhysicalType::Solid;
                    m_cell_type[idx_proxy(i - 1, j, k)] = PhysicalType::Solid;
                    m_cell_type[idx_proxy(i - 1, j - 1, k)] = PhysicalType::Solid;
                }
                total_error_num += error_num;
            }
        }
    }
    if (total_error_num > 0)
    {
        Log::warn("error trans node num:{}", total_error_num);
        return false;
    }
    return true;
}
bool GridFNFactoryZaran::CheckTransNode3D()
{
    IdProxyStruct &idx_proxy = *m_idx_proxy;
    int total_error_num = 0;
    for (int iNode = 0; iNode < m_node_type.size(); iNode++)
    {
        index_type i, j, k;
        idx_proxy.GetIdxStruct(iNode, i, j, k);
        if (m_node_type[iNode] != PhysicalType::FluidSolid)
            continue;
        auto type_ip = m_node_type[idx_proxy(i + 1, j, k)];
        auto type_im = m_node_type[idx_proxy(i - 1, j, k)];
        auto type_jp = m_node_type[idx_proxy(i, j + 1, k)];
        auto type_jm = m_node_type[idx_proxy(i, j - 1, k)];
        auto type_kp = m_node_type[idx_proxy(i, j, k + 1)];
        auto type_km = m_node_type[idx_proxy(i, j, k - 1)];
        int error_num = 0;
        if (type_ip == PhysicalType::Solid && type_im == PhysicalType::Solid)
        {
            // Log::warn("The trans node:({},{},{}) i+1 and i-1 are solid", i, j, k);
            error_num++;
        }
        if (type_jp == PhysicalType::Solid && type_jm == PhysicalType::Solid)
        {
            // Log::warn("The trans node:({},{},{}) j+1 and j-1 are solid", i, j, k);
            error_num++;
        }
        if (type_kp == PhysicalType::Solid && type_km == PhysicalType::Solid)
        {
            // Log::warn("The trans node:({},{},{}) k+1 and k-1 are solid", i, j, k);
            error_num++;
        }
        if (error_num > 0)
        {
            m_cell_type[idx_proxy(i, j, k)] = PhysicalType::Solid;
            m_cell_type[idx_proxy(i, j, k - 1)] = PhysicalType::Solid;
            m_cell_type[idx_proxy(i, j - 1, k)] = PhysicalType::Solid;
            m_cell_type[idx_proxy(i, j - 1, k - 1)] = PhysicalType::Solid;
            m_cell_type[idx_proxy(i - 1, j, k)] = PhysicalType::Solid;
            m_cell_type[idx_proxy(i - 1, j, k - 1)] = PhysicalType::Solid;
            m_cell_type[idx_proxy(i - 1, j - 1, k)] = PhysicalType::Solid;
            m_cell_type[idx_proxy(i - 1, j - 1, k - 1)] = PhysicalType::Solid;
        }
        total_error_num += error_num;
    }
    //if (total_error_num == 0)
    //{
    //    for (int iNode = 0; iNode < m_node_type.size(); iNode++)
    //    {
    //        auto grid = GetBlockGrid();
    //        auto node = grid->GetNode();
    //        index_type i, j, k;
    //        idx_proxy.GetIdxStruct(iNode, i, j, k);
    //        if (m_node_type[iNode] != PhysicalType::FluidSolid)
    //            continue;
    //        double dist = 0;
    //        auto coord = node->GetCoord(i, j, k);
    //        double wall_coord[3];
    //        m_model_manager->GetClosestPoint(coord, wall_coord);
    //        dist = DistanceOfTwoPoints(coord, wall_coord);
    //        double tol_factor = GlobalData::GetDouble("tol_factor");
    //        double tol = tol_factor * sqrt(grid->GetDx() * grid->GetDx() + grid->GetDy() * grid->GetDy() +
    //                                       grid->GetDz() * grid->GetDz());
    //        if (dist < tol)
    //        {
    //            m_cell_type[idx_proxy(i, j, k)] = PhysicalType::Solid;
    //            m_cell_type[idx_proxy(i, j, k - 1)] = PhysicalType::Solid;
    //            m_cell_type[idx_proxy(i, j - 1, k)] = PhysicalType::Solid;
    //            m_cell_type[idx_proxy(i, j - 1, k - 1)] = PhysicalType::Solid;
    //            m_cell_type[idx_proxy(i - 1, j, k)] = PhysicalType::Solid;
    //            m_cell_type[idx_proxy(i - 1, j, k - 1)] = PhysicalType::Solid;
    //            m_cell_type[idx_proxy(i - 1, j - 1, k)] = PhysicalType::Solid;
    //            m_cell_type[idx_proxy(i - 1, j - 1, k - 1)] = PhysicalType::Solid;
    //            total_error_num++;
    //        }
    //        else
    //        {
    //            index_type neighbor_i[6] = {i - 1, i + 1, i, i, i, i};
    //            index_type neighbor_j[6] = {j, j, j - 1, j + 1, j, j};
    //            index_type neighbor_k[6] = {k, k, k, k, k - 1, k + 1};
    //            double coord[6][3];
    //            int direction[6] = {1, 1, 1, 1, 1, 1};
    //            for (index_type iNeighbor = 0; iNeighbor < 6; iNeighbor++)
    //            {
    //                for (index_type iDim = 0; iDim < 3; iDim++)
    //                {
    //                    coord[iNeighbor][iDim] =
    //                        node->GetCoord(neighbor_i[iNeighbor], neighbor_j[iNeighbor], neighbor_k[iNeighbor])[iDim];
    //                }
    //            }
    //            for (index_type iNeighbor = 0; iNeighbor < 6; iNeighbor++)
    //            {
    //                if (m_node_type[idx_proxy(neighbor_i[iNeighbor], neighbor_j[iNeighbor], neighbor_k[iNeighbor])] ==
    //                    PhysicalType::Solid)
    //                {
    //                    direction[iNeighbor] = 0;
    //                    for (index_type iDim = 0; iDim < 3; iDim++)
    //                    {
    //                        coord[iNeighbor][iDim] = wall_coord[iDim];
    //                    }
    //                }
    //            }
    //            Eigen::Vector3d vec1, vec2, vec3;
    //            for (int iDim = 0; iDim < 3; iDim++)
    //            {
    //                vec1[iDim] = coord[1][iDim] - coord[0][iDim];
    //                vec2[iDim] = coord[3][iDim] - coord[2][iDim];
    //                vec3[iDim] = coord[5][iDim] - coord[4][iDim];
    //            }
    //            double coef_x[3], coef_y[3], coef_z[3];
    //            double volume = (vec1.cross(vec2)).dot(vec3);
    //            double cross[3];
    //            CrossProduct(vec1.data(), vec2.data(), cross);
    //            double volume2 = DotProduct(cross, vec3.data());
    //            if (abs(volume) < 1e-5 || std::isnan(fabs(volume)) || std::isinf(fabs(volume)))
    //            {
    //                m_cell_type[idx_proxy(i, j, k)] = PhysicalType::Solid;
    //                m_cell_type[idx_proxy(i, j, k - 1)] = PhysicalType::Solid;
    //                m_cell_type[idx_proxy(i, j - 1, k)] = PhysicalType::Solid;
    //                m_cell_type[idx_proxy(i, j - 1, k - 1)] = PhysicalType::Solid;
    //                m_cell_type[idx_proxy(i - 1, j, k)] = PhysicalType::Solid;
    //                m_cell_type[idx_proxy(i - 1, j, k - 1)] = PhysicalType::Solid;
    //                m_cell_type[idx_proxy(i - 1, j - 1, k)] = PhysicalType::Solid;
    //                m_cell_type[idx_proxy(i - 1, j - 1, k - 1)] = PhysicalType::Solid;
    //                total_error_num++;
    //                // Log::info("dist:{}, volume:{}, modify:{},{},{},{},{},{}", dist,
    //                // volume, direction[0], direction[1], direction[2], direction[3],
    //                // direction[4], direction[5]); Log::info("vec1:{},{},{}", vec1[0],
    //                // vec1[1], vec1[2]); Log::info("vec2:{},{},{}", vec2[0], vec2[1],
    //                // vec2[2]); Log::info("vec3:{},{},{}", vec3[0], vec3[1], vec3[2]);
    //                // Log::info("cross:{},{},{}", cross[0], cross[1], cross[2]);
    //                // Log::info("volume2:{}", volume2);
    //                // Log::info("coord0:{},{},{}", coord[0][0], coord[0][1],
    //                // coord[0][2]); Log::info("coord1:{},{},{}", coord[1][0],
    //                // coord[1][1], coord[1][2]); Log::info("coord2:{},{},{}",
    //                // coord[2][0], coord[2][1], coord[2][2]);
    //                // Log::info("coord3:{},{},{}", coord[3][0], coord[3][1],
    //                // coord[3][2]); Log::info("coord4:{},{},{}", coord[4][0],
    //                // coord[4][1], coord[4][2]); Log::info("coord5:{},{},{}",
    //                // coord[5][0], coord[5][1], coord[5][2]);
    //            }
    //        }
    //    }
    //}
    if (total_error_num > 0)
    {
        Log::warn("error trans node num:{}", total_error_num);
        return false;
    }
    return true;
}
void GridFNFactoryZaran::CheckTransFace()
{
    for (int iFace = 0; iFace < m_trans_face.size(); iFace++)
    {
        auto &face_node_idx = m_trans_face[iFace].idx_block;
        for (int iNode = 0; iNode < face_node_idx.size(); iNode++)
        {
            if (m_node_type[face_node_idx[iNode]] != PhysicalType::FluidSolid)
            {
                Log::error("The trans face is not correct");
                Log::error("node:{} is not fluid solid", face_node_idx[iNode]);
            }
        }
    }
}
void GridFNFactoryZaran::BuildWallNode()
{
    m_fn_info.node[m_layer_num].resize(m_fn_info.node[1].size());
    for (int iNode = 0; iNode < m_fn_info.node[1].size(); iNode++)
    {
        int i, j, k;
        auto trans_coord = m_fn_info.node[1][iNode].coord;
        double wall_coord[3];
        m_model_manager->GetClosestPoint(trans_coord, wall_coord);
        m_fn_info.node[m_layer_num][iNode].coord[0] = wall_coord[0];
        m_fn_info.node[m_layer_num][iNode].coord[1] = wall_coord[1];
        m_fn_info.node[m_layer_num][iNode].coord[2] = wall_coord[2];
        m_fn_info.node[m_layer_num][iNode].idx =
            iNode + m_fn_info.node[0].size() + m_fn_info.node[1].size() * (m_layer_num - 1);
    }
}

void GridFNFactoryZaran::BuildProjectNode()
{
    for (int iLayer = 2; iLayer < m_layer_num; iLayer++)
    {
        m_fn_info.node[iLayer].resize(m_fn_info.node[1].size());
    }
    for (int iNode = 0; iNode < m_fn_info.node[1].size(); iNode++)
    {
        int i, j, k;
        auto trans_coord = m_fn_info.node[1][iNode].coord;
        auto proj_coord = m_fn_info.node[m_layer_num][iNode].coord;
        for (int iLayer = 2; iLayer < m_layer_num; iLayer++)
        {
            double ratio = double(iLayer - 1) / (m_layer_num - 1);
            m_fn_info.node[iLayer][iNode].coord[0] = ratio * (proj_coord[0] - trans_coord[0]) + trans_coord[0];
            m_fn_info.node[iLayer][iNode].coord[1] = ratio * (proj_coord[1] - trans_coord[1]) + trans_coord[1];
            m_fn_info.node[iLayer][iNode].coord[2] = ratio * (proj_coord[2] - trans_coord[2]) + trans_coord[2];
            m_fn_info.node[iLayer][iNode].idx =
                iNode + m_fn_info.node[0].size() + m_fn_info.node[1].size() * (iLayer - 1);
        }
    }
}
void GridFNFactoryZaran::OptimizeWallNode()
{
    auto grid = GetFNGrid();
    if (grid->GetDim() == 2)
    {
        OptimizeWallNode2D();
    }
    else if (grid->GetDim() == 3)
    {
        OptimizeWallNode3D();
    }
    else
    {
        Log::error("Grid dimension is not supported: {}", grid->GetDim());
    }
}
void GridFNFactoryZaran::OptimizeWallNode2D()
{
    // 实现思路：对于每个过渡面元投影在物面的节点，构成的表面面元
    //  计算面元中心点到物面的距离
    //   如果距离小于一个阈值，则认为该表面面元是靠近物面的，
    //   如果距离较大，说明附近物面存在凹陷
    auto grid = GetBlockGrid();
    const double dx = grid->GetDx();
    const double dy = grid->GetDy();
    const double dz = grid->GetDz();
    const double tol_factor = GlobalData::GetDouble("tol_factor");
    double tol = tol_factor * sqrt(dx * dx + dy * dy) * 0.25;
    tol = 0.007;
    Log::info("tol: {}", tol);
    // 优化次数
    int optimize_time = 0;
    while (true)
    {
        int optimize_num = 0;
        for (int iFace = 0; iFace < m_trans_face.size(); iFace++)
        {
            // 二维面元为线段，获取两个节点坐标，平均
            auto &face_node_idx = m_trans_face[iFace].idx_slave;
            double coord[2][3] = {{0, 0, 0}, {0, 0, 0}};
            for (int iNode = 0; iNode < face_node_idx.size(); iNode++)
            {
                int idx = face_node_idx[iNode];
                coord[iNode][0] = m_fn_info.node[m_layer_num][idx].coord[0];
                coord[iNode][1] = m_fn_info.node[m_layer_num][idx].coord[1];
            }
            double face_coord[3] = {0, 0, 0};
            face_coord[0] = (coord[0][0] + coord[1][0]) / 2.0;
            face_coord[1] = (coord[0][1] + coord[1][1]) / 2.0;
            double wall_coord[3] = {0, 0, 0};
            m_model_manager->GetClosestPoint(face_coord, wall_coord);
            double dist = DistanceOfTwoPoints(face_coord, wall_coord);
            if (dist < tol)
                continue; // 如果距离小于阈值，则认为该面元是靠近物面的
                          // 如果距离较大，说明附近物面存在凹陷
                          // 取面心与节点连线的5等分点，向物面投影，该点作为新的物面节点

            coord[0][0] = (4 * coord[0][0] + 1 * face_coord[0]) / 5;
            coord[0][1] = (4 * coord[0][1] + 1 * face_coord[1]) / 5;
            double wall_coord1[3] = {0, 0, 0};
            m_model_manager->GetClosestPoint(coord[0], wall_coord1);
            double dist1 = DistanceOfTwoPoints(coord[0], wall_coord1);
            coord[1][0] = (4 * coord[1][0] + 1 * face_coord[0]) / 5;
            coord[1][1] = (4 * coord[1][1] + 1 * face_coord[1]) / 5;
            double wall_coord2[3] = {0, 0, 0};
            m_model_manager->GetClosestPoint(coord[1], wall_coord2);
            double dist2 = DistanceOfTwoPoints(coord[1], wall_coord2);
            if (dist1 > dist2)
            {
                m_fn_info.node[m_layer_num][face_node_idx[0]].coord[0] = wall_coord1[0];
                m_fn_info.node[m_layer_num][face_node_idx[0]].coord[1] = wall_coord1[1];
            }
            else
            {
                m_fn_info.node[m_layer_num][face_node_idx[1]].coord[0] = wall_coord2[0];
                m_fn_info.node[m_layer_num][face_node_idx[1]].coord[1] = wall_coord2[1];
            }
            optimize_num++;
        }
        if (optimize_num == 0)
        {
            Log::info("OptimizeWallNode2D: No more optimization needed.");
            break; // 如果没有优化，则退出循环
        }
        else if (optimize_time > 30)
        {
            Log::info("OptimizeWallNode2D: optimize_time>30.");
            break; 
        }
        else
        {
            Log::info("OptimizeWallNode2D: Optimized {} face.", optimize_num);
            optimize_time++;
        }
    }
}
void GridFNFactoryZaran::OptimizeWallNode3D()
{

}
void GridFNFactoryZaran::BuildNodeNeighbor()
{
    Log::info("BuildProjectNodeNeighbor");
    BuildProjectNodeNeighbor();
    Log::info("BuildTransNodeNeighbor");
    BuildTransNodeNeighbor();
    // Log::info("ReorderProjectNodeNeighbor");
    // ReorderProjectNodeNeighbor();
    Log::info("CheckProjectNodeNeighbor");
    CheckProjectNodeNeighbor();
}
void GridFNFactoryZaran::BuildProjectNodeNeighbor()
{
    auto grid = GetFNGrid();
    if (grid->GetDim() == 2)
    {
        BuildProjectNodeNeighbor2D();
    }
    else if (grid->GetDim() == 3)
    {
        BuildProjectNodeNeighbor3D();
    }
    else
    {
        Log::error("Grid dimension is not supported: {}", grid->GetDim());
    }
}

void GridFNFactoryZaran::BuildProjectNodeNeighbor2D()
{
    // direct neighbor
    dynamic_array<std::set<int>> node_neighbor_origin;
    node_neighbor_origin.resize(m_fn_info.node[1].size());
    for (int iFace = 0; iFace < m_trans_face.size(); iFace++)
    {
        // 二维面元为线段
        auto &face_node_idx = m_trans_face[iFace].idx_slave;
        node_neighbor_origin[face_node_idx[0]].insert(face_node_idx[1]);
        node_neighbor_origin[face_node_idx[1]].insert(face_node_idx[0]);
    }
    // extend neighbor, if the neighbor node is less than 3, extend the neighbor
    // node
    dynamic_array<std::set<int>> node_neighbor_extend;
    node_neighbor_extend.resize(m_fn_info.node[1].size());
    for (int iNode = 0; iNode < node_neighbor_origin.size(); iNode++)
    {

        node_neighbor_extend[iNode] = node_neighbor_origin[iNode];

        if (node_neighbor_origin[iNode].size() <= 1)
        {
            for (auto &neighbor : node_neighbor_origin[iNode])
            {
                for (auto &neighbor_neighbor : node_neighbor_origin[neighbor])
                {
                    node_neighbor_extend[iNode].insert(neighbor_neighbor);
                }
            }
            node_neighbor_extend[iNode].erase(iNode);
        }
    }
    // build the neighbor node
    for (int iNode = 0; iNode < m_fn_info.node[1].size(); iNode++)
    {
        for (int iLayer = 2; iLayer < m_layer_num; iLayer++)
        {
            int idx = m_fn_info.node[iLayer][iNode].idx;
            m_fn_info.node[iLayer][iNode].neighbor_node.resize(node_neighbor_extend[iNode].size() + 2);
            m_fn_info.node[iLayer][iNode].neighbor_node[0] = m_fn_info.node[iLayer - 1][iNode].idx;
            m_fn_info.node[iLayer][iNode].neighbor_node[1] = m_fn_info.node[iLayer + 1][iNode].idx;
            int id = 2;
            for (auto &neighbor : node_neighbor_extend[iNode])
            {
                m_fn_info.node[iLayer][iNode].neighbor_node[id++] = m_fn_info.node[iLayer][neighbor].idx;
            }
        }
        m_fn_info.node[m_layer_num][iNode].neighbor_node.resize(node_neighbor_extend[iNode].size() + 2);
        m_fn_info.node[m_layer_num][iNode].neighbor_node[0] = m_fn_info.node[m_layer_num - 1][iNode].idx;
        m_fn_info.node[m_layer_num][iNode].neighbor_node[1] = -1;
        int id = 2;
        for (auto &neighbor : node_neighbor_extend[iNode])
        {
            m_fn_info.node[m_layer_num][iNode].neighbor_node[id++] = m_fn_info.node[m_layer_num][neighbor].idx;
        }
    }
}

void GridFNFactoryZaran::BuildProjectNodeNeighbor3D()
{
    // direct neighbor
    dynamic_array<std::set<int>> node_neighbor_origin;
    node_neighbor_origin.resize(m_fn_info.node[1].size());
    for (int iFace = 0; iFace < m_trans_face.size(); iFace++)
    {
        // auto face_node_idx = m_trans_face[iFace].idx_master;
        // for (int iNode = 0; iNode < face_node_idx.size(); iNode++)
        // {
        //     for (auto &slave : m_trans_node)
        //     {
        //         if (slave.idx_block == face_node_idx[iNode])
        //         {
        //             face_node_idx[iNode] = slave.idx_local_layer;
        //             break;
        //         }
        //     }
        // }
        auto &face_node_idx = m_trans_face[iFace].idx_slave;
        int prior_node, next_node;
        for (int iNode = 0; iNode < face_node_idx.size(); iNode++)
        {
            int current_node = face_node_idx[iNode];
            prior_node = face_node_idx[(iNode + face_node_idx.size() - 1) % face_node_idx.size()];
            next_node = face_node_idx[(iNode + 1) % face_node_idx.size()];
            node_neighbor_origin[current_node].insert(prior_node);
            node_neighbor_origin[current_node].insert(next_node);
        }
    }
    // extend neighbor, if the neighbor node is less than 3, extend the neighbor
    // node
    dynamic_array<std::set<int>> node_neighbor_extend;
    node_neighbor_extend.resize(m_fn_info.node[1].size());
    for (int iNode = 0; iNode < node_neighbor_origin.size(); iNode++)
    {

        node_neighbor_extend[iNode] = node_neighbor_origin[iNode];

        if (node_neighbor_origin[iNode].size() <= 3)
        {
            for (auto &neighbor : node_neighbor_origin[iNode])
            {
                for (auto &neighbor_neighbor : node_neighbor_origin[neighbor])
                {
                    node_neighbor_extend[iNode].insert(neighbor_neighbor);
                }
            }
            node_neighbor_extend[iNode].erase(iNode);
        }
    }
    // build the neighbor node
    for (int iNode = 0; iNode < m_fn_info.node[1].size(); iNode++)
    {
        for (int iLayer = 2; iLayer < m_layer_num; iLayer++)
        {
            int idx = m_fn_info.node[iLayer][iNode].idx;
            m_fn_info.node[iLayer][iNode].neighbor_node.resize(node_neighbor_extend[iNode].size() + 2);
            m_fn_info.node[iLayer][iNode].neighbor_node[0] = m_fn_info.node[iLayer - 1][iNode].idx;
            m_fn_info.node[iLayer][iNode].neighbor_node[1] = m_fn_info.node[iLayer + 1][iNode].idx;
            int id = 2;
            for (auto &neighbor : node_neighbor_extend[iNode])
            {
                m_fn_info.node[iLayer][iNode].neighbor_node[id++] = m_fn_info.node[iLayer][neighbor].idx;
            }
        }
        m_fn_info.node[m_layer_num][iNode].neighbor_node.resize(node_neighbor_extend[iNode].size() + 2);
        m_fn_info.node[m_layer_num][iNode].neighbor_node[0] = m_fn_info.node[m_layer_num - 1][iNode].idx;
        m_fn_info.node[m_layer_num][iNode].neighbor_node[1] = -1;
        int id = 2;
        for (auto &neighbor : node_neighbor_extend[iNode])
        {
            m_fn_info.node[m_layer_num][iNode].neighbor_node[id++] = m_fn_info.node[m_layer_num][neighbor].idx;
        }
    }
}
void GridFNFactoryZaran::ReorderProjectNodeNeighbor()
{
    auto grid = GetFNGrid();
    auto node = grid->GetNode();
    struct node_pair
    {
        int node1, node2;
    };
    int node_num = m_fn_info.node[1].size();
    std::map<double, node_pair> node_pair_map;
    for (int iLayer = 2; iLayer < m_layer_num; iLayer++)
    {
        for (int iNode = 0; iNode < node_num; iNode++)
        {
            int idx = m_fn_info.node[iLayer][iNode].idx;
            auto &neighbor = m_fn_info.node[iLayer][iNode].neighbor_node;
            if (neighbor.size() == 6)
            {
                continue;
            }
            // {
            //     DVector3D vec1, vec2, vec3;
            //     for (int i = 0; i < 3; ++i)
            //     {
            //         vec1[i] = node->GetCoord(neighbor[0])[i] -
            //         node->GetCoord(neighbor[1])[i]; vec2[i] =
            //         node->GetCoord(neighbor[2])[i] -
            //         node->GetCoord(neighbor[3])[i]; vec3[i] =
            //         node->GetCoord(neighbor[4])[i] -
            //         node->GetCoord(neighbor[5])[i];
            //     }
            //     double volume = vec1.dot(vec2.cross(vec3));
            //     if (volume < 0)
            //     {
            //         std::swap(neighbor[2], neighbor[3]);
            //         std::swap(neighbor[4], neighbor[5]);
            //     }
            //     if (abs(volume) < EPSILON_NUMBER)
            //     {
            //         std::swap(neighbor[3], neighbor[4]);
            //     }
            //     for (int i = 0; i < 3; ++i)
            //     {
            //         vec1[i] = node->GetCoord(neighbor[1])[i] -
            //         node->GetCoord(neighbor[0])[i]; vec2[i] =
            //         node->GetCoord(neighbor[3])[i] -
            //         node->GetCoord(neighbor[2])[i]; vec3[i] =
            //         node->GetCoord(neighbor[5])[i] -
            //         node->GetCoord(neighbor[4])[i];
            //     }
            //     volume = vec1.dot(vec2.cross(vec3));
            //     Log::info("volume:{}", volume);
            // }

            node_pair main_pair;
            main_pair.node1 = neighbor[0];
            main_pair.node2 = neighbor[1];
            Eigen::Vector3d main_vec;
            for (int i = 0; i < 3; ++i)
            {
                main_vec[i] = node.GetCoord(main_pair.node2)[i] - node.GetCoord(main_pair.node1)[i];
            }
            main_vec.normalize();
            neighbor.erase(std::find(neighbor.begin(), neighbor.end(), main_pair.node1));
            neighbor.erase(std::find(neighbor.begin(), neighbor.end(), main_pair.node2));
            // 求出所有邻居节点在以主方向向量为法向量，经过当地节点的平面上的投影
            map<int, Eigen::Vector3d> node_proj_map;
            Eigen::Vector3d vec;
            for (int i = 0; i < neighbor.size(); ++i)
            {
                for (int k = 0; k < 3; ++k)
                {
                    vec[k] = node.GetCoord(neighbor[i])[k] - node.GetCoord(idx)[k];
                }
                vec -= vec.dot(main_vec) * main_vec;
                double vec_norm = vec.norm();
                // 如果投影向量的模长小�?1e-6，删除该邻居节点
                if (vec_norm < EPSILON_NUMBER)
                {
                    neighbor.erase(neighbor.begin() + i);
                    --i;
                    continue;
                }
                node_proj_map[neighbor[i]] = vec;
            }
            // 以第一个邻居节点投影向量为基准向量，求出基准向量以法向量为旋转轴旋转到其他投影向量的角�?,
            // 0~2pi
            map<double, int> node_angle_map;
            for (int i = 0; i < neighbor.size(); ++i)
            {
                if (i == 0)
                {
                    node_angle_map[0] = neighbor[i];
                    continue;
                }
                Eigen::Vector3d vec = node_proj_map[neighbor[i]];
                double angle = AngleOfTwoArray3D(node_proj_map[neighbor[0]].data(), vec.data());
                if ((node_proj_map[neighbor[0]].cross(vec).dot(main_vec) < 0))
                    angle = 2 * PI - angle;
                // 如果map中已经有这个角度，比较两个角度对应的邻居节点的距离，删除距离大的邻居节点
                if (node_angle_map.find(angle) != node_angle_map.end())
                {
                    if (node_proj_map[neighbor[i]].norm() > node_proj_map[node_angle_map[angle]].norm())
                    {
                        node_angle_map[angle] = neighbor[i];
                        neighbor.erase(std::find(neighbor.begin(), neighbor.end(), node_angle_map[angle]));
                    }
                    else
                    {
                        neighbor.erase(std::find(neighbor.begin(), neighbor.end(), neighbor[i]));
                    }
                    --i;
                    continue;
                }
                node_angle_map[angle] = neighbor[i];
            }
            // 根据角度排序后的邻居节点
            neighbor.clear();
            for (auto &i : node_angle_map)
            {
                neighbor.push_back(i.second);
            }
            // 求出邻居节点与当地节点之间的距离
            map<int, double> node_dis_map;
            for (int i = 0; i < neighbor.size(); ++i)
            {
                node_dis_map[neighbor[i]] = node_proj_map[neighbor[i]].norm();
            }
            node_pair_map.clear();
            // 获取下一个点的lamda表达�?
            auto get_next_node = [&](index_type iNode, dynamic_array<index_type> neiborNode) -> index_type {
                if (iNode == neiborNode.size() - 1)
                    return 0;
                else
                    return iNode + 1;
            };
            // 获取上一个点的lamda表达�?
            auto get_last_node = [&](index_type iNode, dynamic_array<index_type> neiborNode) -> index_type {
                if (iNode == 0)
                    return neiborNode.size() - 1;
                else
                    return iNode - 1;
            };
            for (int i = 0; i < neighbor.size(); ++i)
            {
                node_pair temp;
                temp.node1 = neighbor[i];
                temp.node2 = neighbor[get_next_node(i, neighbor)];
                Eigen::Vector3d vec1, vec2;
                vec1 = node_proj_map[temp.node1];
                vec2 = node_proj_map[temp.node2];
                double angle = AngleOfTwoArray3D(vec1.data(), vec2.data());
                if (vec1.cross(vec2).dot(main_vec) < 0)
                    angle = 2 * PI - angle;
                node_pair_map[angle] = temp;
            }

            if (node_pair_map.size() != neighbor.size())
                Log::info("node_pair_map.size()!=currentNeibor.size()");
            // 取出map中第一个点对，即夹角最小的点对
            // 删除这个点对中距离最大的�?
            while (neighbor.size() > 4)
            {

                auto &temp_pair = node_pair_map.begin()->second;
                index_type remove_node, remove_index;
                node_pair temp;
                if (node_dis_map[temp_pair.node1] > node_dis_map[temp_pair.node2])
                {
                    remove_node = temp_pair.node1;
                    remove_index = std::find(neighbor.begin(), neighbor.end(), remove_node) - neighbor.begin();
                    for (auto &i : node_pair_map)
                    {
                        if (i.second.node1 == neighbor[get_last_node(remove_index, neighbor)] &&
                            i.second.node2 == neighbor[remove_index])
                        {
                            node_pair_map.erase(i.first);
                            break;
                        }
                    }
                }
                else
                {
                    remove_node = temp_pair.node2;
                    remove_index = std::find(neighbor.begin(), neighbor.end(), remove_node) - neighbor.begin();
                    for (auto &i : node_pair_map)
                    {
                        if (i.second.node1 == neighbor[remove_index] &&
                            i.second.node2 == neighbor[get_next_node(remove_index, neighbor)])
                        {
                            node_pair_map.erase(i.first);
                            break;
                        }
                    }
                }
                node_pair_map.erase(node_pair_map.begin());
                temp.node1 = neighbor[get_last_node(remove_index, neighbor)];
                temp.node2 = neighbor[get_next_node(remove_index, neighbor)];
                neighbor.erase(std::find(neighbor.begin(), neighbor.end(), remove_node));
                node_dis_map.erase(remove_node);
                Eigen::Vector3d vec1, vec2;
                vec1 = node_proj_map[temp.node1];
                vec2 = node_proj_map[temp.node2];
                double angle = AngleOfTwoArray3D(vec1.data(), vec2.data()) + GetRand(0.0, 1.0) * EPSILON_NUMBER;
                if (vec1.cross(vec2).dot(main_vec) < 0)
                    angle = 2 * PI - angle;
                node_pair_map[angle] = temp;
            }
            //! 前四个点没有排序
            neighbor.push_back(main_pair.node1);
            neighbor.push_back(main_pair.node2);
            std::swap(neighbor[0], neighbor[4]);
            std::swap(neighbor[1], neighbor[5]);
        }
    }
}
void GridFNFactoryZaran::CheckProjectNodeNeighbor()
{
    auto grid = GetBlockGrid();
    if (grid->GetDim() == 2)
    {
        return CheckProjectNodeNeighbor2D();
    }
    else if (grid->GetDim() == 3)
    {
        return CheckProjectNodeNeighbor3D();
    }
    else
    {
        Log::error("Grid dimension is not supported: {}", grid->GetDim());
    }
}
void GridFNFactoryZaran::CheckProjectNodeNeighbor2D()
{
    auto grid = GetFNGrid();
    NodeFN &node = grid->GetNode();
    int node_num = m_fn_info.node[1].size();
    double delta = 15 * PI / 180;
    double angle;
    for (int iLayer = 1; iLayer < m_layer_num; iLayer++)
    {
        for (int iNode = 0; iNode < node_num; iNode++)
        {
            int idx = m_fn_info.node[iLayer][iNode].idx;
            auto &neighbor = m_fn_info.node[iLayer][iNode].neighbor_node;
            dynamic_array<Eigen::Vector3d> vec(3);
            for (int i = 0; i < 3; i++)
            {
                vec[0][i] = node.GetCoord(neighbor[1])[i] - node.GetCoord(neighbor[0])[i];
                vec[1][i] = node.GetCoord(neighbor[3])[i] - node.GetCoord(neighbor[2])[i];
            }
            vec[2] = {0, 0, 1};
            double volume = vec[0].cross(vec[1]).dot(vec[2]);
            angle = AngleOfTwoArray3D(vec[1].data(), vec[2].data());
            // 检查是否是右手坐标系
            if (volume < 0)
            {
                std::swap(neighbor[0], neighbor[1]);
            }
            angle = AngleOfTwoArray3D(vec[1].data(), vec[2].data());
            // 检查是否是右手坐标系
            for (int i = 0; i < 3; i++)
            {
                vec[0][i] = node.GetCoord(neighbor[1])[i] - node.GetCoord(neighbor[0])[i];
                vec[1][i] = node.GetCoord(neighbor[3])[i] - node.GetCoord(neighbor[2])[i];
            }
            if (vec[0].cross(vec[1]).dot(vec[2]) < 0)
            {
                std::swap(neighbor[0], neighbor[1]);
            }
            // 检查是否是右手坐标系
            for (int i = 0; i < 3; i++)
            {
                vec[0][i] = node.GetCoord(neighbor[1])[i] - node.GetCoord(neighbor[0])[i];
                vec[1][i] = node.GetCoord(neighbor[3])[i] - node.GetCoord(neighbor[2])[i];
            }
            volume = vec[0].cross(vec[1]).dot(vec[2]);
            if (volume < 0 || std::isnan(abs(volume)) || std::isinf(abs(volume)))
            {
                Log::error("The node:{} is not right hand", idx);
            }
        }
    }
}

void GridFNFactoryZaran::CheckProjectNodeNeighbor3D()
{
    auto grid = GetFNGrid();
    NodeFN &node = grid->GetNode();
    int node_num = m_fn_info.node[1].size();
    double delta = 15 * PI / 180;
    double angle;
    for (int iLayer = 1; iLayer < m_layer_num; iLayer++)
    {
        for (int iNode = 0; iNode < node_num; iNode++)
        {
            int idx = m_fn_info.node[iLayer][iNode].idx;
            auto &neighbor = m_fn_info.node[iLayer][iNode].neighbor_node;
            dynamic_array<Eigen::Vector3d> vec(3);
            for (int i = 0; i < 3; i++)
            {
                vec[0][i] = node.GetCoord(neighbor[1])[i] - node.GetCoord(neighbor[0])[i];
                vec[1][i] = node.GetCoord(neighbor[3])[i] - node.GetCoord(neighbor[2])[i];
                vec[2][i] = node.GetCoord(neighbor[5])[i] - node.GetCoord(neighbor[4])[i];
            }
            double volume = vec[0].cross(vec[1]).dot(vec[2]);
            angle = AngleOfTwoArray3D(vec[1].data(), vec[2].data());
            // 检查是否是右手坐标系
            if (volume < 0)
            {
                std::swap(neighbor[0], neighbor[1]);
            }
            angle = AngleOfTwoArray3D(vec[1].data(), vec[2].data());
            // j,k方向平行
            if (abs(angle) < delta)
            {
                std::swap(neighbor[3], neighbor[5]);
            }
            else if (abs(angle - PI) < delta)
            {
                std::swap(neighbor[3], neighbor[4]);
            }
            // double angle = AngleOfTwoArray3D(vec[0].data(), vec[1].data());
            // // i,j方向平行
            // if (abs(angle) < delta)
            // {
            //     std::swap(neighbor[1], neighbor[3]);
            // }
            // else if (abs(angle - PI) < delta)
            // {
            //     std::swap(neighbor[1], neighbor[2]);
            // }

            // angle = AngleOfTwoArray3D(vec[0].data(), vec[2].data());
            // // i,k方向平行
            // if (abs(angle) < delta)
            // {
            //     std::swap(neighbor[1], neighbor[5]);
            // }
            // else if (abs(angle - PI) < delta)
            // {
            //     std::swap(neighbor[1], neighbor[4]);
            // }

            // 检查是否是右手坐标系
            for (int i = 0; i < 3; i++)
            {
                vec[0][i] = node.GetCoord(neighbor[1])[i] - node.GetCoord(neighbor[0])[i];
                vec[1][i] = node.GetCoord(neighbor[3])[i] - node.GetCoord(neighbor[2])[i];
                vec[2][i] = node.GetCoord(neighbor[5])[i] - node.GetCoord(neighbor[4])[i];
            }
            if (vec[0].cross(vec[1]).dot(vec[2]) < 0)
            {
                std::swap(neighbor[0], neighbor[1]);
            }
            // 检查是否是右手坐标系
            for (int i = 0; i < 3; i++)
            {
                vec[0][i] = node.GetCoord(neighbor[1])[i] - node.GetCoord(neighbor[0])[i];
                vec[1][i] = node.GetCoord(neighbor[3])[i] - node.GetCoord(neighbor[2])[i];
                vec[2][i] = node.GetCoord(neighbor[5])[i] - node.GetCoord(neighbor[4])[i];
            }
            volume = vec[0].cross(vec[1]).dot(vec[2]);
            if (volume < 0 || std::isnan(abs(volume)) || std::isinf(abs(volume)))
            {
                Log::error("The node:{} is not right hand", idx);
            }
        }
    }
}
void GridFNFactoryZaran::BuildTransNodeNeighbor()
{
    auto grid = GetBlockGrid();
    if (grid->GetDim() == 2)
    {
        return BuildTransNodeNeighbor2D();
    }
    else if (grid->GetDim() == 3)
    {
        return BuildTransNodeNeighbor3D();
    }
    else
    {
        Log::error("Grid dimension is not supported: {}", grid->GetDim());
    }
}
void GridFNFactoryZaran::BuildTransNodeNeighbor2D()
{
    IdProxyStruct &idx_proxy = *m_idx_proxy;
    for (int iNode = 0; iNode < m_fn_info.node[1].size(); iNode++)
    {
        m_fn_info.node[1][iNode].neighbor_node.resize(4);
    }
    for (auto &nodes : m_trans_node)
    {
        index_type i, j, k;
        idx_proxy.GetIdxStruct(nodes.idx_block, i, j, k);
        m_fn_info.node[1][nodes.idx_local_layer].neighbor_node[0] = idx_proxy(i - 1, j, k);
        m_fn_info.node[1][nodes.idx_local_layer].neighbor_node[1] = idx_proxy(i + 1, j, k);
        m_fn_info.node[1][nodes.idx_local_layer].neighbor_node[2] = idx_proxy(i, j - 1, k);
        m_fn_info.node[1][nodes.idx_local_layer].neighbor_node[3] = idx_proxy(i, j + 1, k);
    }
    for (int iNode = 0; iNode < m_fn_info.node[1].size(); iNode++)
    {
        int next_layer_node = m_fn_info.node[2][iNode].idx;
        for (int iNeighbor = 0; iNeighbor < 4; iNeighbor++)
        {
            bool find = false;
            int solid_num = 0;
            int idx_master = m_fn_info.node[1][iNode].neighbor_node[iNeighbor];
            if (m_node_type[idx_master] == PhysicalType::FluidSolid)
            {
                for (auto &nodes : m_trans_node)
                {
                    if (nodes.idx_block == idx_master)
                    {
                        m_fn_info.node[1][iNode].neighbor_node[iNeighbor] =
                            m_fn_info.node[1][nodes.idx_local_layer].idx;
                        find = true;
                        break;
                    }
                }
            }
            else if (m_node_type[idx_master] == PhysicalType::Fluid)
            {
                for (auto &nodes : m_ref_node)
                {
                    if (nodes.idx_block == idx_master)
                    {
                        m_fn_info.node[1][iNode].neighbor_node[iNeighbor] =
                            m_fn_info.node[0][nodes.idx_local_layer].idx;
                        find = true;
                        break;
                    }
                }
            }
            else if (m_node_type[idx_master] == PhysicalType::Solid)
            {
                m_fn_info.node[1][iNode].neighbor_node[iNeighbor] = next_layer_node;
                find = true;
                solid_num++;
            }
            if (!find)
            {
                index_type i, j, k;
                idx_proxy.GetIdxStruct(idx_master, i, j, k);
                Log::error("Invalid node:{}, {},{},{}, type:{}", idx_master, i, j, k, int(m_node_type[idx_master]));
            }
            if (solid_num > 3)
            {
                Log::error("Invalid node:{}, type:{}", idx_master, int(m_node_type[idx_master]));
            }
        }
    }
}

void GridFNFactoryZaran::BuildTransNodeNeighbor3D()
{
    IdProxyStruct &idx_proxy = *m_idx_proxy;
    for (int iNode = 0; iNode < m_fn_info.node[1].size(); iNode++)
    {
        m_fn_info.node[1][iNode].neighbor_node.resize(6);
    }
    for (auto &nodes : m_trans_node)
    {
        index_type i, j, k;
        idx_proxy.GetIdxStruct(nodes.idx_block, i, j, k);
        m_fn_info.node[1][nodes.idx_local_layer].neighbor_node[0] = idx_proxy(i - 1, j, k);
        m_fn_info.node[1][nodes.idx_local_layer].neighbor_node[1] = idx_proxy(i + 1, j, k);
        m_fn_info.node[1][nodes.idx_local_layer].neighbor_node[2] = idx_proxy(i, j - 1, k);
        m_fn_info.node[1][nodes.idx_local_layer].neighbor_node[3] = idx_proxy(i, j + 1, k);
        m_fn_info.node[1][nodes.idx_local_layer].neighbor_node[4] = idx_proxy(i, j, k - 1);
        m_fn_info.node[1][nodes.idx_local_layer].neighbor_node[5] = idx_proxy(i, j, k + 1);
    }
    for (int iNode = 0; iNode < m_fn_info.node[1].size(); iNode++)
    {
        int next_layer_node = m_fn_info.node[2][iNode].idx;
        for (int iNeighbor = 0; iNeighbor < 6; iNeighbor++)
        {
            bool find = false;
            int solid_num = 0;
            int idx_master = m_fn_info.node[1][iNode].neighbor_node[iNeighbor];
            if (m_node_type[idx_master] == PhysicalType::FluidSolid)
            {
                for (auto &nodes : m_trans_node)
                {
                    if (nodes.idx_block == idx_master)
                    {
                        m_fn_info.node[1][iNode].neighbor_node[iNeighbor] =
                            m_fn_info.node[1][nodes.idx_local_layer].idx;
                        find = true;
                        break;
                    }
                }
            }
            else if (m_node_type[idx_master] == PhysicalType::Fluid)
            {
                for (auto &nodes : m_ref_node)
                {
                    if (nodes.idx_block == idx_master)
                    {
                        m_fn_info.node[1][iNode].neighbor_node[iNeighbor] =
                            m_fn_info.node[0][nodes.idx_local_layer].idx;
                        find = true;
                        break;
                    }
                }
            }
            else if (m_node_type[idx_master] == PhysicalType::Solid)
            {
                m_fn_info.node[1][iNode].neighbor_node[iNeighbor] = next_layer_node;
                find = true;
                solid_num++;
            }
            if (!find)
            {
                index_type i, j, k;
                idx_proxy.GetIdxStruct(idx_master, i, j, k);
                Log::error("Invalid node:{}, {},{},{}, type:{}", idx_master, i, j, k, int(m_node_type[idx_master]));
            }
            if (solid_num > 3)
            {
                Log::error("Invalid node:{}, type:{}", idx_master, int(m_node_type[idx_master]));
            }
        }
    }
}
void GridFNFactoryZaran::BuildFNCell()
{
    auto grid = GetBlockGrid();
    if (grid->GetDim() == 2)
    {
        BuildFNCell2D();
    }
    else if (grid->GetDim() == 3)
    {
        BuildFNCell3D();
    }
    else
    {
        Log::error("Grid dimension is not supported: {}", grid->GetDim());
    }
}
void GridFNFactoryZaran::BuildFNCell2D()
{
    m_fn_info.cell.resize(m_trans_face.size() * (m_layer_num - 1));
    for (int iFace = 0; iFace < m_trans_face.size(); iFace++)
    {
        auto &face_node_idx = m_trans_face[iFace].idx_slave;
        for (int iLevel = 1; iLevel < m_layer_num; iLevel++)
        {
            int iCell = iFace * (m_layer_num - 1) + iLevel - 1;
            m_fn_info.cell[iCell].resize(4);
            m_fn_info.cell[iCell][0] = m_fn_info.node[iLevel][face_node_idx[0]].idx;
            m_fn_info.cell[iCell][1] = m_fn_info.node[iLevel][face_node_idx[1]].idx;
            m_fn_info.cell[iCell][2] = m_fn_info.node[iLevel + 1][face_node_idx[1]].idx;
            m_fn_info.cell[iCell][3] = m_fn_info.node[iLevel + 1][face_node_idx[0]].idx;
        }
    }
}
void GridFNFactoryZaran::BuildFNCell3D()
{
    m_fn_info.cell.resize(m_trans_face.size() * (m_layer_num - 1));
    for (int iFace = 0; iFace < m_trans_face.size(); iFace++)
    {
        // auto face_node_idx = m_trans_face[iFace].idx_master;
        // for (int iNode = 0; iNode < face_node_idx.size(); iNode++)
        // {
        //     for (auto &slave : m_trans_node)
        //     {
        //         if (slave.idx_block == face_node_idx[iNode])
        //         {
        //             face_node_idx[iNode] = slave.idx_local_layer;
        //             break;
        //         }
        //     }
        // }
        auto &face_node_idx = m_trans_face[iFace].idx_slave;
        for (int iLevel = 1; iLevel < m_layer_num; iLevel++)
        {
            int iCell = iFace * (m_layer_num - 1) + iLevel - 1;
            m_fn_info.cell[iCell].resize(8);
            for (int iNode = 0; iNode < face_node_idx.size(); iNode++)
            {
                m_fn_info.cell[iCell][iNode] = m_fn_info.node[iLevel][face_node_idx[iNode]].idx;
            }
            for (int iNode = 0; iNode < face_node_idx.size(); iNode++)
            {
                m_fn_info.cell[iCell][iNode + 4] = m_fn_info.node[iLevel + 1][face_node_idx[iNode]].idx;
            }
        }
    }
}
} // namespace zaran