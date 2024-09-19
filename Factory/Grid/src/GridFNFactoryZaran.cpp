#include "GridFNFactoryZaran.h"
#include "Log.h"
#include "MathBasic.h"
#include <omp.h>
namespace zaran
{
    void GridFNFactoryZaran::CreateGrid(GridBlock *grid_master, GridFN *&grid, ModelManager *model_manager)
    {
        m_model_manager = model_manager;
        m_block_grid = grid_master;
        m_fn_grid = grid;
        m_idx_proxy = new StructIdxProxy(m_block_grid);
        m_layer_num = 2;
        TagBlockGrid();
        BuildFNNodeCoord();
        SetFNGridNode();
        WriteProjectNode();
        WriteModelSurface();
        WriteTransFace();
        BuildCell();
        SetFNGridCell();
        WriteSlaveGrid();
        BuildNodeNeighbor();
        SetFNGridNodeNeighbor();
        SetFNGridBoundary();
        SetFNGridBoundaryFace();
    }

    void GridFNFactoryZaran::TagBlockGrid()
    {
        TagCells();
        TagNodes();
        CheckTransFace();
        WriteNodeTag();
        SetNodeTag();
    }

    void GridFNFactoryZaran::TagCells()
    {
        auto grid = GetBlockGrid();
        int ni, nj, nk;
        ni = grid->GetNi();
        nj = grid->GetNj();
        nk = grid->GetNk();
        int ghost_size = grid->GetGhostLevel();
        m_cell_type.resize(ni * nj * nk);
        for (int iCell = 0; iCell < ni * nj * nk; iCell++)
        {
            m_cell_type[iCell] = PhysicalType::Unset;
        }
        double dx = grid->GetDx();
        double dy = grid->GetDy();
        double dz = grid->GetDz();
        double tol = 0.5 * sqrt(dx * dx + dy * dy + dz * dz);
        auto model_manager = GetModelManager();
        auto box = grid->GetBoundBox();
        auto model_box = model_manager->GetBox();
        auto idx_proxy = GetIdxProxy();
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
        //             if (cell_center[0] < model_box.x_min - tol || cell_center[0] > model_box.x_max + tol || cell_center[1] < model_box.y_min - tol || cell_center[1] > model_box.y_max + tol || cell_center[2] < model_box.z_min - tol || cell_center[2] > model_box.z_max + tol)
        //             {
        //                 m_cell_type[iCell] = PhysicalType::Fluid;
        //             }
        //             else
        //             {
        //                 double dist = model_manager->GetClosestDistance(cell_center);
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
        m_cell_type[m_idx_proxy->GetIdx(0, 0, 0)] = PhysicalType::Fluid;
        auto IsValidCell = [&](int i, int j, int k) -> bool
        {
            return i >= 0 && i < ni && j >= 0 && j < nj && k >= 0 && k < nk;
        };
        for (int k = 0; k < nk; k++)
        {
            for (int j = 0; j < nj; j++)
            {
                for (int i = 0; i < ni; i++)
                {
                    int idx = m_idx_proxy->GetIdx(i, j, k);
                    if (m_cell_type[idx] == PhysicalType::Fluid)
                    {
                        if (IsValidCell(i + 1, j, k) && m_cell_type[m_idx_proxy->GetIdx(i + 1, j, k)] == PhysicalType::Unset)
                        {
                            m_cell_type[m_idx_proxy->GetIdx(i + 1, j, k)] = PhysicalType::Fluid;
                        }
                        if (IsValidCell(i, j + 1, k) && m_cell_type[m_idx_proxy->GetIdx(i, j + 1, k)] == PhysicalType::Unset)
                        {
                            m_cell_type[m_idx_proxy->GetIdx(i, j + 1, k)] = PhysicalType::Fluid;
                        }
                        if (IsValidCell(i, j, k + 1) && m_cell_type[m_idx_proxy->GetIdx(i, j, k + 1)] == PhysicalType::Unset)
                        {
                            m_cell_type[m_idx_proxy->GetIdx(i, j, k + 1)] = PhysicalType::Fluid;
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
                    int idx = m_idx_proxy->GetIdx(i, j, k);
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
                        int idx = m_idx_proxy->GetIdx(i, j, k);
                        if (m_cell_type[idx] != PhysicalType::Fluid)
                            continue;
                        if (m_cell_type[m_idx_proxy->GetIdx(i + 1, j, k)] == PhysicalType::Solid && m_cell_type[m_idx_proxy->GetIdx(i - 1, j, k)] == PhysicalType::Solid)
                        {
                            m_cell_type[idx] = PhysicalType::Solid;
                            new_solid_num++;
                        }
                        else if (m_cell_type[m_idx_proxy->GetIdx(i, j + 1, k)] == PhysicalType::Solid && m_cell_type[m_idx_proxy->GetIdx(i, j - 1, k)] == PhysicalType::Solid)
                        {
                            m_cell_type[idx] = PhysicalType::Solid;
                            new_solid_num++;
                        }
                        else if (m_cell_type[m_idx_proxy->GetIdx(i, j, k + 1)] == PhysicalType::Solid && m_cell_type[m_idx_proxy->GetIdx(i, j, k - 1)] == PhysicalType::Solid)
                        {
                            m_cell_type[idx] = PhysicalType::Solid;
                            new_solid_num++;
                        }
                    }
                }
            }
            Log::info("new_solid_num={}", new_solid_num);
        }

        // tag the fluid-solid cells
        for (int k = 1; k < nk - 1; k++)
        {
            for (int j = 1; j < nj - 1; j++)
            {
                for (int i = 1; i < ni - 1; i++)
                {
                    int idx = m_idx_proxy->GetIdx(i, j, k);
                    if (m_cell_type[idx] != PhysicalType::Fluid)
                        continue;
                    if (m_cell_type[m_idx_proxy->GetIdx(i + 1, j, k)] == PhysicalType::Solid || m_cell_type[m_idx_proxy->GetIdx(i + 1, j, k)] == PhysicalType::Unset)
                    {
                        m_cell_type[idx] = PhysicalType::FluidSolid;
                    }
                    if (m_cell_type[m_idx_proxy->GetIdx(i - 1, j, k)] == PhysicalType::Solid || m_cell_type[m_idx_proxy->GetIdx(i - 1, j, k)] == PhysicalType::Unset)
                    {
                        m_cell_type[idx] = PhysicalType::FluidSolid;
                    }
                    if (m_cell_type[m_idx_proxy->GetIdx(i, j + 1, k)] == PhysicalType::Solid || m_cell_type[m_idx_proxy->GetIdx(i, j + 1, k)] == PhysicalType::Unset)
                    {
                        m_cell_type[idx] = PhysicalType::FluidSolid;
                    }
                    if (m_cell_type[m_idx_proxy->GetIdx(i, j - 1, k)] == PhysicalType::Solid || m_cell_type[m_idx_proxy->GetIdx(i, j - 1, k)] == PhysicalType::Unset)
                    {
                        m_cell_type[idx] = PhysicalType::FluidSolid;
                    }
                    if (m_cell_type[m_idx_proxy->GetIdx(i, j, k + 1)] == PhysicalType::Solid || m_cell_type[m_idx_proxy->GetIdx(i, j, k + 1)] == PhysicalType::Unset)
                    {
                        m_cell_type[idx] = PhysicalType::FluidSolid;
                    }
                    if (m_cell_type[m_idx_proxy->GetIdx(i, j, k - 1)] == PhysicalType::Solid || m_cell_type[m_idx_proxy->GetIdx(i, j, k - 1)] == PhysicalType::Unset)
                    {
                        m_cell_type[idx] = PhysicalType::FluidSolid;
                    }
                }
            }
        }
    }
    void GridFNFactoryZaran::TagNodes()
    {
        auto grid = GetBlockGrid();
        auto node = grid->GetNode();
        int ni, nj, nk;
        ni = grid->GetNi();
        nj = grid->GetNj();
        nk = grid->GetNk();
        auto IsValidNode = [&](int i, int j, int k) -> bool
        {
            return i >= 0 && i < ni && j >= 0 && j < nj && k >= 0 && k < nk;
        };
        auto IsValidCell = [&](int i, int j, int k) -> bool
        {
            return i >= 0 && i < ni - 1 && j >= 0 && j < nj - 1 && k >= 0 && k < nk - 1;
        };
        m_node_type.resize(ni * nj * nk);
        for (int idx = 0; idx < ni * nj * nk; idx++)
        {
            m_node_type[idx] = PhysicalType::Unset;
        }
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        std::set<int> trans_node_set;
        std::set<TransFace> trans_face_set;
        // find the fluid-solid nodes
        for (int k = 0; k < nk; k++)
        {
            for (int j = 0; j < nj; j++)
            {
                for (int i = 0; i < ni; i++)
                {
                    int idx = m_idx_proxy->GetIdx(i, j, k);
                    if (IsValidCell(i - 1, j, k))
                    {
                        int idx_left = m_idx_proxy->GetIdx(i - 1, j, k);
                        if (m_cell_type[idx_left] == PhysicalType::FluidSolid && m_cell_type[idx] == PhysicalType::Fluid || m_cell_type[idx_left] == PhysicalType::Fluid && m_cell_type[idx] == PhysicalType::FluidSolid)
                        {
                            trans_face_set.insert(TransFace{{m_idx_proxy->GetIdx(i, j, k), m_idx_proxy->GetIdx(i, j + 1, k), m_idx_proxy->GetIdx(i, j + 1, k + 1), m_idx_proxy->GetIdx(i, j, k + 1)}});
                            m_node_type[m_idx_proxy->GetIdx(i, j, k)] = PhysicalType::FluidSolid;
                            m_node_type[m_idx_proxy->GetIdx(i, j + 1, k)] = PhysicalType::FluidSolid;
                            m_node_type[m_idx_proxy->GetIdx(i, j, k + 1)] = PhysicalType::FluidSolid;
                            m_node_type[m_idx_proxy->GetIdx(i, j + 1, k + 1)] = PhysicalType::FluidSolid;
                        }
                    }
                    if (IsValidCell(i, j - 1, k))
                    {
                        int idx_down = m_idx_proxy->GetIdx(i, j - 1, k);
                        if (m_cell_type[idx_down] == PhysicalType::FluidSolid && m_cell_type[idx] == PhysicalType::Fluid || m_cell_type[idx_down] == PhysicalType::Fluid && m_cell_type[idx] == PhysicalType::FluidSolid)
                        {
                            trans_face_set.insert(TransFace{{m_idx_proxy->GetIdx(i, j, k), m_idx_proxy->GetIdx(i + 1, j, k), m_idx_proxy->GetIdx(i + 1, j, k + 1), m_idx_proxy->GetIdx(i, j, k + 1)}});
                            m_node_type[m_idx_proxy->GetIdx(i, j, k)] = PhysicalType::FluidSolid;
                            m_node_type[m_idx_proxy->GetIdx(i + 1, j, k)] = PhysicalType::FluidSolid;
                            m_node_type[m_idx_proxy->GetIdx(i, j, k + 1)] = PhysicalType::FluidSolid;
                            m_node_type[m_idx_proxy->GetIdx(i + 1, j, k + 1)] = PhysicalType::FluidSolid;
                        }
                    }
                    if (IsValidCell(i, j, k - 1))
                    {
                        int idx_back = m_idx_proxy->GetIdx(i, j, k - 1);
                        if (m_cell_type[idx_back] == PhysicalType::FluidSolid && m_cell_type[idx] == PhysicalType::Fluid || m_cell_type[idx_back] == PhysicalType::Fluid && m_cell_type[idx] == PhysicalType::FluidSolid)
                        {
                            trans_face_set.insert(TransFace{{m_idx_proxy->GetIdx(i, j, k), m_idx_proxy->GetIdx(i + 1, j, k), m_idx_proxy->GetIdx(i + 1, j + 1, k), m_idx_proxy->GetIdx(i, j + 1, k)}});
                            m_node_type[m_idx_proxy->GetIdx(i, j, k)] = PhysicalType::FluidSolid;
                            m_node_type[m_idx_proxy->GetIdx(i + 1, j, k)] = PhysicalType::FluidSolid;
                            m_node_type[m_idx_proxy->GetIdx(i, j + 1, k)] = PhysicalType::FluidSolid;
                            m_node_type[m_idx_proxy->GetIdx(i + 1, j + 1, k)] = PhysicalType::FluidSolid;
                        }
                    }
                }
            }
        }
        m_trans_face.resize(trans_face_set.size());
        int idx = 0;
        for (auto &trans : trans_face_set)
        {
            m_trans_face[idx++].idx_master = trans.idx_master;
        }

        m_node_type[m_idx_proxy->GetIdx(0, 0, 0)] = PhysicalType::Fluid;
        // find the fluid nodes
        for (int k = 0; k < nk; k++)
        {
            for (int j = 0; j < nj; j++)
            {
                for (int i = 0; i < ni; i++)
                {
                    int idx = m_idx_proxy->GetIdx(i, j, k);
                    if (m_node_type[idx] != PhysicalType::Fluid)
                        continue;
                    if (IsValidNode(i + 1, j, k))
                    {
                        if (m_node_type[m_idx_proxy->GetIdx(i + 1, j, k)] == PhysicalType::Unset)
                            m_node_type[m_idx_proxy->GetIdx(i + 1, j, k)] = PhysicalType::Fluid;
                    }
                    if (IsValidNode(i, j + 1, k))
                    {
                        if (m_node_type[m_idx_proxy->GetIdx(i, j + 1, k)] == PhysicalType::Unset)
                            m_node_type[m_idx_proxy->GetIdx(i, j + 1, k)] = PhysicalType::Fluid;
                    }
                    if (IsValidNode(i, j, k + 1))
                    {
                        if (m_node_type[m_idx_proxy->GetIdx(i, j, k + 1)] == PhysicalType::Unset)
                            m_node_type[m_idx_proxy->GetIdx(i, j, k + 1)] = PhysicalType::Fluid;
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

    void GridFNFactoryZaran::ProcessCell(int start_i, int end_i, int start_j, int end_j, int start_k, int end_k)
    {
        // Log::info("start_i={}, end_i={}, start_j={}, end_j={}, start_k={}, end_k={}", start_i, end_i, start_j, end_j, start_k, end_k);
        auto grid = GetBlockGrid();
        auto &grid_box = grid->GetBoundBox();
        auto model_manager = GetModelManager();
        auto model_box = model_manager->GetBox();
        int mid_i = (start_i + end_i) / 2;
        int mid_j = (start_j + end_j) / 2;
        int mid_k = (start_k + end_k) / 2;
        double x_min, x_max, y_min, y_max, z_min, z_max;
        x_min = grid_box.x_min + (start_i - grid->GetGhostLevel() + 0.5) * grid->GetDx();
        x_max = grid_box.x_min + (end_i - 1 - grid->GetGhostLevel() + 0.5) * grid->GetDx();
        y_min = grid_box.y_min + (start_j - grid->GetGhostLevel() + 0.5) * grid->GetDy();
        y_max = grid_box.y_min + (end_j - 1 - grid->GetGhostLevel() + 0.5) * grid->GetDy();
        z_min = grid_box.z_min + (start_k - grid->GetGhostLevel() + 0.5) * grid->GetDz();
        z_max = grid_box.z_min + (end_k - 1 - grid->GetGhostLevel() + 0.5) * grid->GetDz();
        double cell_center[3];
        cell_center[0] = 0.5 * (x_min + x_max);
        cell_center[1] = 0.5 * (y_min + y_max);
        cell_center[2] = 0.5 * (z_min + z_max);
        double tol_x = x_max - x_min;
        double tol_y = y_max - y_min;
        double tol_z = z_max - z_min;
        double tol = 0.5 * sqrt(tol_x * tol_x + tol_y * tol_y + tol_z * tol_z);
        double tol1 = 0.5 * sqrt(grid->GetDx() * grid->GetDx() + grid->GetDy() * grid->GetDy() + grid->GetDz() * grid->GetDz());
        tol += tol1;
        PhysicalType cell_type = PhysicalType::Unset;
        double dist = 0;
        if (x_max < model_box.x_min - tol1 || x_min > model_box.x_max + tol1 || y_max < model_box.y_min - tol1 || y_min > model_box.y_max + tol1 || z_max < model_box.z_min - tol1 || z_min > model_box.z_max + tol1)
        {
            cell_type = PhysicalType::Unset;
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
        if (end_i - start_i <= 1 && end_j - start_j <= 1 && end_k - start_k <= 1)
        {
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
                    int idx = m_idx_proxy->GetIdx(start_i, start_j, start_k);
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

    void GridFNFactoryZaran::SetNodeTag()
    {
        auto grid = GetBlockGrid();
        auto node = grid->GetNode();
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        for (int k = ks; k <= ke; k++)
        {
            for (int j = js; j <= je; j++)
            {
                for (int i = is; i <= ie; i++)
                {
                    int idx = m_idx_proxy->GetIdx(i, j, k);
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
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        out << "variables = x, y, z, type" << std::endl;
        out << "Zone T=test\n I=" << ie - is + 1 << " J=" << je - js + 1 << " K=" << ke - ks + 1 << " F=POINT" << std::endl;
        for (int k = ks; k <= ke; k++)
        {
            for (int j = js; j <= je; j++)
            {
                for (int i = is; i <= ie; i++)
                {
                    int idx = m_idx_proxy->GetIdx(i, j, k);
                    out << node->GetCoord(i, j, k)[0] << " " << node->GetCoord(i, j, k)[1] << " " << node->GetCoord(i, j, k)[2] << " " << (int)m_node_type[idx] << std::endl;
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
        out << " N=" << m_fn_info.node[m_layer_num].size() << " E=" << m_trans_face.size() << " ZONETYPE=FEQuadrilateral" << std::endl;
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
        out << " N=" << m_fn_info.node[1].size() << " E=" << m_trans_face.size() << " ZONETYPE=FEQuadrilateral" << std::endl;
        out << "DATAPACKING=POINT" << std::endl;
        for (auto &coord : m_fn_info.node[1])
        {
            out << coord.coord[0] << " " << coord.coord[1] << " " << coord.coord[2] << std::endl;
        }
        for (int iFace = 0; iFace < m_trans_face.size(); iFace++)
        {
            // auto face_node_idx = m_trans_face[iFace].idx_master;
            // std::vector<int> idx_layer;
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
        out << "variables = x, y, z" << std::endl;
        out << "Zone T=test" << std::endl;
        int node_num = 0;
        for (int iLayer = 0; iLayer < m_fn_info.node.size(); iLayer++)
        {
            node_num += m_fn_info.node[iLayer].size();
        }
        int cell_num = 0;
        cell_num = m_fn_info.cell.size();
        out << " N=" << node_num << " E=" << cell_num << " F=FEPOINT, ET=Brick" << std::endl;
        for (int iLayer = 0; iLayer < m_fn_info.node.size(); iLayer++)
        {
            for (int iNode = 0; iNode < m_fn_info.node[iLayer].size(); iNode++)
            {
                out << m_fn_info.node[iLayer][iNode].coord[0] << " " << m_fn_info.node[iLayer][iNode].coord[1] << " " << m_fn_info.node[iLayer][iNode].coord[2] << std::endl;
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

    void GridFNFactoryZaran::BuildFNGridInfo()
    {
        BuildFNNodeCoord();
        WriteProjectNode();
        WriteModelSurface();
        WriteTransFace();
        BuildCell();
        WriteSlaveGrid();
        BuildNodeNeighbor();
    }

    void GridFNFactoryZaran::BuildFNNodeCoord()
    {

        m_fn_info.node.resize(m_layer_num + 1);
        BuildRefNode();
        BuildTransNode();
        BuildWallNode();
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
        NodeFN *node = new NodeFN(node_num);
        for (int iLayer = 0; iLayer < m_fn_info.node.size(); iLayer++)
        {
            for (int iNode = 0; iNode < m_fn_info.node[iLayer].size(); iNode++)
            {
                int idx = m_fn_info.node[iLayer][iNode].idx;
                auto coord = m_fn_info.node[iLayer][iNode].coord;
                node->SetCoord(idx, coord);
                if (iLayer == 0)
                {
                    node->SetType(idx, NodeType::ghost);
                }
                else if (iLayer == m_layer_num)
                {
                    node->SetType(idx, NodeType::wall);
                }
                else
                {
                    node->SetType(idx, NodeType::inner);
                }
            }
        }
        grid->SetNode(node);
    }

    void GridFNFactoryZaran::SetFNGridNodeNeighbor()
    {
        auto grid = GetFNGrid();
        auto node = grid->GetNode();
        int node_num = node->GetCount();
        std::vector<int> neighbor_node_num;
        neighbor_node_num.resize(node_num);
        for (int iLayer = 0; iLayer < m_fn_info.node.size(); iLayer++)
        {
            for (int iNode = 0; iNode < m_fn_info.node[iLayer].size(); iNode++)
            {
                int idx = m_fn_info.node[iLayer][iNode].idx;
                neighbor_node_num[idx] = m_fn_info.node[iLayer][iNode].neighbor_node.size();
            }
        }
        int total_neighbor_num = 0;
        for (int i = 0; i < node_num; i++)
        {
            total_neighbor_num += neighbor_node_num[i];
        }
        std::vector<int> neighbor_node_idx;
        neighbor_node_idx.reserve(total_neighbor_num);
        for (int iLayer = 0; iLayer < m_fn_info.node.size(); iLayer++)
        {
            for (int iNode = 0; iNode < m_fn_info.node[iLayer].size(); iNode++)
            {
                neighbor_node_idx.insert(neighbor_node_idx.end(), m_fn_info.node[iLayer][iNode].neighbor_node.begin(), m_fn_info.node[iLayer][iNode].neighbor_node.end());
            }
        }
        node->SetNeighborNode(node_num, neighbor_node_num.data(), neighbor_node_idx.data());
    }

    void GridFNFactoryZaran::SetFNGridCell()
    {
        auto grid = GetFNGrid();
        int cell_num = m_fn_info.cell.size();
        CellFN *cell = new CellFN(cell_num);
        cell->SetNode(m_fn_info.cell);
        grid->SetCell(cell);
    }

    void GridFNFactoryZaran::SetFNGridBoundary()
    {
        auto grid = GetFNGrid();
        BoundMapFN *bound = new BoundMapFN();
        bound->CreateBoundary("wall");
        auto &wall = bound->GetBoundary("wall");
        int node_num = m_fn_info.node[m_layer_num].size();
        wall.resize(node_num);
        for (int iNode = 0; iNode < node_num; iNode++)
        {
            int bound_idx = m_fn_info.node[m_layer_num][iNode].idx;
            int ref_idx = m_fn_info.node[1][iNode].idx;
            double normal[3];
            auto bound_coord = m_fn_info.node[m_layer_num][iNode].coord;
            auto ref_coord = m_fn_info.node[1][iNode].coord;
            normal[0] = bound_coord[0] - ref_coord[0];
            normal[1] = bound_coord[1] - ref_coord[1];
            normal[2] = bound_coord[2] - ref_coord[2];
            double len = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
            normal[0] /= len;
            normal[1] /= len;
            normal[2] /= len;
            wall[iNode].SetIdxBound(bound_idx);
            wall[iNode].SetIdxRef(ref_idx);
            wall[iNode].SetNorm(normal);
        }
        grid->SetBoundaryMap(bound);
    }

    void GridFNFactoryZaran::SetFNGridBoundaryFace()
    {
        auto grid = GetFNGrid();
        FaceFN *face = new FaceFN();
        int face_num = m_trans_face.size();
        std::vector<int> face_node_num;
        face_node_num.resize(face_num);
        for (int iFace = 0; iFace < face_num; iFace++)
        {
            face_node_num[iFace] = m_trans_face[iFace].idx_master.size();
        }
        face->Allocate(face_num, face_node_num.data());
        double normal[3];
        double area;
        for (int iFace = 0; iFace < face_num; iFace++)
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

            for (int iNode = 0; iNode < face_node_idx.size(); iNode++)
            {
                face_node_idx[iNode] = m_fn_info.node[m_layer_num][face_node_idx[iNode]].idx;
            }
            face->SetFace2Node(iFace, face_node_idx.data(), face_node_idx.size());
            face->SetNormal(iFace, normal);
            face->SetArea(iFace, area);
        }

        grid->SetFace(face);
    }

    void GridFNFactoryZaran::BuildRefNode()
    {
        auto grid = GetBlockGrid();
        auto node = grid->GetNode();
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        std::set<int> ref_node_idx_set;
        for (int k = ks; k <= ke; k++)
        {
            for (int j = js; j <= je; j++)
            {
                for (int i = is; i <= ie; i++)
                {
                    int idx = m_idx_proxy->GetIdx(i, j, k);
                    if (m_node_type[idx] == PhysicalType::FluidSolid)
                    {
                        if (m_node_type[m_idx_proxy->GetIdx(i + 1, j, k)] == PhysicalType::Fluid)
                        {
                            ref_node_idx_set.insert(m_idx_proxy->GetIdx(i + 1, j, k));
                        }
                        if (m_node_type[m_idx_proxy->GetIdx(i - 1, j, k)] == PhysicalType::Fluid)
                        {
                            ref_node_idx_set.insert(m_idx_proxy->GetIdx(i - 1, j, k));
                        }
                        if (m_node_type[m_idx_proxy->GetIdx(i, j + 1, k)] == PhysicalType::Fluid)
                        {
                            ref_node_idx_set.insert(m_idx_proxy->GetIdx(i, j + 1, k));
                        }
                        if (m_node_type[m_idx_proxy->GetIdx(i, j - 1, k)] == PhysicalType::Fluid)
                        {
                            ref_node_idx_set.insert(m_idx_proxy->GetIdx(i, j - 1, k));
                        }
                        if (m_node_type[m_idx_proxy->GetIdx(i, j, k + 1)] == PhysicalType::Fluid)
                        {
                            ref_node_idx_set.insert(m_idx_proxy->GetIdx(i, j, k + 1));
                        }
                        if (m_node_type[m_idx_proxy->GetIdx(i, j, k - 1)] == PhysicalType::Fluid)
                        {
                            ref_node_idx_set.insert(m_idx_proxy->GetIdx(i, j, k - 1));
                        }
                    }
                }
            }
        }
        m_fn_info.node[0].resize(ref_node_idx_set.size());
        int idx = 0;
        for (auto &ref_node_idx : ref_node_idx_set)
        {
            int i, j, k;
            m_idx_proxy->GetIdxStruct(ref_node_idx, i, j, k);
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
        // int is, ie, js, je, ks, ke;
        // grid->GetRange(is, ie, js, je, ks, ke);
        std::set<int> trans_node_idx_set;
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
            auto face_node_idx = m_trans_face[iFace].idx_master;
            for (int iNode = 0; iNode < face_node_idx.size(); iNode++)
            {
                trans_node_idx_set.insert(face_node_idx[iNode]);
            }
        }

        m_fn_info.node[1].resize(trans_node_idx_set.size());
        int idx = 0;
        for (auto &trans_idx : trans_node_idx_set)
        {
            int i, j, k;
            m_idx_proxy->GetIdxStruct(trans_idx, i, j, k);
            m_fn_info.node[1][idx].coord[0] = node->GetCoord(i, j, k)[0];
            m_fn_info.node[1][idx].coord[1] = node->GetCoord(i, j, k)[1];
            m_fn_info.node[1][idx].coord[2] = node->GetCoord(i, j, k)[2];
            m_fn_info.node[1][idx].idx = idx + m_fn_info.node[0].size();
            m_trans_node.insert(ConnectInfo{trans_idx, 1, idx});
            idx++;
        }
        for (int iFace = 0; iFace < m_trans_face.size(); iFace++)
        {
            auto idx_master = m_trans_face[iFace].idx_master;
            auto &idx_slave = m_trans_face[iFace].idx_slave;
            idx_slave.resize(idx_master.size());
            for (int iNode = 0; iNode < idx_slave.size(); iNode++)
            {
                for (auto &slave : m_trans_node)
                {
                    if (slave.idx_block == idx_master[iNode])
                    {
                        idx_slave[iNode] = slave.idx_local_layer;
                        break;
                    }
                }
            }
        }
    }
    void GridFNFactoryZaran::CheckTransFace()
    {
        for (int iFace = 0; iFace < m_trans_face.size(); iFace++)
        {
            auto face_node_idx = m_trans_face[iFace].idx_master;
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
            double proj_coord[3];
            m_model_manager->GetClosestPoint(trans_coord, proj_coord);
            m_fn_info.node[m_layer_num][iNode].coord[0] = proj_coord[0];
            m_fn_info.node[m_layer_num][iNode].coord[1] = proj_coord[1];
            m_fn_info.node[m_layer_num][iNode].coord[2] = proj_coord[2];
            m_fn_info.node[m_layer_num][iNode].idx = iNode + m_fn_info.node[0].size() + m_fn_info.node[1].size() * (m_layer_num - 1);
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
                m_fn_info.node[iLayer][iNode].idx = iNode + m_fn_info.node[0].size() + m_fn_info.node[1].size() * (iLayer - 1);
            }
        }
    }
    void GridFNFactoryZaran::BuildNodeNeighbor()
    {
        BuildProjectNodeNeighbor();
        BuildTransNodeNeighbor();
        ReorderProjectNodeNeighbor();
        CheckProjectNodeNeighbor();
    }
    void GridFNFactoryZaran::BuildProjectNodeNeighbor()
    {
        // direct neighbor
        std::vector<std::set<int>> node_neighbor_origin;
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
        // extend neighbor, if the neighbor node is less than 3, extend the neighbor node
        std::vector<std::set<int>> node_neighbor_extend;
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
            if (iNode == 508)
            {
                int a = 0;
            }
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
                //         vec1[i] = node->GetCoord(neighbor[0])[i] - node->GetCoord(neighbor[1])[i];
                //         vec2[i] = node->GetCoord(neighbor[2])[i] - node->GetCoord(neighbor[3])[i];
                //         vec3[i] = node->GetCoord(neighbor[4])[i] - node->GetCoord(neighbor[5])[i];
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
                //         vec1[i] = node->GetCoord(neighbor[1])[i] - node->GetCoord(neighbor[0])[i];
                //         vec2[i] = node->GetCoord(neighbor[3])[i] - node->GetCoord(neighbor[2])[i];
                //         vec3[i] = node->GetCoord(neighbor[5])[i] - node->GetCoord(neighbor[4])[i];
                //     }
                //     volume = vec1.dot(vec2.cross(vec3));
                //     Log::info("volume:{}", volume);
                // }

                node_pair main_pair;
                main_pair.node1 = neighbor[0];
                main_pair.node2 = neighbor[1];
                DVector3D main_vec;
                for (int i = 0; i < 3; ++i)
                {
                    main_vec[i] = node->GetCoord(main_pair.node2)[i] - node->GetCoord(main_pair.node1)[i];
                }
                main_vec.normalize();
                neighbor.erase(std::find(neighbor.begin(), neighbor.end(), main_pair.node1));
                neighbor.erase(std::find(neighbor.begin(), neighbor.end(), main_pair.node2));
                // 求出所有邻居节点在以主方向向量为法向量，经过当地节点的平面上的投影
                map<int, DVector3D> node_proj_map;
                DVector3D vec;
                for (int i = 0; i < neighbor.size(); ++i)
                {
                    for (int k = 0; k < 3; ++k)
                    {
                        vec[k] = node->GetCoord(neighbor[i])[k] - node->GetCoord(idx)[k];
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
                // 以第一个邻居节点投影向量为基准向量，求出基准向量以法向量为旋转轴旋转到其他投影向量的角�?, 0~2pi
                map<double, int> node_angle_map;
                for (int i = 0; i < neighbor.size(); ++i)
                {
                    if (i == 0)
                    {
                        node_angle_map[0] = neighbor[i];
                        continue;
                    }
                    DVector3D vec = node_proj_map[neighbor[i]];
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
                auto get_next_node = [&](int iNode, IArray neiborNode) -> int
                {
                    if (iNode == neiborNode.size() - 1)
                        return 0;
                    else
                        return iNode + 1;
                };
                // 获取上一个点的lamda表达�?
                auto get_last_node = [&](int iNode, IArray neiborNode) -> int
                {
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
                    DVector3D vec1, vec2;
                    vec1 = node_proj_map[temp.node1];
                    vec2 = node_proj_map[temp.node2];
                    double angle = AngleOfTwoArray3D(vec1.data(), vec2.data()) + GetRand(0.0, 1.0) * EPSILON_NUMBER;
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
                    int remove_node, remove_index;
                    node_pair temp;
                    if (node_dis_map[temp_pair.node1] > node_dis_map[temp_pair.node2])
                    {
                        remove_node = temp_pair.node1;
                        remove_index = std::find(neighbor.begin(), neighbor.end(), remove_node) - neighbor.begin();
                        for (auto i : node_pair_map)
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
                        for (auto i : node_pair_map)
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
                    DVector3D vec1, vec2;
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
            }
        }
    }
    void GridFNFactoryZaran::CheckProjectNodeNeighbor()
    {
        auto grid = GetFNGrid();
        auto node = grid->GetNode();
        int node_num = m_fn_info.node[1].size();
        double delta = 15*PI / 180;
        for (int iLayer = 1; iLayer < m_layer_num; iLayer++)
        {
            for (int iNode = 0; iNode < node_num; iNode++)
            {
                int idx = m_fn_info.node[iLayer][iNode].idx;
                if (iNode == 508)
                {
                    int a = 0;
                }
                auto &neighbor = m_fn_info.node[iLayer][iNode].neighbor_node;
                Array<DVector3D> vec(3);
                for (int i = 0; i < 3; i++)
                {
                    vec[0][i] = node->GetCoord(neighbor[1])[i] - node->GetCoord(neighbor[0])[i];
                    vec[1][i] = node->GetCoord(neighbor[3])[i] - node->GetCoord(neighbor[2])[i];
                    vec[2][i] = node->GetCoord(neighbor[5])[i] - node->GetCoord(neighbor[4])[i];
                }
                double volume = vec[0].cross(vec[1]).dot(vec[2]);
                // 检查是否是右手坐标系
                if (vec[0].cross(vec[1]).dot(vec[2]) < 0)
                {
                    std::swap(neighbor[0], neighbor[1]);
                }
                double angle = AngleOfTwoArray3D(vec[0].data(), vec[1].data());
                // i,j方向平行
                if (abs(angle) < delta)
                {
                    std::swap(neighbor[1], neighbor[3]);
                }
                else if (abs(angle - PI) < delta)
                {
                    std::swap(neighbor[1], neighbor[2]);
                }

                angle = AngleOfTwoArray3D(vec[0].data(), vec[2].data());
                // i,k方向平行
                if (abs(angle) < delta)
                {
                    std::swap(neighbor[1], neighbor[5]);
                }
                else if (abs(angle - PI) < delta)
                {
                    std::swap(neighbor[1], neighbor[4]);
                }
                angle = AngleOfTwoArray3D(vec[1].data(), vec[2].data());
                double angle1 =angle*180/PI;
                // j,k方向平行
                if (abs(angle) < delta)
                {
                    std::swap(neighbor[3], neighbor[5]);
                }
                else if (abs(angle - PI) < delta)
                {
                    std::swap(neighbor[3], neighbor[4]);
                }
                // 检查是否是右手坐标系
                for (int i = 0; i < 3; i++)
                {
                    vec[0][i] = node->GetCoord(neighbor[1])[i] - node->GetCoord(neighbor[0])[i];
                    vec[1][i] = node->GetCoord(neighbor[3])[i] - node->GetCoord(neighbor[2])[i];
                    vec[2][i] = node->GetCoord(neighbor[5])[i] - node->GetCoord(neighbor[4])[i];
                }
                if (vec[0].cross(vec[1]).dot(vec[2]) < 0)
                {
                    std::swap(neighbor[4], neighbor[5]);
                }
                // 检查是否是右手坐标系
                for (int i = 0; i < 3; i++)
                {
                    vec[0][i] = node->GetCoord(neighbor[1])[i] - node->GetCoord(neighbor[0])[i];
                    vec[1][i] = node->GetCoord(neighbor[3])[i] - node->GetCoord(neighbor[2])[i];
                    vec[2][i] = node->GetCoord(neighbor[5])[i] - node->GetCoord(neighbor[4])[i];
                }
                if (idx == 11400)
                {
                    Log::error("node:{}, iLayer;{}, iNode:{}, coord:{}, {}, {}, volume:{}", idx, iLayer, iNode, node->GetCoord(idx)[0], node->GetCoord(idx)[1], node->GetCoord(idx)[2], vec[0].cross(vec[1]).dot(vec[2]));
                }
            }
        }
    }
    void GridFNFactoryZaran::BuildTransNodeNeighbor()
    {
        for (int iNode = 0; iNode < m_fn_info.node[1].size(); iNode++)
        {
            m_fn_info.node[1][iNode].neighbor_node.resize(6);
        }
        for (auto &nodes : m_trans_node)
        {
            int i, j, k;
            m_idx_proxy->GetIdxStruct(nodes.idx_block, i, j, k);
            m_fn_info.node[1][nodes.idx_local_layer].neighbor_node[0] = m_idx_proxy->GetIdx(i - 1, j, k);
            m_fn_info.node[1][nodes.idx_local_layer].neighbor_node[1] = m_idx_proxy->GetIdx(i + 1, j, k);
            m_fn_info.node[1][nodes.idx_local_layer].neighbor_node[2] = m_idx_proxy->GetIdx(i, j - 1, k);
            m_fn_info.node[1][nodes.idx_local_layer].neighbor_node[3] = m_idx_proxy->GetIdx(i, j + 1, k);
            m_fn_info.node[1][nodes.idx_local_layer].neighbor_node[4] = m_idx_proxy->GetIdx(i, j, k - 1);
            m_fn_info.node[1][nodes.idx_local_layer].neighbor_node[5] = m_idx_proxy->GetIdx(i, j, k + 1);
        }
        for (int iNode = 0; iNode < m_fn_info.node[1].size(); iNode++)
        {
            if (m_fn_info.node[1][iNode].idx == 1914)
            {
                int a = 0;
            }
            int next_layer_node = m_fn_info.node[2][iNode].idx;
            for (int iNeighbor = 0; iNeighbor < 6; iNeighbor++)
            {
                int idx_master = m_fn_info.node[1][iNode].neighbor_node[iNeighbor];
                if (m_node_type[idx_master] == PhysicalType::FluidSolid)
                {
                    for (auto &nodes : m_trans_node)
                    {
                        if (nodes.idx_block == idx_master)
                        {
                            m_fn_info.node[1][iNode].neighbor_node[iNeighbor] = m_fn_info.node[1][nodes.idx_local_layer].idx;
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
                            m_fn_info.node[1][iNode].neighbor_node[iNeighbor] = m_fn_info.node[0][nodes.idx_local_layer].idx;
                            break;
                        }
                    }
                }
                else if (m_node_type[idx_master] == PhysicalType::Solid)
                {
                    m_fn_info.node[1][iNode].neighbor_node[iNeighbor] = next_layer_node;
                }
                else
                {
                    Log::error("Invalid node type: {}", idx_master);
                }
            }
        }
    }
    void GridFNFactoryZaran::BuildCell()
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
}