#include "GridStructFactory.h"
#include "Log.h"
#include <fstream>
namespace zaran
{
    GridStructFactory::GridStructFactory()
    {
        // m_node_file_name = GlobalData::GetString("NodeFile");
        // m_bnd_file_name = GlobalData::GetString("BoundFile");
        m_node_file_name = "1.dat";
        m_bnd_file_name = "1.inp";

    }

    GridStruct* GridStructFactory::CreateGrid()
    {
        m_ghost_size = 1;
        ReadNodeFile();
        ReadBoundFile();
        GridStruct* grid = new GridStruct("Structured", 1, m_dim);
        ConvertToGrid(grid);
        SetBoundInfo(grid);
        return grid;
    }

    void GridStructFactory::ConvertToGrid(GridStruct* grid)
    {
        grid->Allocate(m_ni, m_nj, m_nk, m_ghost_size);
        auto node = grid->GetNode();
        int is, ie, js, je, ks, ke;
        grid->GetRange(is, ie, js, je, ks, ke);
        for (int idx_k = ks;idx_k <= ke;++idx_k)
        {
            for (int idx_j = js;idx_j <= je;++idx_j)
            {
                for (int idx_i = is;idx_i <= ie;++idx_i)
                {
                    node->SetCoord(idx_i, idx_j, idx_k, m_node_coord[idx_i - m_ghost_size][idx_j - m_ghost_size][idx_k - m_ghost_size].coord);
                }
            }
        }
        //Set Ghost Node Coord
        const double* ref_coord;
        int ref_idx_i, ref_idx_j, ref_idx_k;
        int di, dj, dk;
        double coord[3];
        int ni, nj, nk;
        ni = grid->GetNi();
        nj = grid->GetNj();
        nk = grid->GetNk();
        for (int idx_k = 0;idx_k < nk;++idx_k)
        {
            if (idx_k < m_ghost_size)
            {
                ref_idx_k = 2 * m_ghost_size - idx_k;
                dk = idx_k - m_ghost_size;
            }
            else if (idx_k >= nk - m_ghost_size)
            {
                ref_idx_k = 2 * (nk - m_ghost_size - 1) - idx_k;
                dk = idx_k - (nk - 1 - m_ghost_size);
            }
            else
            {
                ref_idx_k = idx_k;
                dk = 0;
            }
            for (int idx_j = 0;idx_j < nj;++idx_j)
            {
                if (idx_j < m_ghost_size)
                {
                    ref_idx_j = 2 * m_ghost_size - idx_j;
                    dj = idx_j - m_ghost_size;
                }
                else if (idx_j >= nj - m_ghost_size)
                {
                    ref_idx_j = 2 * (nj - m_ghost_size - 1) - idx_j;
                    dj = idx_j - (nj - 1 - m_ghost_size);
                }
                else
                {
                    ref_idx_j = idx_j;
                    dj = 0;
                }
                for (int idx_i = 0;idx_i < ni;++idx_i)
                {
                    if (idx_i < m_ghost_size)
                    {
                        ref_idx_i = 2 * m_ghost_size - idx_i;
                        di = idx_i - m_ghost_size;
                    }
                    else if (idx_i >= ni - m_ghost_size)
                    {
                        ref_idx_i = 2 * (ni - m_ghost_size - 1) - idx_i;
                        di = idx_i - (ni - 1 - m_ghost_size);
                    }
                    else
                    {
                        ref_idx_i = idx_i;
                        di = 0;
                    }
                    ref_coord = node->GetCoord(ref_idx_i, ref_idx_j, ref_idx_k);
                    coord[0] = ref_coord[0] + 2 * (node->GetCoord(ref_idx_i + di, ref_idx_j, ref_idx_k)[0] - ref_coord[0]);
                    coord[1] = ref_coord[1] + 2 * (node->GetCoord(ref_idx_i, ref_idx_j + dj, ref_idx_k)[1] - ref_coord[1]);
                    coord[2] = ref_coord[2] + 2 * (node->GetCoord(ref_idx_i, ref_idx_j, ref_idx_k + dk)[2] - ref_coord[2]);
                    node->SetCoord(idx_i, idx_j, idx_k, coord);
                }
            }
        }
    }

    void GridStructFactory::ReadNodeFile()
    {
        std::ifstream node_file(m_node_file_name);
        if (!node_file.is_open())
        {
            Log::error("Can't open node file:{}, Please Check!", m_node_file_name);
            system("pause");
        }
        int block_num;
        node_file >> block_num;
        node_file >> m_ni >> m_nj >> m_nk;
        // 设置网格的维度
        if (m_nk == 1)
        {
            m_dim = 2;
        }
        else
        {
            m_dim = 3;
        }
        m_node_coord.resize(m_ni);
        for (int idx_i = 0; idx_i < m_ni; ++idx_i)
        {
            m_node_coord[idx_i].resize(m_nj);
            for (int idx_j = 0; idx_j < m_nj; ++idx_j)
            {
                m_node_coord[idx_i][idx_j].resize(m_nk);
            }
        }
        for (int idx_k = 0;idx_k < m_nk;++idx_k)
        {
            for (int idx_j = 0;idx_j < m_nj;++idx_j)
            {
                for (int idx_i = 0;idx_i < m_ni;++idx_i)
                {
                    node_file >> m_node_coord[idx_i][idx_j][idx_k].coord[0];
                }
            }
        }
        for (int idx_k = 0;idx_k < m_nk;++idx_k)
        {
            for (int idx_j = 0;idx_j < m_nj;++idx_j)
            {
                for (int idx_i = 0;idx_i < m_ni;++idx_i)
                {
                    node_file >> m_node_coord[idx_i][idx_j][idx_k].coord[1];
                }
            }
        }
        for (int idx_k = 0;idx_k < m_nk;++idx_k)
        {
            for (int idx_j = 0;idx_j < m_nj;++idx_j)
            {
                for (int idx_i = 0;idx_i < m_ni;++idx_i)
                {
                    node_file >> m_node_coord[idx_i][idx_j][idx_k].coord[2];
                }
            }
        }
        node_file.close();
    }
    void GridStructFactory::ReadBoundFile()
    {
        // read inp file
        std::ifstream bnd_file(m_bnd_file_name);
        if (!bnd_file.is_open())
        {
            Log::error("Can't open bound file:{}, Please Check!", m_bnd_file_name);
            system("pause");
        }
        int block_num;
        bnd_file >> block_num;
        for (int idx_block = 0;idx_block < block_num;++idx_block)
        {
            int block_idx;
            bnd_file >> block_idx;
            int ni, nj, nk;
            bnd_file >> ni >> nj >> nk;
            std::string block_name;
            bnd_file >> block_name;
            int bound_num;
            bnd_file >> bound_num;
            for (int idx_bound = 0;idx_bound < bound_num;++idx_bound)
            {
                BoundInfo bound_info;
                bound_info.block_indx = block_idx;
                bnd_file >> bound_info.i_start >> bound_info.i_end;
                bnd_file >> bound_info.j_start >> bound_info.j_end;
                bnd_file >> bound_info.k_start >> bound_info.k_end;
                bnd_file >> bound_info.bound_type;
                m_bound_info.push_back(bound_info);
            }
        }
    }
    void GridStructFactory::SetBoundInfo(GridStruct* grid)
    {
        std::map<int, string> gridgen_bound = { {0,"none"},{1,"interblock_connection"},{2,"wall"},{3,"symmetry"},{4,"farfield"},{5,"inlet"},{6,"outlet"} };
        auto bound_map = grid->GetBoundMap();
        auto node = grid->GetNode();
        int i_bound, j_bound, k_bound;
        int i_ghost, j_ghost, k_ghost;
        int i_ref, j_ref, k_ref;
        double bound_norm[3];
        for (auto& bound_info : m_bound_info)
        {
            for (int idx_k = bound_info.k_start;idx_k <= bound_info.k_end;++idx_k)
            {
                for (int idx_j = bound_info.j_start;idx_j <= bound_info.j_end;++idx_j)
                {
                    for (int idx_i = bound_info.i_start;idx_i <= bound_info.i_end;++idx_i)
                    {
                        i_ghost = i_ref = i_bound = idx_i + m_ghost_size - 1;
                        j_ghost = j_ref = j_bound = idx_j + m_ghost_size - 1;
                        k_ghost = k_ref = k_bound = idx_k + m_ghost_size - 1;
                        if (bound_info.i_start == bound_info.i_end)
                        {
                            if (idx_i == 1)
                            {
                                i_ghost = i_bound - 1;
                                i_ref = i_bound + 1;
                            }
                            else
                            {
                                i_ghost = i_bound + 1;
                                i_ref = i_bound - 1;
                            }
                        }
                        if (bound_info.j_start == bound_info.j_end)
                        {
                            if (idx_j == 1)
                            {
                                j_ghost = j_bound - 1;
                                j_ref = j_bound + 1;
                            }
                            else
                            {
                                j_ghost = j_bound + 1;
                                j_ref = j_bound - 1;
                            }
                        }
                        if (bound_info.k_start == bound_info.k_end)
                        {
                            if (idx_k == 1)
                            {
                                k_ghost = k_bound - 1;
                                k_ref = k_bound + 1;
                            }
                            else
                            {
                                k_ghost = k_bound + 1;
                                k_ref = k_bound - 1;
                            }
                        }
                        bound_norm[0] = node->GetCoord(i_ref, j_ref, k_ref)[0] - node->GetCoord(i_bound, j_bound, k_bound)[0];
                        bound_norm[1] = node->GetCoord(i_ref, j_ref, k_ref)[1] - node->GetCoord(i_bound, j_bound, k_bound)[1];
                        bound_norm[2] = node->GetCoord(i_ref, j_ref, k_ref)[2] - node->GetCoord(i_bound, j_bound, k_bound)[2];
                        double normal = sqrt(bound_norm[0] * bound_norm[0] + bound_norm[1] * bound_norm[1] + bound_norm[2] * bound_norm[2]);
                        for (int i = 0;i < 3;++i)
                        {
                            bound_norm[i] /= normal;
                        }
                        BoundStruct bound(i_bound, j_bound, k_bound, i_ref, j_ref, k_ref, i_ghost, j_ghost, k_ghost, bound_norm);
                        bound_map->AddBoundary(gridgen_bound[bound_info.bound_type], bound);
                    }
                }
            }
        }
    }
} // namespace zaran