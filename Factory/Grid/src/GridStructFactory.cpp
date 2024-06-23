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
        m_bnd_file_name = "1.fvbnd";

    }

    GridStruct* GridStructFactory::CreateGrid()
    {
        m_ghost_size = 1;
        ReadNodeFile();
        ReadBoundFile();
        GridStruct* grid = new GridStruct("Structured", 1, m_dim);
        ConvertToGrid(grid);
        return grid;
    }

    void GridStructFactory::ConvertToGrid(GridStruct*& grid)
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

        // print node coord for debug
        for (int idx_k = 0;idx_k < m_nk + 2 * m_ghost_size;++idx_k)
        {
            for (int idx_j = 0;idx_j < m_nj + 2 * m_ghost_size;++idx_j)
            {
                for (int idx_i = 0;idx_i < m_ni + 2 * m_ghost_size;++idx_i)
                {
                    const double* coord = node->GetCoord(idx_i, idx_j, idx_k);
                    Log::info("Node Coord:({},{},{})=({},{},{})", idx_i, idx_j, idx_k, coord[0], coord[1], coord[2]);
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
        //TODO Read Bound File Here

    }
} // namespace zaran