#include "GridStructFactory.h"
#include "Log.h"
#include <fstream>
namespace zaran
{
    GridStructFactory::GridStructFactory()
    {
        m_node_file_name = GlobalData::GetString("NodeFile");
        m_bnd_file_name = GlobalData::GetString("BoundFile");
    }

    GridStruct* GridStructFactory::CreateGrid()
    {
        ReadNodeFile();
        ReadBoundFile();
        GridStruct* grid = new GridStruct("Structured", 1, m_dim);
        return grid;
    }

    void GridStructFactory::ConvertToGrid(GridStruct*& grid)
    {
        auto node = grid->GetNode();
        if(node != nullptr)
        {
            delete node;
            node = nullptr;
        }
        node = new NodeStruct(m_ni, m_nj, m_nk);
        for (int k = 0;k < m_nk;++k)
        {
            for (int j = 0;j < m_nj;++j)
            {
                for (int i = 0;i < m_ni;++i)
                {
                    node->SetCoord(i, j, k, m_node_coord[i][j][k].coord);
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
        int i_num, j_num, k_num;
        node_file >> i_num >> j_num >> k_num;
        // 设置网格的维度
        if(k_num == 1)
        {
            m_dim = 2;
        }
        else
        {
            m_dim = 3;
        }
        m_node_coord.resize(i_num);
        for (int i = 0; i < i_num; ++i)
        {
            m_node_coord[i].resize(j_num);
            for (int j = 0; j < j_num; ++j)
            {
                m_node_coord[i][j].resize(k_num);
            }
        }
        for (int k = 0;k < k_num;++k)
        {
            for (int j = 0;j < j_num;++j)
            {
                for (int i = 0;i < i_num;++i)
                {
                    node_file >> m_node_coord[i][j][k].coord[0];
                }
            }
        }
        for (int k = 0;k < k_num;++k)
        {
            for (int j = 0;j < j_num;++j)
            {
                for (int i = 0;i < i_num;++i)
                {
                    node_file >> m_node_coord[i][j][k].coord[1];
                }
            }
        }
        for (int k = 0;k < k_num;++k)
        {
            for (int j = 0;j < j_num;++j)
            {
                for (int i = 0;i < i_num;++i)
                {
                    node_file >> m_node_coord[i][j][k].coord[2];
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