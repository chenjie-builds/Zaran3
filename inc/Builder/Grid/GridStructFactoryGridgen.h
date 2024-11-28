#pragma once
#include "GridFactory.h"
#include "GridStruct.h"
#include "NodeStruct.h"
#include "FaceStruct.h"
#include "CellStruct.h"
#include "GlobalData.h"
#include <fstream>
namespace zaran
{
    /// @brief Generate structured grid using gridgen generic
    /// Multi-block structured grid
    class GridStructFactoryGridgen : public GridBuilder
    {
    public:
        GridStructFactoryGridgen();
        void CreateGrid(GridBase **&grid_list, int &grid_num) override;
        void AllocateGridMemory(GridBase **&grid_list, int &grid_num);

    private:
        void ReadNodeFile();
        void ReadBoundFile();
        void SetMultiBlockInfo();
        void SetBoundInfo(GridBase **grid_list, int &grid_num);
        void SetNodeCoord(GridBase **grid_list, int &grid_num);
        void SetGhostNodeCoord(GridBase **grid_list, int &grid_num);
        void SetGhostNodeCoord3D(GridBase **grid_list, int &grid_num);
        void SetGhostNodeCoord2D(GridBase **grid_list, int &grid_num);
        size_t GetBlockNum() { return m_block.size(); };
        void WriteGridTest(GridBase **grid_list, int &grid_num);
        struct NodeCoord
        {
            double coord[3];
        };
        struct BoundInfo
        {
            // source block info
            //  source block index
            int idx_block_source;
            // boundary type
            int bound_type;
            // start and end index of i, j, k direction
            int is_s, ie_s, js_s, je_s, ks_s, ke_s;
            // boundary direction -1 left boundary, 1 right boundary, 0 not boundary
            int dir_s[3];

            // target block info
            // target block index
            int block_indx_target;
            // start and end index of i, j, k direction
            int is_t, ie_t, js_t, je_t, ks_t, ke_t;
            // boundary direction -1 left boundary, 1 right boundary, 0 not boundary
            int dir_t[3];

            // block connection info, block_conn_info[i]: direction i in source block is connected to target block
            // 0: connect to i direction, 1: connect to j direction, 2: connect to k direction
            int conn_info[3];
        };
        struct BlockInfo
        {
            int ni, nj, nk;
            std::vector<std::vector<std::vector<NodeCoord>>> node_coord;
            std::vector<BoundInfo> bound_info;
        };

    private:
        int m_dim;
        int m_ghost_size;
        std::vector<BlockInfo> m_block;
        string m_node_file_name;
        string m_bnd_file_name;
    };
}