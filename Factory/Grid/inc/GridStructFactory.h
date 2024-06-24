#pragma once
#include"GridFactory.h"
#include"GridStruct.h"
#include"NodeStruct.h"
#include"FaceStruct.h"
#include"CellStruct.h"
#include"GlobalData.h"
#include<vector>
namespace zaran
{
  
    class GridStructFactory :public GridFactory
    {
    public:
        GridStructFactory();
        GridStruct* CreateGrid();
        void ConvertToGrid(GridStruct* grid);
    private:
        void ReadNodeFile();
        void ReadBoundFile();
        void SetBoundInfo(GridStruct* grid);
        struct NodeCoord
        {
            double coord[3];
        };
        struct BoundInfo
        {
            int block_indx;
            int bound_type;
            int i_start, i_end, j_start, j_end, k_start, k_end;
        };
    private:
        int m_dim;
        int m_ni, m_nj, m_nk;
        int m_ghost_size;
        std::vector<std::vector<std::vector<NodeCoord>>> m_node_coord;
        std::vector<BoundInfo> m_bound_info;
        string m_node_file_name;
        string m_bnd_file_name;
    };
}