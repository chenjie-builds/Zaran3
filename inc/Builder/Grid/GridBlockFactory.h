#pragma once
#include "CommonPara.h"
#include "GridFactory.h"
#include "GridBlock.h"
namespace zaran
{
    struct GridBlockInfo
    {
        Box bound_box;
        int ni, nj, nk;
        std::string bound_type_i_minus, bound_type_i_plus;
        std::string bound_type_j_minus, bound_type_j_plus;
        std::string bound_type_k_minus, bound_type_k_plus;
    };
    class GridBlockFactory : public GridBuilder
    {
    public:
        GridBlockFactory() {}
        void CreateGrid(GridBlock *&grid, GridBlockInfo &grid_info);

    private:
        void AllocateGrid();
        void CreateNode();
        void CreateBound();
        GridBlock *GetGrid() { return m_grid; }

    private:
        GridBlockInfo m_grid_info;
        GridBlock *m_grid;
    };
}