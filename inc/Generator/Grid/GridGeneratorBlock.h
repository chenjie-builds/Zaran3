#pragma once
#include "CommonPara.h"
#include "GridGenerator.h"
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
	// 网格生成器：块结构网格，均匀正交结构网格
    class GridGeneratorBlock : public GridGenerator
    {
    public:
        GridGeneratorBlock() {}
       void CreateGrid(shared_ptr<GridBlock> grid, GridBlockInfo &grid_info);

	protected:
        void AllocateGrid();
        void CreateNode();
        void CreateBound();
        shared_ptr<GridBlock> GetGrid() { return m_grid; }

    private:
        GridBlockInfo m_grid_info;
        shared_ptr<GridBlock> m_grid;
    };
}