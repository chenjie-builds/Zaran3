#pragma once
#include "GridGeneratorBlock.h"
namespace zaran {
	class GridGeneratorBlockExplosion : public GridGeneratorBlock
	{
	public:
		GridGeneratorBlockExplosion();
		~GridGeneratorBlockExplosion();
		void CreateGrid(shared_ptr<GridBlock> grid, GridBlockInfo& grid_info);
	private:
		//根据建筑物外形设置iBlank
		void SetIBlankByBuilding();
		//输出测试
		void ExportTecplot();
	};
} // namespace zaran