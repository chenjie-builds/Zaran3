/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file GridGeneratorFlexibleZaran.h
 * \brief GridGeneratorFlexibleZaran class, used to generate flexible grid in Zaran.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "FaceFNFDM.h"
#include "GridFNFDM.h"
#include "NodeFNFDM.h"
#include "GridBlock.h"
#include "GridGenerator.h"
#include "ModelManager.h"
#include "StructIdxProxy.h"
#include <set>
namespace zaran
{
    enum class PhysicalType
    {
        Fluid = 1,      // 流体单元
        Solid = 2,      // 固体单元
        FluidSolid = 3, // 流体单元，邻居单元有固体单元
        SolidFluid = 4, // 固体单元，邻居单元有流体单元
        Unset = 5,      // 未知单元
    };
    //过渡面元
    struct TransFace
    {
        dynamic_array<index_type> idx_block;
        bool operator<(const TransFace &rhs) const
        {
            return idx_block < rhs.idx_block;
        }
        dynamic_array<index_type> idx_slave;
    };
    struct SlaveNode
    {
        // 节点索引
        index_type idx;
        // 坐标
        double coord[3];
        // 邻居节点记录其投影层数索引
        dynamic_array<index_type> neighbor_node;
    };
    struct ConnectInfo
    {
        // 在Block中的索引
        index_type idx_block;
        // 在slave grid中的索引,第几层
        index_type idx_n_layers;
        // 在slave grid中的索引,该层的第几个节点
        index_type idx_local_layer;
        bool operator<(const ConnectInfo &rhs) const
        {
            return idx_block < rhs.idx_block;
        }
    };
    struct FNGridInfo
    {
        dynamic_array<dynamic_array<SlaveNode>> node;
        dynamic_array<dynamic_array<index_type>> cell;
    };

    class GridFNFactoryZaran : public GridGenerator
    {
    public:
        GridFNFactoryZaran() = default;
        void CreateGrid(const shared_ptr<GridBlock> &block, const shared_ptr<GridFN> &grid, const shared_ptr<ModelManager> &model_manager);
        ~GridFNFactoryZaran() override = default;

    public:
        shared_ptr<GridBlock> GetBlockGrid() { return m_block_grid; }
        shared_ptr<GridFN> GetFNGrid() { return m_fn_grid; }
        shared_ptr<ModelManager> GetModelManager() { return m_model_manager; }
        shared_ptr<IdProxyStruct> GetIdxProxy() { return m_idx_proxy; }
        dynamic_array<PhysicalType> &GetCellType() { return m_cell_type; }
        dynamic_array<PhysicalType> &GetNodeType() { return m_node_type; }
        dynamic_array<TransFace> &GetTransFace() { return m_trans_face; }
        std::set<ConnectInfo> &GetRefNode() { return m_ref_node; }
        std::set<ConnectInfo> &GetTransNode() { return m_trans_node; }
        FNGridInfo &GetFNGridInfo() { return m_fn_info; }

    private:
        void TagBlockGrid();
        void TagCells();
        void ProcessCell(index_type start_i, index_type end_i, index_type start_j, index_type end_j, index_type start_k, index_type end_k);
		//标记流体单元(i,j,k)对应的邻居节点
		void TagFluidCells(int i, int j, int k);
        void TagNodes3D();
		void TagNodes();
        void TagNodes2D();
        void ReTagBlockGrid();
        void ReTagCells();
        void SetNodeTag();
        void WriteNodeTag();
        void WriteProjectNode();
        void WriteModelSurface();
        void WriteTransFace();
        void WriteSlaveGrid();
        void WriteCellTag();

        void BuildFNGridInfo();
        void BuildFNNodeCoord();
        void SetFNGrid();
        void SetFNGridNode();
        void SetFNGridNodeNeighbor();
        void SetFNGridCell();
        void SetFNGridBoundary();
        void SetFNGridBoundaryFace();

        // 生成参考节点,即过渡节点计算需要的内部点,在SlaveGrid第0层
        void BuildRefNode();
        // 生成过渡节点,即过渡节点计算需要的内部点,在SlaveGrid第1层
        void BuildTransNode();
        // 检查过渡节点,不能同时i+1和i-1都是固体节点
        bool CheckTransNode();
        bool CheckTransNode2D();
		bool CheckTransNode3D();
        void CheckTransFace();
        // 生成物面节点,即过渡节点计算需要的内部点,在SlaveGrid最后一层
        void BuildWallNode();
        // 生成中间节点,即过渡节点到物面节点之间的投影点
        void BuildProjectNode();
        // 优化物面节点位置，对于凹坑等，对于物面节点的坐标进行优化
        void OptimizeWallNode();
        void OptimizeWallNode2D();
        void OptimizeWallNode3D();
        void BuildNodeNeighbor();
        // 生成节点邻居信息,投影节点
        void BuildProjectNodeNeighbor();
        void BuildProjectNodeNeighbor2D();
        void BuildProjectNodeNeighbor3D();
        void ReorderProjectNodeNeighbor();

        void CheckProjectNodeNeighbor();
        void CheckProjectNodeNeighbor2D();
        void CheckProjectNodeNeighbor3D();
        // 生成节点邻居信息,过渡节点
        void BuildTransNodeNeighbor();
        void BuildTransNodeNeighbor2D();
        void BuildTransNodeNeighbor3D();
        void BuildFNCell();
        void BuildFNCell2D();
        void BuildFNCell3D();

    private:
        shared_ptr<GridBlock> m_block_grid;
        shared_ptr<GridFN> m_fn_grid;
        shared_ptr<ModelManager> m_model_manager;
        shared_ptr<IdProxyStruct> m_idx_proxy;
        dynamic_array<PhysicalType> m_cell_type;
        dynamic_array<PhysicalType> m_node_type;
        dynamic_array<TransFace> m_trans_face;
        // 保存参考节点, 与过渡节点连接的流体节点
        std::set<ConnectInfo> m_ref_node;
        // 保存过渡节点, 与物面节点连接的流体节点
        std::set<ConnectInfo> m_trans_node;
        FNGridInfo m_fn_info;

        int m_layer_num{};
    };
}