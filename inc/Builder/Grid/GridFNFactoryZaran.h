#pragma once
#include "FaceFNFDM.h"
#include "GridFNFDM.h"
#include "NodeFNFDM.h"
#include "GridBlock.h"
#include "GridFactory.h"
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
    struct TransFace
    {
        std::vector<int> idx_block;
        bool operator<(const TransFace &rhs) const
        {
            return idx_block < rhs.idx_block;
        }
        std::vector<int> idx_slave;
    };
    struct SlaveNode
    {
        // 节点索引
        int idx;
        // 坐标
        double coord[3];
        // 邻居节点记录其投影层数索引
        std::vector<int> neighbor_node;
    };
    struct ConnectInfo
    {
        // 在Block中的索引
        int idx_block;
        // 在FN中的索引,第几层
        int idx_n_layers;
        // 在FN中的索引,该层的第几个节点
        int idx_local_layer;
        bool operator<(const ConnectInfo &rhs) const
        {
            return idx_block < rhs.idx_block;
        }
    };
    struct FNGridInfo
    {
        std::vector<std::vector<SlaveNode>> node;
        std::vector<std::vector<int>> cell;
    };

    class GridFNFactoryZaran : public GridBuilder
    {
    public:
        GridFNFactoryZaran() {}
        void CreateGrid(GridBlock *block, GridFN *&grid, ModelManager *model_manager);
        ~GridFNFactoryZaran() {}

    public:
        GridBlock *GetBlockGrid() { return m_block_grid; }
        GridFN *GetFNGrid() { return m_fn_grid; }
        ModelManager *GetModelManager() { return m_model_manager; }
        IdxStruct *GetIdxProxy() { return m_idx_proxy; }
        std::vector<PhysicalType> &GetCellType() { return m_cell_type; }
        std::vector<PhysicalType> &GetNodeType() { return m_node_type; }
        std::vector<TransFace> &GetTransFace() { return m_trans_face; }
        std::set<ConnectInfo> &GetRefNode() { return m_ref_node; }
        std::set<ConnectInfo> &GetTransNode() { return m_trans_node; }
        FNGridInfo &GetFNGridInfo() { return m_fn_info; }

    private:
        void TagBlockGrid();
        void TagCells();
        void ProcessCell(int start_i, int end_i, int start_j, int end_j, int start_k, int end_k);
        void TagNodes();
        void ReTagBlockGrid();
        void ReTagCells();
        void SetNodeTag();
        void WriteNodeTag();
        void WriteProjectNode();
        void WriteModelSurface();
        void WriteTransFace();
        void WriteSlaveGrid();

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
        void CheckTransFace();
        // 生成物面节点,即过渡节点计算需要的内部点,在SlaveGrid最后一层
        void BuildWallNode();
        // 生成中间节点,即过渡节点到物面节点之间的投影点
        void BuildProjectNode();

        void BuildNodeNeighbor();
        // 生成节点邻居信息,投影节点
        void BuildProjectNodeNeighbor();
        void ReorderProjectNodeNeighbor();
        void CheckProjectNodeNeighbor();
        // 生成节点邻居信息,过渡节点
        void BuildTransNodeNeighbor();
        void BuildCell();

    private:
        GridBlock *m_block_grid;
        GridFN *m_fn_grid;
        ModelManager *m_model_manager;
        IdxStruct *m_idx_proxy;
        std::vector<PhysicalType> m_cell_type;
        std::vector<PhysicalType> m_node_type;
        std::vector<TransFace> m_trans_face;
        // 保存参考节点, 与过渡节点连接的流体节点
        std::set<ConnectInfo> m_ref_node;
        // 保存过渡节点, 与物面节点连接的流体节点
        std::set<ConnectInfo> m_trans_node;
        FNGridInfo m_fn_info;

        int m_layer_num;
    };
}