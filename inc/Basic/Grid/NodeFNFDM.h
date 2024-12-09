//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	NodeTopoInfo.h														||
//*	@brief	节点的拓扑结构														||
//*			包含节点类型，更新标记，邻居模板，坐标变换系数，面元单元连接信息			||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "BasicType.h"
#include "NodeBase.h"
namespace zaran
{
  enum class NodeType // 节点类型
  {
    inner = -1,      // 内部点，不是边界
    inlet = 1,       // 入口
    outlet = 2,      // 出口
    wall = 3,        // 壁面
    hole = 4,        // 洞边界
    symmetry = 5,    // 对称面
    farfield = 6,    // 远场
    corner = 7,      // 角点
    ghost = 8,       // ghonst point
    userDefined = 9, // 用户自定义
    undefined = -2   // 未定义
  };
  class NodeFN : public NodeBase
  {
  public:
    NodeFN(Id node_num);
    NodeFN(Id node_num, const Id*neighbor_node_num, const Id*neighbor_face_num,
           const Id*neighbor_cell_num);
    virtual ~NodeFN();
    // 设置节点坐标
    void SetCoord(Id idx, const double *coord);
    // 设置节点类型
    void SetType(Id idx, NodeType type);
    void SetNeighborNode(Id node_num, const Id*neighbor_node_num, const Id*neighbor_node);
    /// @brief 设置邻居节点
    /// @param idx 需要设置的节点编号
    /// @param neighbor_num 邻居节点个数
    /// @param neighbor_node_idx 邻居节点编号
    void SetNeighborNode(Id idx, Id neighbor_num, const Id*neighbor_node_idx);
    // 设置邻居面元
    void SetNeighborFace(Id idx, Id neighbor_num, const Id*neighbor_face);
    // 设置邻居单元
    void SetNeighborCell(Id idx, Id neighbor_num, const Id*neighbor_cell);

  public:
    // 返回节点类型
    const NodeType &GetType(const Id&index) const;
    const double *GetCoord(const Id&index) const;
    const Id*GetNeighborNode(const Id&idx) const;
    Id GetNeighborNodeNum(const Id&idx) const;
    const Id*GetNeighborFace(const Id&idx) const;
    Id GetNeighborFaceNum(const Id&idx) const;
    const Id*GetNeighborCell(const Id&idx) const;
    Id GetNeighborCellNum(const Id&idx) const;

  private:
    // 节点坐标
    double *m_coordinate;
    // 节点类型
    NodeType *m_type;
    /// @brief m_node_num[idx] 第iNode节点包含的邻居节点个数
    Id *m_node_num;
    /// @brief m_node[idx] 第iNode节点包含的邻居节点在m_neighor_node中的起始位置
    Id*m_node_offset;
    /// @brief m_node_idx[m_node[idx]] 第iNode节点的邻居节点
    /// @brief 前6个为差分模板, 剩下为其他邻居
    Id*m_node_idx;
    /// @brief m_face_num[idx] 第iNode节点包含的面元个数
    Id*m_face_num;
    /// @brief m_face[m_face_num[idx]] 第iNode节点包含的面元编号
    Id*m_face_idx;
    /// @brief m_neighor_face[m_face[idx]] 第iNode节点包含的面元的邻居面元
    Id*m_neighor_face;
    /// @brief m_cell_num[idx] 第iNode节点包含的单元个数
    Id*m_cell_num;
    /// @brief m_cell_index[m_cell_num[idx]] 第iNode节点包含的单元编号
    Id*m_neighor_cell_index;
    /// @brief m_neighor_cell[m_cell_index[idx]] 第iNode节点包含的单元的邻居单元
    Id*m_neighor_cell;
  };
} // namespace zaran