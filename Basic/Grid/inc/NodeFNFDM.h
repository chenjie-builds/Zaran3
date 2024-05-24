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
    inner,       // 内部点，不是边界
    inlet,       // 入口
    outlet,      // 出口
    wall,    // 壁面
    hole,        // 洞边界
    symmetry,    // 对称面
    farfield,    // 远场
    corner,      // 角点
    ghost,       // ghonst point
    userDefined, // 用户自定义
    undefined    // 未定义
};
class NodeFN:public NodeBase
{
  public:
    NodeFN(int node_num);
    NodeFN(int node_num, const int *neighbor_node_num, const int *neighbor_face_num,
             const int *neighbor_cell_num);
    virtual ~NodeFN();
    // 设置节点坐标
    void SetCoord(int iNode, const double *coord);
    // 设置节点类型
    void SetType(int iNode, NodeType &type);
    // 设置邻居点云
    void SetNeighborNode(int iNode, int neighbor_num, const int *neighbor_cloud);
    // 设置邻居面元
    void SetNeighborFace(int iNode, int neighbor_num, const int *neighbor_face);
    // 设置邻居单元
    void SetNeighborCell(int iNode, int neighbor_num, const int *neighbor_cell);
  public:
    // 返回节点类型
    const NodeType &GetType(const int &index) const;
    const double *GetCoord(const int &index) const;
    const int *GetNeighborNode(const int &iNode) const;
    int GetNeighborNodeNum(const int &iNode) const;
    const int *GetNeighborFace(const int &iNode) const;
    int GetNeighborFaceNum(const int &iNode) const;
    const int *GetNeighborCell(const int &iNode) const;
    int GetNeighborCellNum(const int &iNode) const;
  private:
    // 节点坐标
    double *m_coordinate;
    // 节点类型
    NodeType *m_type;
    /// @brief m_node_num[iNode] 第iNode节点包含的邻居节点个数
    int *m_neighor_node_num;
    /// @brief m_node[iNode] 第iNode节点包含的邻居节点在m_neighor_node中的起始位置
    int *m_neighor_node_index;
    /// @brief m_neighor_node[m_node[iNode]] 第iNode节点的邻居节点
    /// @brief 前6个为差分模板, 剩下为其他邻居
    int *m_neighor_node;
    /// @brief m_face_num[iNode] 第iNode节点包含的面元个数
    int *m_neighor_face_num;
    /// @brief m_face[m_face_num[iNode]] 第iNode节点包含的面元编号
    int *m_neighor_face_index;
    /// @brief m_neighor_face[m_face[iNode]] 第iNode节点包含的面元的邻居面元
    int *m_neighor_face;
    /// @brief m_cell_num[iNode] 第iNode节点包含的单元个数
    int *m_neighor_cell_num;
    /// @brief m_cell_index[m_cell_num[iNode]] 第iNode节点包含的单元编号
    int *m_neighor_cell_index;
    /// @brief m_neighor_cell[m_cell_index[iNode]] 第iNode节点包含的单元的邻居单元
    int *m_neighor_cell;
};
} // namespace zaran