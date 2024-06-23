//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	NodeTopoInfoStruct.h												||
//*	@brief	结构网格节点拓扑信息, 还未完成, 后期可能删除							||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "NodeBase.h"
namespace zaran
{
class NodeStruct : public NodeBase
{

  public:
    NodeStruct();
    ~NodeStruct();
    void Allocate(int i_num, int j_num, int k_num);
    /// @brief 设置节点坐标
    /// @param idx_i i方向上的索引
    /// @param idx_j j方向上的索引
    /// @param idx_k k方向上的索引
    /// @param coord 输入的坐标，长度为3
    void SetCoord(int idx_i, int idx_j, int idx_k, const double *coord);
    /// @brief 获取节点坐标
    /// @param idx_i i方向上的索引
    /// @param idx_j j方向上的索引
    /// @param idx_k k方向上的索引
    /// @return 返回节点坐标，长度为3，不可修改
    const double *GetCoord(int idx_i, int idx_j, int idx_k) const;
    int GetINum() const { return m_i_num; }
    int GetJNum() const { return m_j_num; }
    int GetKNum() const { return m_k_num; }
  protected:
    /// @brief 将i,j,k转换为一维数组的索引，用于访问坐标数组
    /// @return 一维数组的索引
    int GetIdx(int i, int j, int k) const;
  protected:
    /// @brief 节点个数，分别为i,j,k方向上的节点个数
    int m_i_num, m_j_num, m_k_num;
    /// @brief 节点坐标，按照每个节点x,y,z的顺序存储
    double *m_coord;
};
} // namespace zaran