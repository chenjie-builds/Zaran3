//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	NodeInfoFNFDM.h														||
//*	@brief	自由节点有限差分节点信息														||
//*			包含节点类型，更新标记，邻居模板，坐标变换系数，面元单元连接信息		||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
namespace zaran
{
class NodeBase
{
  public:
    NodeBase();
    virtual ~NodeBase();
  public:
    // 返回节点总数
    const int &GetNodeNum() const;
  protected:
  void SetNodeNum(int node_num);
  private:
    int m_node_num;
};
} // namespace zaran