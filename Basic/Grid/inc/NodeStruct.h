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
    NodeStruct(int ni, int nj, int nk);
    ~NodeStruct();
    void SetCoord(int i, int j, int k, const double *coord);
    const double *GetCoord(int i, int j, int k) const;
  protected:
    int GetIndex(int i, int j, int k) const;
  protected:
    int m_ni, m_nj, m_nk;
    double *m_coord;
};
} // namespace zaran