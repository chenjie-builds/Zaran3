//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	GridFNFDM.h															||
//*	@brief	自由节点有限差分网格													||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "GridBase.h"
#include "NodeFNFDM.h"
#include "FaceFNFDM.h"
#include "CellFNFDM.h"
#include"BoundMapFN.h"
namespace zaran
{
  class GridFN : public GridBase
  {
  public:
    GridFN(const string& name, int index, int dim);
    virtual ~GridFN();

  public:
    void SetNode(NodeFN* node);
    void SetFace(FaceFN* face);
    void SetCell(CellFN* cell);
    void SetBoundaryMap(BoundMapFN* boundaryMap);
  public:
    int GetTotalNodeNum() const;
    int GetInnerNodeNum() const;
    int GetBoundNodeNum() const;
    NodeFN* GetNode();
    FaceFN* GetFace();
    CellFN* GetCell();
    BoundMapFN* GetBoundaryMap();
    int* GetInnerNode();
    int* GetBoundNode();
  protected:
    void InitNode();
  private:
    NodeFN* m_node;
    FaceFN* m_face;
    CellFN* m_cell;
    BoundMapFN* m_boundary_map;
    int m_inner_node_num;
    int* m_inner_node_index;
    int m_bound_node_num;
    int* m_bound_node_index;
    int m_total_node_num;
  };
} // namespace zaran
