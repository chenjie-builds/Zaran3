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
#include"BoundaryMap.h"
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
    void SetBoundaryMap(BoundaryMap* boundaryMap);
  public:
    int GetTotalNodeNum() const;
    int GetInnerNodeNum() const;
    int GetBoundNodeNum() const;
    NodeFN* GetNodeTopo();
    FaceFN* GetFaceTopo();
    CellFN* GetCellTopo();
    BoundaryMap* GetBoundaryMap();
  private:
    NodeFN* m_node;
    FaceFN* m_face;
    CellFN* m_cell;
    BoundaryMap* m_boundary_map;
  };
} // namespace zaran
