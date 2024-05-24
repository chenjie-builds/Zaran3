//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	FNFDM3D.h															||
//*	@brief	读取三维自由节点有限差分网格											||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "FaceFNFDM.h"
#include "GridBuilder.h"
#include "NodeFNFDM.h"
#include <string>
#include <unordered_map>
#include <vector>
#include"GridFNFDM.h"
namespace zaran
{
  class FNGridBuilder : public GridCreater
  {
  public:
    FNGridBuilder(const string& node_file_name = "node.dat", const string& ele_file_name = "cell.dat", const string& bnd_file_name = "bound.dat");
    GridFN* CreateGrid() override;
  private:
    void ReadNodeFile();
    void SortNeiborNode();
    void ExtendNeighborNode();
    void ReadCellFile();
    void ReadBoundFile();
    void CheckNode();
    void CheckUnkownNode();
    void AddSelfToNeighbor();
    void ConvertToGrid(GridFN*& grid);
  private:
    struct BoundNode
    {
      std::string type;
      int bound_index;
      int ref_index;
      double normal[3];
    };
    struct BoundFace
    {
      std::vector<int> face_node;
      std::vector<int> face_cell;
      std::vector<double> normal;
      double area;
    };


  private:
    string m_node_file_name;
    string m_ele_file_name;
    string m_bnd_file_name;
    std::vector<std::vector<double>> m_node_coord;
    std::vector<std::vector<int>> m_node_neibor;
    std::vector<NodeType> m_node_type;
    std::vector<BoundNode> m_bound_node;
    std::vector<std::vector<int>> m_cell_node;
    std::vector<BoundFace> m_bound_face;
  };
} // namespace zaran