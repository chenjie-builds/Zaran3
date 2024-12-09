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
#include "GridGenerator.h"
#include "GridFNFDM.h"
#include <string>
#include <unordered_map>
#include <vector>
namespace zaran
{
  /// @brief 读取中山大学网格数据
  /// @details 该数据分为三个文件，分别是node.dat, cell.dat, bound.dat
  /// @note 具体文件格式在文档中有详细说明
  class GridBuilderSYSU_FN : public GridGenerator
  {
  public:
    GridBuilderSYSU_FN(const string &node_file_name, const string &ele_file_name, const string &bnd_file_name);
    void CreateGrid(Array<std::shared_ptr<GridBase>>& grid_list) override;

  private:
    void ReadNodeFile();
    void SortNeiborNode();
    void ExtendNeighborNode();
    void ReadCellFile();
    void ReadBoundFile();
    void CheckNode();
    void CheckUnkownNode();
    void CheckNeighborNum();
    void AddSelfToNeighbor();
    // 将内部节点的邻居节点加入到边界节点的邻居节点中
    void SetBoundNeighbor();
    void ConvertToGrid(std::shared_ptr<GridFN>grid);

  private:
    struct BoundNode
    {
      std::string type;
      Id bound_index;
      Id ref_index;
      double normal[3];
    };
    struct BoundFace
    {
      std::vector<Id> face_node;
      std::vector<Id> face_cell;
      std::vector<double> normal;
      double area;
    };

  private:
    string m_node_file_name;
    string m_ele_file_name;
    string m_bnd_file_name;
    std::vector<std::vector<double>> m_node_coord;
    std::vector<std::vector<Id>> m_node_neibor;
    std::vector<NodeType> m_node_type;
    std::vector<BoundNode> m_bound_node;
    std::vector<std::vector<Id>> m_cell_node;
    std::vector<BoundFace> m_bound_face;
  };
} // namespace zaran