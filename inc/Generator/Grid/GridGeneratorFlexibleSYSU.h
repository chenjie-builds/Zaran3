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
    void CreateGrid(dynamic_array<shared_ptr<GridBase>>& grid_list) override;

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
    void ConvertToGrid(shared_ptr<GridFN>grid);

  private:
    struct BoundNode
    {
      std::string type;
      index_type bound_index;
      index_type ref_index;
      double normal[3];
    };
    struct BoundFace
    {
      dynamic_array<index_type> face_node;
      dynamic_array<index_type> face_cell;
      dynamic_array<double> normal;
      double area;
    };

  private:
    string m_node_file_name;
    string m_ele_file_name;
    string m_bnd_file_name;
    dynamic_array<dynamic_array<double>> m_node_coord;
    dynamic_array<dynamic_array<index_type>> m_node_neibor;
    dynamic_array<NodeType> m_node_type;
    dynamic_array<BoundNode> m_bound_node;
    dynamic_array<dynamic_array<index_type>> m_cell_node;
    dynamic_array<BoundFace> m_bound_face;
  };
} // namespace zaran