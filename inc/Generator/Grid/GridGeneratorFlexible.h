/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file GridGeneratorFlexible.h
 * \brief GridGeneratorFlexible class, used to generate flexible grid.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "GridGenerator.h"
#include "GridFNFDM.h"

namespace zaran
{
	class GridGeneratorFlexible : public GridGenerator
	{
	public:
		GridGeneratorFlexible(const string& node_file_name, const string& ele_file_name, const string& bnd_file_name);
		void CreateGrid(dynamic_array<shared_ptr<GridBase>>& grid_list) override;
	private:
		void ReadNodeFile();
		void SortNeiborNode();
		void ReadCellFile();
		void ReadBoundFile();
		void CheckNode();
		void CheckUnkownNode();
		void CheckNeighborNum();
		void AddSelfToNeighbor();
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


}