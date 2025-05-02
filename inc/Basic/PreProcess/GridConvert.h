/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file GridConvertPiflow2Flexible.h
 * \brief GridConvertPiflow2Flexible class, used to convert Piflow grid to flexible grid.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "BasicType.h"
#include <fstream>
namespace zaran
{
	/// @brief 将piflow的网格数据转换为Flexible的网格数据
	/// 输出为文件node.dat, cell.dat
	class GridConvertPiflow2Flexible
	{
	public:
		void SetFileName(const string& file_name) { m_file_name = file_name; }
		void ReadFile();
		void WriteFile();
	private:
		void ReadNodeCoord(std::ifstream& file);
		void ReadInnerNode(std::ifstream& file);
		void ReadGradNode(std::ifstream& file);
		void ReadBound(std::ifstream& file);
		void ReadCell(std::ifstream& file);
	private:
		void WriteNodeCoord(std::ofstream& file);
		void WriteInnerNode(std::ofstream& file);
		void WriteBound(std::ofstream& file);
		void WriteCell(std::ofstream& file);
	private:
		struct Coord
		{
			double x, y, z;
		};
		struct InnerNode
		{
			index_type idx;
			index_type neighbor[6];
		};
		struct GradNode
		{
			index_type idx;
			dynamic_array<index_type> neighbor;
		};
		struct Bound
		{
			index_type idx;
			index_type node_type;
			index_type ref_node;
			double bound_normal[3];
		};
		struct Cell
		{
			index_type node[8];
		};
	private:
		string m_file_name;
		dynamic_array<Coord> m_coord;
		dynamic_array<InnerNode> m_inner_node;
		dynamic_array<GradNode> m_grad_node;
		dynamic_array<Bound> m_bound;
		dynamic_array<Cell> m_cell;

	};
}