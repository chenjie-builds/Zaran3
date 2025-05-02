/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file FaceFN.h
 * \brief Face Topology Class for FNFDM.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include"BasicType.h"
#include "FaceBase.h"
namespace zaran
{
	class FaceFN :public FaceBase
	{
	public:
		FaceFN();
		FaceFN(index_type face_num);
		~FaceFN();
		void SetFaceNum(index_type face_num);
		void Allocate(index_type nFace, index_type* face_node_num);
		void SetLeftCell(index_type iFace, index_type iCell) { m_face2cell[2 * iFace] = iCell; }
		void SetRightCell(index_type iFace, index_type iCell) { m_face2cell[2 * iFace + 1] = iCell; }
		void SetArea(index_type iFace, double area) { m_area[iFace] = area; }
		void SetNormal(index_type iFace, double* normal) { m_normal[3 * iFace] = normal[0]; m_normal[3 * iFace + 1] = normal[1]; m_normal[3 * iFace + 2] = normal[2]; }
		void SetFace2Node(index_type iFace, index_type* node, index_type nNode);
		index_type GetFaceNodeNum(index_type iFace) { return m_face_node_num[iFace]; }
		span<index_type> GetFace2Node(index_type iFace) { return span<index_type>(m_face2node.data() + m_node_id[iFace], m_face_node_num[iFace]); }
		index_type GetLeftCell(index_type iFace) { return m_face2cell[2 * iFace]; }
		index_type GetRightCell(index_type iFace) { return m_face2cell[2 * iFace + 1]; }
		span<double> GetNormal(index_type iFace) { return span<double>(m_normal.data() + 3 * iFace, 3); }
		double GetArea(index_type iFace) { return m_area[iFace]; }
	private:
		// 面元包含的节点，所有的节点都记录在这里
		dynamic_array<index_type> m_face2node;
		// 面元包含的节点个数
		dynamic_array<index_type> m_face_node_num;
		// 面元记录节点起始位置在m_face2node中的索引
		dynamic_array<index_type> m_node_id;
		// 面元左右单元
		dynamic_array<index_type> m_face2cell;
		//面元面积
		dynamic_array<double> m_area;
		//面元法向
		dynamic_array<double> m_normal;
	};
}