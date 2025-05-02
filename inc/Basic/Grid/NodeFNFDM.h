/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file NodeFN.h
 * \brief Node information for FNFDM.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "BasicType.h"
#include "NodeBase.h"
#include "span.hpp"
namespace zaran
{
	enum class NodeType // 节点类型
	{
		inner = -1,      // 内部点，不是边界
		inlet = 1,       // 入口
		outlet = 2,      // 出口
		wall = 3,        // 壁面
		hole = 4,        // 洞边界
		symmetry = 5,    // 对称面
		farfield = 6,    // 远场
		corner = 7,      // 角点
		ghost = 8,       // ghonst point
		userDefined = 9, // 用户自定义
		undefined = -2   // 未定义
	};
	class NodeFN : public NodeBase
	{
	public:
		NodeFN();
		NodeFN(index_type node_num);
		NodeFN(index_type node_num, const index_type* neighbor_node_num, const index_type* neighbor_face_num,
			const index_type* neighbor_cell_num);
		void SetNodeNum(index_type node_num);
		//设置邻居信息:二维数组，第一维为节点编号，第二维为对应的邻居节点编号
		void SetNeighborNode(dynamic_array<dynamic_array<index_type>>& neighbor_node);
		void SetNeighborFace(dynamic_array<dynamic_array<index_type>>& neighbor_face);
		void SetNeighborCell(dynamic_array<dynamic_array<index_type>>& neighbor_cell);
		virtual ~NodeFN();
		// 设置节点坐标，第一维为节点编号，第二维为坐标（x,y,z）
		void SetCoord(const dynamic_array< dynamic_array<double>>& coord);
		void SetCoord(index_type idx, const double* coord);
		// 设置节点类型
		void SetType(const dynamic_array<NodeType>& type);
		void SetType(index_type idx, NodeType type);
		/// @brief 设置邻居节点
		/// @param idx 需要设置的节点编号
		/// @param neighbor_num 邻居节点个数
		/// @param neighbor_node_idx 邻居节点编号
		void SetNeighborNode(index_type idx, index_type neighbor_num, const index_type* neighbor_node_idx);
		// 设置邻居面元
		void SetNeighborFace(index_type idx, index_type neighbor_num, const index_type* neighbor_face);
		// 设置邻居单元
		void SetNeighborCell(index_type idx, index_type neighbor_num, const index_type* neighbor_cell);

	public:
		// 返回节点类型
		const NodeType& GetType(const index_type& index) const;
		span<double> GetCoord(const index_type& index);
		span<index_type> GetNeighborNode(index_type idx);
		index_type GetNeighborNodeNum(const index_type& idx) const;
		span<index_type> GetNeighborFace(index_type idx);
		index_type GetNeighborFaceNum(const index_type& idx) const;
		span<index_type> GetNeighborCell(index_type idx);
		index_type GetNeighborCellNum(const index_type& idx) const;
	private:
		// 节点坐标
		dynamic_array<double> m_coordinate;
		// 节点类型
		dynamic_array<NodeType> m_type;
		/// @brief m_node_num[idx] 第iNode节点包含的邻居节点个数
		dynamic_array<index_type> m_node_num;
		/// @brief m_node[idx] 第iNode节点包含的邻居节点在m_neighor_node中的起始位置
		dynamic_array<index_type> m_node_offset;
		/// @brief m_node_idx[m_node[idx]] 第iNode节点的邻居节点
		/// @brief 前6个为差分模板, 剩下为其他邻居
		dynamic_array<index_type> m_node_idx;
		/// @brief m_face_num[idx] 第iNode节点包含的面元个数
		dynamic_array<index_type> m_face_num;
		/// @brief m_face[m_face_num[idx]] 第iNode节点包含的面元编号
		dynamic_array<index_type> m_face_idx;
		/// @brief m_neighor_face[m_face[idx]] 第iNode节点包含的面元的邻居面元
		dynamic_array<index_type> m_neighor_face;
		/// @brief m_cell_num[idx] 第iNode节点包含的单元个数
		dynamic_array<index_type> m_cell_num;
		/// @brief m_cell_index[m_cell_num[idx]] 第iNode节点包含的单元编号
		dynamic_array<index_type> m_neighor_cell_index;
		/// @brief m_neighor_cell[m_cell_index[idx]] 第iNode节点包含的单元的邻居单元
		dynamic_array<index_type> m_neighor_cell;
	};
} // namespace zaran