/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file GridStruct.h
 * \brief Struct Grid
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "GridBase.h"
#include "NodeStruct.h"
#include "FaceStruct.h"
#include "CellStruct.h"
#include "BoundStructManager.h"
#include "StructIdxProxy.h"
#include <memory>
namespace zaran
{
	class GridStruct : public GridBase
	{
	public:
		GridStruct(const string& name, index_type index, index_type dim);
		virtual ~GridStruct();
		/// @brief 开辟内存空间
		/// @param ni i方向上的节点总数
		/// @param nj j方向上的节点总数
		/// @param nk k方向上的节点总数
		/// @param ghost_level ghost节点层数
		virtual void Allocate(count_type ni, count_type nj, count_type nk, index_type ghost_level);

	public:
		count_type GetNi() const;
		count_type GetNj() const;
		count_type GetNk() const;
		void GetNodeNum(count_type& ni, count_type& nj, count_type& nk) const;
		int GetTotalNodeNum();
		NodeStruct* GetNode() { return m_node.get(); }
		FaceStruct* GetFace() { return m_face.get(); }
		CellStruct* GetCell() { return m_cell.get(); }
		BoundStructManager* GetBoundMap() { return m_bound_map.get(); }
		/// @brief 获取ghost节点层数
		/// @return  ghost节点层数
		int GetGhostLevel() const { return m_ghost_level; }
		IdProxyStruct& GetIdxProxy() { return *m_idx_proxy; }
	public:
		/// @brief 返回用于计算的节点范围,不包含ghost节点
		/// @param iStart i方向上的起始索引
		/// @param iEnd i方向上的结束索引
		/// @param jStart j方向上的起始索引
		/// @param jEnd j方向上的结束索引
		/// @param kStart k方向上的起始索引
		/// @param kEnd k方向上的结束索引
		void GetRange(index_type& iStart, index_type& iEnd, index_type& jStart, index_type& jEnd, index_type& kStart, index_type& kEnd) const;

	private:
		/// @brief ghost节点层数
		index_type m_ghost_level;
		/// @brief 结构网格节点
		unique_ptr<NodeStruct> m_node;
		/// @brief 	结构网格面
		unique_ptr<FaceStruct> m_face;
		/// @brief 	结构网格单元
		unique_ptr<CellStruct>m_cell;
		/// @brief 	边界条件
		unique_ptr<BoundStructManager> m_bound_map;
		shared_ptr<IdProxyStruct> m_idx_proxy;
	};
}