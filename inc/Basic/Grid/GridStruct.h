//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Grid_Struct_3D.h													||
//*	@brief	三维结构网格类，根据结构网格数据结构特点存储信息						||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "GridBase.h"
#include "NodeStruct.h"
#include "FaceStruct.h"
#include "CellStruct.h"
#include "BoundMapStruct.h"
#include "StructIdxProxy.h"
#include <memory>
namespace zaran
{
	class GridStruct : public GridBase
	{
	public:
		GridStruct(const string &name, Id index, Id dim);
		virtual ~GridStruct();
		/// @brief 开辟内存空间
		/// @param ni i方向上的节点总数
		/// @param nj j方向上的节点总数
		/// @param nk k方向上的节点总数
		/// @param ghost_level ghost节点层数
		virtual void Allocate(Id ni, Id nj, Id nk, Id ghost_level);

	public:
		int GetNi() const;
		int GetNj() const;
		int GetNk() const;
		void GetNodeNum(Id& ni, Id& nj, Id& nk);
		int GetTotalNodeNum();
		std::shared_ptr<NodeStruct> GetNode() { return m_node; }
		std::shared_ptr<FaceStruct> GetFace() { return m_face; }
		std::shared_ptr<CellStruct> GetCell() { return m_cell; }
		std::shared_ptr<BoundManagerStruct> GetBoundMap() { return m_bound_map; }
		/// @brief 获取ghost节点层数
		/// @return  ghost节点层数
		int GetGhostLevel() const { return m_ghost_level; }
		std::shared_ptr<IdProxyStruct> GetIdxProxy() { return m_idx_proxy; }
	public:
		/// @brief 返回用于计算的节点范围,不包含ghost节点
		/// @param iStart i方向上的起始索引
		/// @param iEnd i方向上的结束索引
		/// @param jStart j方向上的起始索引
		/// @param jEnd j方向上的结束索引
		/// @param kStart k方向上的起始索引
		/// @param kEnd k方向上的结束索引
		void GetRange(int &iStart, int &iEnd, int &jStart, int &jEnd, int &kStart, int &kEnd);

	private:
		/// @brief ghost节点层数
		Id m_ghost_level;
		/// @brief 结构网格节点
		std::shared_ptr<NodeStruct> m_node;
		/// @brief 	结构网格面
		std::shared_ptr<FaceStruct> m_face;
		/// @brief 	结构网格单元
		std::shared_ptr<CellStruct>m_cell;	
		/// @brief 	边界条件
		std::shared_ptr<BoundManagerStruct> m_bound_map;
		std::shared_ptr<IdProxyStruct> m_idx_proxy;
	};
}