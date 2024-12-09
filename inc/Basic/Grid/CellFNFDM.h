//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	CellTopoInfo.h														||
//*	@brief	Cell Topology setting												||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include"BasicType.h"
#include "CellBase.h"
namespace zaran
{
	// 单元拓扑类，记录单元包含的节点和面元
	class CellFN:public CellBase
	{
	public:
		CellFN(int cell_num);
		~CellFN();
		void SetNode(std::vector<std::vector<Id>>& node_index);
		/// @brief 设置单元包含的面元
		/// @param face_index 面元索引
		void SetFace(std::vector<std::vector<Id>>& face_index);
		/// @brief 设置单元中心坐标
		/// @param iCell 单元编号
		/// @param center 中心坐标
		void SetCenter(Id iCell, const double* center);
		/// @brief 获取单元包含的节点个数
		Id GetNodeNum(Id iCell)const;
		/// @brief 获取单元包含的面元个数
		Id GetFaceNum(Id iCell)const;
		/// @brief 获取单元包含的节点
		const Id* GetNode(Id iCell);
		/// @brief 获取单元包含的面元
		const Id* GetFace(Id iCell);
		/// @brief 获取单元中心坐标
		const double* GetCenterCoord(Id iCell)const;
	private:
		/// @brief 单元包含的节点个数
		Id* m_node_num;
		/// @brief 记录单元包含节点在m_node中的索引
		Id* m_node_id;
		// 单元包含的节点下标
		// m_node[m_node_id[id]] = {node1,node2,node3,...} 表示单元id包含的节点
		Id* m_node;
		// 单元包含的面元个数
		Id* m_cell_face_num;
		// 记录单元包含面元在m_face中的索引
		Id* m_face_id;
		// 单元包含的面元
		// m_face[id] = {face1,face2,face3,...} 表示单元id包含的面元
		Id* m_face;
		// 单元中心坐标
		double* m_center;
	};
}