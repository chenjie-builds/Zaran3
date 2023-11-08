//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Grid.h																||
//*	@brief	网格类, 存储点线面几何信息												||
//*	@author	Chen Jie.															||
//==============================================================================||

#pragma once
#include "GridBase.h"
#include "InterNodeInfo.h"
#include "NodeTopoInfo.h"
#include "FaceTopoInfo.h"
#include "CellTopoInfo.h"
#include "BoundaryMap.h"
#include "CommonPara.h"
#include <Eigen/Dense>
#include"BasicType.h"
// grid
namespace zaran
{
	using BoundaryMapPtr = Ptr<BoundaryMap>;
	using InterNodeInfoPtr = Ptr<InterNodeInfo>;
	class Grid : public GridBase
	{
	public:
		Grid();
		Grid(const Grid& grid);
	public:
		void SetLevel(const int& level);
		void SetInnerNodeNum(const int& innerNodeNum);
		void SetBoundNodeNum(const int& boundNodeNum);
		void SetTotalNodeNum(const int& totalNodeNum);
		void SetNodeTopo(NodeTopoInfo& nodeTopo);
		void SetFaceTopo(FaceTopoInfo& faceTopo);
		void SetCellTopo(CellTopoInfo& cellTopo);
		void SetInterNodeInfo(Ptr<InterNodeInfo> interNodeInfo);
		void SetBoundaryMap(Ptr<BoundaryMap> boundaryMap);
	public:
		const int& GetLevel()const;
		const int& GetTotalNodeNum()const;
		const int& GetInnerNodeNum()const;
		const int& GetBoundNodeNum()const;
		Ptr<CellTopoInfo>& GetCellTopo();
		Ptr<FaceTopoInfo>& GetFaceTopo();
		Ptr<NodeTopoInfo>& GetNodeTopo();
		InterNodeInfoPtr& GetInterNode();
		BoundaryMapPtr& GetBoundaryMap();
	protected:
		/// 网格的级别
		int m_level;
		/// 内部节点，正常计算的点
		int m_num_inner_node;
		/// 边界节点数目
		int m_num_bound_node;
		/// 节点总数
		int m_num_total_node;
		// 节点拓扑信息
		Ptr<NodeTopoInfo> m_node_topo;
		// 单元拓扑信息
		Ptr<CellTopoInfo> m_cell_topo;
		// 面元拓扑信息
		Ptr<FaceTopoInfo> m_face_topo;
		// 插值节点信息
		InterNodeInfoPtr  m_inter_node_info;
		// 边界信息
		BoundaryMapPtr m_bound_map;
	};
}