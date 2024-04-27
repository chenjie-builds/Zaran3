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
	class Grid : public GridBase
	{
	public:
		Grid();
		Grid(const Grid& grid);
		~Grid();
	public:
		void SetLevel(const int& level);
		void SetInnerNodeNum(const int& innerNodeNum);
		void SetBoundNodeNum(const int& boundNodeNum);
		void SetTotalNodeNum(const int& totalNodeNum);
		void SetNodeTopo(NodeTopo* nodeTopo);
		void SetFaceTopo(FaceTopo* faceTopo);
		void SetCellTopo(CellTopo* cellTopo);
		void SetBoundaryMap(BoundaryMap* boundaryMap);
	public:
		const int& GetLevel()const;
		const int& GetTotalNodeNum()const;
		const int& GetInnerNodeNum()const;
		const int& GetBoundNodeNum()const;
		virtual CellTopo* GetCellTopo();
		virtual FaceTopo* GetFaceTopo();
		virtual NodeTopo* GetNodeTopo();
		BoundaryMap* GetBoundaryMap();
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
		NodeTopo* m_node_topo;
		// 单元拓扑信息
		CellTopo* m_cell_topo;
		// 面元拓扑信息
		FaceTopo* m_face_topo;
		// 边界信息
		BoundaryMap* m_bound_map;
	};
}