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
	using NodeTopoInfoVec = vector<NodeTopoInfo>;
	using CellTopoInfoVec = vector<CellTopoInfo>;
	using FaceTopoInfoVec = vector<FaceTopoInfo>;
	using BoundaryMapPtr = shared_ptr<BoundaryMap>;
	using InterNodeInfoPtr = shared_ptr<InterNodeInfo>;
	class Grid : public GridBase
	{
	public:
		Grid();
		Grid(const Grid& grid);
	public:
		void SetLevel(const int& level);
		void SetInnerNodeNumber(const int& innerNodeNum);
		void SetBoundNodeNumber(const int& boundNodeNum);
		void SetTotalNodeNumber(const int& totalNodeNum);
		void SetNodeTopo(NodeTopoInfoVec& nodeTopo);
		void SetFaceTopo(FaceTopoInfoVec& faceTopo);
		void SetCellTopoInfo(CellTopoInfoVec& cellTopo);
		void SetInterNodeInfo(shared_ptr<InterNodeInfo> interNodeInfo);
		void SetBoundaryMap(shared_ptr<BoundaryMap> boundaryMap);
	public:
		const int& GetLevel()const;
		const int& GetTotalNodeNum()const;
		const int& GetInnerNodeNum()const;
		const int& GetBoundNodeNum()const;
		CellTopoInfoVec& GetCellTopoInfo();
		FaceTopoInfoVec& GetFaceTopoInfo();
		NodeTopoInfoVec& GetNodeTopoInfo();
		InterNodeInfoPtr& GetInterNodeInfo();
		BoundaryMapPtr& GetBoundaryMap();
	private:
		/// 网格的级别
		int level_;
		/// 内部节点，正常计算的点
		int nInnerNode_;
		/// 边界节点数目
		int nBoundNode_;
		/// 节点总数
		int nTotalNode_;
		// 节点拓扑信息
		NodeTopoInfoVec nodeTopo_;
		// 单元拓扑信息
		CellTopoInfoVec cellTopo_;
		// 面元拓扑信息
		FaceTopoInfoVec faceTopo_;
		// 插值节点信息
		InterNodeInfoPtr  interNodeInfo_;
		// 边界信息
		BoundaryMapPtr boundaryMap_;
	};
}