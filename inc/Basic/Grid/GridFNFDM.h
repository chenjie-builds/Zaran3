//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	GridFNFDM.h															||
//*	@brief	自由节点有限差分网格													||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "GridBase.h"
#include "NodeFNFDM.h"
#include "FaceFNFDM.h"
#include "CellFNFDM.h"
#include"BoundFNManager.h"
namespace zaran
{
	class GridFN : public GridBase
	{
	public:
		GridFN(const string& name, int index, dimension_type dim);
		virtual ~GridFN();

	public:
		int GetTotalNodeNum() const;
		int GetInnerNodeNum() const;
		int GetBoundNodeNum() const;
		int* GetInnerNode();
		int* GetBoundNode();
		//get node
		NodeFN& GetNode() ;
		//get face
		FaceFN& GetFace() ;
		//get cell
		CellFN& GetCell() ;
		//get boundary map
		BoundManagerFN& GetBoundaryMap() ;

	protected:
		void InitNode();
	private:
		std::unique_ptr<NodeFN> m_node;
		std::unique_ptr<FaceFN> m_face;
		std::unique_ptr<CellFN> m_cell;
		std::unique_ptr<BoundManagerFN> m_boundary_map;
		int m_inner_node_num;
		int* m_inner_node_index;
		int m_bound_node_num;
		int* m_bound_node_index;
		int m_total_node_num;
	};
} // namespace zaran
