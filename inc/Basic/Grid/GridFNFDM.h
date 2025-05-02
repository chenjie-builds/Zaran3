/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file GridFN.h
 * \brief Grid for Flexible Node FDM.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
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
		unique_ptr<NodeFN> m_node;
		unique_ptr<FaceFN> m_face;
		unique_ptr<CellFN> m_cell;
		unique_ptr<BoundManagerFN> m_boundary_map;
		int m_inner_node_num;
		int* m_inner_node_index;
		int m_bound_node_num;
		int* m_bound_node_index;
		int m_total_node_num;
	};
} // namespace zaran
