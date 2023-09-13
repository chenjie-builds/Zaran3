#pragma once
#include"grid.h"
#include "log.h"
namespace zaran
{

	Grid::Grid()
	{
		GridBase();
		m_node_topo = std::make_shared<NodeTopoInfo>();
		m_face_topo = std::make_shared<FaceTopoInfo>();
		m_cell_topo = std::make_shared<CellTopoInfo>();
		m_bound_map = std::make_shared<BoundaryMap>();
		m_inter_node_info = std::make_shared<InterNodeInfo>();
	}

	Grid::Grid(const Grid& grid)
	{

	}

	void Grid::SetLevel(const int& level)
	{
		m_level = level;
	}


	void Grid::SetInnerNodeNum(const int& innerNodeNum)
	{
		m_num_inner_node = innerNodeNum;
	}

	void Grid::SetBoundNodeNum(const int& boundNodeNum)
	{
		m_num_bound_node = boundNodeNum;
	}
	void Grid::SetNodeTopo(NodeTopoInfo& nodeTopo)
	{
		*m_node_topo = nodeTopo;
	}
	void Grid::SetFaceTopo(FaceTopoInfo& faceTopo)
	{
		*m_face_topo = faceTopo;
	}
	void Grid::SetCellTopo(CellTopoInfo& cellTopo)
	{
		*m_cell_topo = cellTopo;
	}
	void Grid::SetTotalNodeNum(const int& totalNodeNum)
	{
		m_num_total_node = totalNodeNum;
	}

	void Grid::SetInterNodeInfo(Ptr<InterNodeInfo> interNodeInfo)
	{
		m_inter_node_info = interNodeInfo;
	}


	void Grid::SetBoundaryMap(Ptr<BoundaryMap> boundaryMap)
	{
		m_bound_map = boundaryMap;
	}

	const int& Grid::GetLevel() const
	{
		return m_level;
	}


	const int& Grid::GetTotalNodeNum() const
	{
		return m_num_total_node;
	}


	const int& Grid::GetInnerNodeNum() const
	{
		return m_num_inner_node;
	}

	const int& Grid::GetBoundNodeNum() const
	{
		return m_num_bound_node;
	}
	Ptr<NodeTopoInfo>& Grid::GetNodeTopo()
	{
		return m_node_topo;
	}
	Ptr<FaceTopoInfo>& Grid::GetFaceTopo()
	{
		return m_face_topo;
	}
	Ptr<CellTopoInfo>& Grid::GetCellTopo()
	{
		return m_cell_topo;
	}
	InterNodeInfoPtr& Grid::GetInterNode()
	{
		return m_inter_node_info;
	}

	BoundaryMapPtr& Grid::GetBoundaryMap()
	{
		return m_bound_map;
	}
}