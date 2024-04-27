#pragma once
#include"grid.h"
#include "log.h"
namespace zaran
{

	Grid::Grid()
	{
		GridBase();
		m_node_topo = new NodeTopo();
		m_face_topo = new FaceTopo();
		m_cell_topo = new CellTopo();
		m_bound_map = new BoundaryMap();
	}

	Grid::Grid(const Grid& grid)
	{
	}

	Grid::~Grid()
	{
		delete[] m_node_topo;
		delete[] m_face_topo;
		delete[] m_cell_topo;
		delete[] m_bound_map;
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
	void Grid::SetNodeTopo(NodeTopo* nodeTopo)
	{
		m_node_topo = nodeTopo;
	}
	void Grid::SetFaceTopo(FaceTopo* faceTopo)
	{
		m_face_topo = faceTopo;
	}
	void Grid::SetCellTopo(CellTopo* cellTopo)
	{
		m_cell_topo = cellTopo;
	}
	void Grid::SetTotalNodeNum(const int& totalNodeNum)
	{
		m_num_total_node = totalNodeNum;
	}


	void Grid::SetBoundaryMap(BoundaryMap* boundaryMap)
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
	NodeTopo* Grid::GetNodeTopo()
	{
		return m_node_topo;
	}
	FaceTopo* Grid::GetFaceTopo()
	{
		return m_face_topo;
	}
	CellTopo* Grid::GetCellTopo()
	{
		return m_cell_topo;
	}

	BoundaryMap* Grid::GetBoundaryMap()
	{
		return m_bound_map;
	}
}