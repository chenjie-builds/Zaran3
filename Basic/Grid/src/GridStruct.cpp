#include "GridStruct.h"
namespace zaran
{
	GridStruct::GridStruct(const string &name, int index, int dim) : GridBase(name, index, dim, GridType::Structured)
	{
		m_ghost_level = 1;
		m_node = nullptr;
		m_face = nullptr;
		m_cell = nullptr;
		m_bound_map = nullptr;
	}
	GridStruct::~GridStruct()
	{
		if (m_node)
		{
			delete[] m_node;
			m_node = nullptr;
		}
		if (m_face)
		{
			delete[] m_face;
			m_face = nullptr;
		}
		if (m_cell)
		{
			delete[] m_cell;
			m_cell = nullptr;
		}
		if (m_bound_map)
		{
			delete[] m_bound_map;
			m_bound_map = nullptr;
		}
	}
	void GridStruct::Allocate(int ni, int nj, int nk, int ghost_level)
	{
		m_ghost_level = ghost_level;
		if (m_node != nullptr)
		{
			delete m_node;
			m_node = nullptr;
		}
		if (m_face != nullptr)
		{
			delete m_face;
			m_face = nullptr;
		}
		if (m_cell != nullptr)
		{
			delete m_cell;
			m_cell = nullptr;
		}
		if (m_bound_map != nullptr)
		{
			delete m_bound_map;
			m_bound_map = nullptr;
		}
		m_node = new NodeStruct();
		m_node->Allocate(ni + ghost_level * 2, nj + ghost_level * 2, nk + ghost_level * 2);
		m_face = new FaceStruct();
		m_cell = new CellStruct();
		m_cell->Allocate(ni + ghost_level * 2 - 1, nj + ghost_level * 2 - 1, nk + ghost_level * 2 - 1);
		m_bound_map = new BoundManagerStruct();
	}
	int GridStruct::GetNi()const
	{
		return m_node->GetINum();
	}

	int GridStruct::GetNj()const
	{
		return m_node->GetJNum();
	}

	int GridStruct::GetNk()const
	{
		return m_node->GetKNum();
	}
	void GridStruct::GetNodeNum(int &ni, int &nj, int &nk)
	{
		ni = GetNi();
		nj = GetNj();
		nk = GetNk();
	}
	int GridStruct::GetTotalNodeNum()
	{
		return m_node->GetCount();
	}
	void GridStruct::GetRange(int &is, int &ie, int &js, int &je, int &ks, int &ke)
	{
		is = GetGhostLevel();
		ie = GetNi() - GetGhostLevel() - 1;
		js = GetGhostLevel();
		je = GetNj() - GetGhostLevel() - 1;
		ks = GetGhostLevel();
		ke = GetNk() - GetGhostLevel() - 1;
	}
}