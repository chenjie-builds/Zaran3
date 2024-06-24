#include "GridStruct.h"
namespace zaran
{
	GridStruct::GridStruct(const string& name, int index, int dim) :GridBase(name, index, dim, GridType::Structured)
	{
		m_ghost_size = 1;
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
	void GridStruct::Allocate(int ni, int nj, int nk, int ghost_size)
	{
		m_ghost_size = ghost_size;
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
		m_node->Allocate(ni + ghost_size * 2, nj + ghost_size * 2, nk + ghost_size * 2);
		m_face = new FaceStruct();
		m_cell = new CellStruct();
		m_cell->Allocate(ni + ghost_size * 2 - 1, nj + ghost_size * 2 - 1, nk + ghost_size * 2 - 1);
		m_bound_map = new BoundMapStruct();
	}
	int GridStruct::GetNi()
	{
		return m_node->GetINum();
	}

	int GridStruct::GetNj()
	{
		return m_node->GetJNum();
	}

	int GridStruct::GetNk()
	{
		return m_node->GetKNum();
	}
	void GridStruct::GetNodeNum(int& ni, int& nj, int& nk)
	{
		ni = GetNi();
		nj = GetNj();
		nk = GetNk();
	}
	int GridStruct::GetTotalNodeNum()
	{
		return m_node->GetNodeNum();
	}
	void GridStruct::GetRange(int& iStart, int& iEnd, int& jStart, int& jEnd, int& kStart, int& kEnd)
	{
		iStart = m_ghost_size;
		iEnd = GetNi() - 2 * m_ghost_size;
		jStart = m_ghost_size;
		jEnd = GetNj() - 2 * m_ghost_size;
		kStart = m_ghost_size;
		kEnd = GetNk() - 2 * m_ghost_size;
	}
}