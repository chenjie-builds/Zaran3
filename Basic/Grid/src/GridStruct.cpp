#include "GridStruct.h"
namespace zaran
{
	GridStruct::GridStruct(const string& name, int index, int dim) :GridBase(name, index, dim, GridType::Structured)
	{
		m_ghost_size = 1;
		m_node = nullptr;
		m_face = nullptr;
		m_cell = nullptr;
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

	}
	int GridStruct::GetNi()
	{
		return m_node->GetNi();
	}

	int GridStruct::GetNj()
	{
		return m_node->GetNj();
	}

	int GridStruct::GetNk()
	{
		return m_node->GetNk();
	}
	void GridStruct::GetNodeNum(int& ni, int& nj, int& nk)
	{
		ni = GetNi();
		nj = GetNj();
		nk = GetNk();
	}
    int GridStruct::GetTotalNodeNum()
    {
        return GetNi() * GetNj() * GetNk();
    }
	void GridStruct::GetRange(int& iStart, int& iEnd, int& jStart, int& jEnd, int& kStart, int& kEnd)
	{
		iStart = m_ghost_size;
		iEnd = GetNi() - m_ghost_size;
		jStart = m_ghost_size;
		jEnd = GetNj() - m_ghost_size;
		kStart = m_ghost_size;
		kEnd = GetNk() - m_ghost_size;
	}
}