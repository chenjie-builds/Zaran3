#include "GridStruct.h"
namespace zaran
{
    GridStruct::GridStruct(const string& name, int index, int dim, GridType type):GridBase(name, index, dim, type)
    {
    }
    GridStruct::~GridStruct()
    {
		if (m_node)
		{
			delete[] m_node;
			m_node = nullptr;
		}
		if(m_face)
		{
			delete[] m_face;
			m_face = nullptr;
		}
		if(m_cell)
		{
			delete[] m_cell;
			m_cell = nullptr;
		}
    }
	int GridStruct::GetNi()
	{
		return m_ni;
	}

	int GridStruct::GetNj()
	{
		return m_nj;
	}

	int GridStruct::GetNk()
	{
		return m_nk;
	}
	void GridStruct::SetNodeNum(int ni, int nj, int nk)
	{
		m_ni = ni;
		m_nj = nj;
		m_nk = nk;
		m_node = new NodeStruct(ni, nj, nk);
	}
	void GridStruct::GetNodeNum(int& ni, int& nj, int& nk)
	{
		ni = GetNi();
		nj = GetNj();
		nk = GetNk();
	}
	void GridStruct::GetRange(int& iStart, int& iEnd, int& jStart, int& jEnd, int& kStart, int& kEnd)
	{
		iStart = 1;
		iEnd = GetNi() - 1;
		jStart = 1;
		jEnd = GetNj() - 1;
		kStart = 1;
		kEnd = GetNk() - 1;
	}
}