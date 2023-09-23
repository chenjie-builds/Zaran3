#include"Grid_Struct_2D.h"
namespace zaran
{
	void Grid_Struct_2D::SetNi(int ni)
	{
		m_ni = ni;
	}
	void Grid_Struct_2D::SetNj(int nj)
	{
		m_nj = nj;
	}
	int Grid_Struct_2D::GetNi()
	{
		return m_ni;
	}
	int Grid_Struct_2D::GetNj()
	{
		return m_nj;
	}
	void Grid_Struct_2D::SetNodeNum(int ni, int nj)
	{
		SetNi(ni);
		SetNj(nj);
		SetTotalNodeNum(ni * nj);
	}
	void Grid_Struct_2D::GetNodeNum(int& ni, int& nj)
	{
		ni = GetNi();
		nj = GetNj();
	}
	int Grid_Struct_2D::GetNodeIndex(int i, int j)
	{
		return j * GetNi() + i;
	}
	void Grid_Struct_2D::GetNodeIndex(int index, int& i, int& j)
	{
		j = index / GetNi();
		i = index - j * GetNi();
	}
	int Grid_Struct_2D::GetCellIndex(int i, int j)
	{
		return (j - 1) * (GetNi() - 1) + i - 1;
	}
	void Grid_Struct_2D::GetCellIndex(int index, int& i, int& j)
	{
		j = index / (GetNi() - 1) + 1;
		i = index - (j - 1) * (GetNi() - 1) + 1;
	}
	void Grid_Struct_2D::GetRange(int& iStart, int& iEnd, int& jStart, int& jEnd)
	{
		iStart = 1;
		iEnd = GetNi() - 1;
		jStart = 1;
		jEnd = GetNj() - 1;
	}
}