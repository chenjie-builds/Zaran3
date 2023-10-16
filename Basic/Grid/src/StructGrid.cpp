#include "Grid_Struct_3D.h"
namespace zaran
{
	void Grid_Struct_3D::SetNi(int ni)
	{
		m_ni = ni;
	}

	void Grid_Struct_3D::SetNj(int nj)
	{
		m_nj = nj;
	}

	void Grid_Struct_3D::SetNk(int nk)
	{
		m_nk = nk;
	}

	int Grid_Struct_3D::GetNi()
	{
		return m_ni;
	}

	int Grid_Struct_3D::GetNj()
	{
		return m_nj;
	}

	int Grid_Struct_3D::GetNk()
	{
		return m_nk;
	}
	void Grid_Struct_3D::SetNodeNum(int ni, int nj, int nk)
	{
		SetNi(ni);
		SetNj(nj);
		SetNk(nk);
		SetTotalNodeNum(ni * nj * nk);
	}
	void Grid_Struct_3D::GetNodeNum(int& ni, int& nj, int& nk)
	{
		ni = GetNi();
		nj = GetNj();
		nk = GetNk();
	}
	int Grid_Struct_3D::GetNodeIndex(int i, int j, int k)
	{
		return k * GetNi() * GetNj() + j * GetNi() + i;
	}
	void Grid_Struct_3D::GetNodeIndex(int index, int& i, int& j, int& k)
	{
		k = index / (GetNi() * GetNj());
		j = (index - k * GetNi() * GetNj()) / GetNi();
		i = index - k * GetNi() * GetNj() - j * GetNi();
	}
	int Grid_Struct_3D::GetCellIndex(int i, int j, int k)
	{
		return k * (GetNi() - 1) * (GetNj() - 1) + j * (GetNi() - 1) + i;
	}
	void Grid_Struct_3D::GetCellIndex(int index, int& i, int& j, int& k)
	{
		k = index / ((GetNi() - 1) * (GetNj() - 1));
		j = (index - k * (GetNi() - 1) * (GetNj() - 1)) / (GetNi() - 1);
		i = index - k * (GetNi() - 1) * (GetNj() - 1) - j * (GetNi() - 1);
	}
	void Grid_Struct_3D::GetRange(int& iStart, int& iEnd, int& jStart, int& jEnd, int& kStart, int& kEnd)
	{
		iStart = 1;
		iEnd = GetNi() - 1;
		jStart = 1;
		jEnd = GetNj() - 1;
		kStart = 1;
		kEnd = GetNk() - 1;
	}
}