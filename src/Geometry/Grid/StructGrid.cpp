#include "StructGrid.h"
namespace zaran
{

	void StructGrid::SetNi(int ni)
	{
		ni_ = ni;
	}

	void StructGrid::SetNj(int nj)
	{
		nj_ = nj;
	}

	void StructGrid::SetNk(int nk)
	{
		nk_ = nk;
	}

	int StructGrid::GetNi()
	{
		return ni_;
	}

	int StructGrid::GetNj()
	{
		return nj_;
	}

	int StructGrid::GetNk()
	{
		return nk_;
	}
	void StructGrid::SetNodeNum(int ni, int nj, int nk)
	{
		SetNi(ni);
		SetNj(nj);
		SetNk(nk);
	}
	void StructGrid::GetNodeNum(int& ni, int& nj, int& nk)
	{
		ni = GetNi();
		nj = GetNj();
		nk = GetNk();
	}
	int StructGrid::GetNodeIndex(int i, int j, int k)
	{
		return k * GetNi() * GetNj() + j * GetNi() + i;
	}
	void StructGrid::GetNodeIndex(int index, int& i, int& j, int& k)
	{
		k = index / (GetNi() * GetNj());
		j = (index - k * GetNi() * GetNj()) / GetNi();
		i = index - k * GetNi() * GetNj() - j * GetNi();
	}
	int StructGrid::GetCellIndex(int i, int j, int k)
	{
		return k * (GetNi() - 1) * (GetNj() - 1) + j * (GetNi() - 1) + i;
	}
	void StructGrid::GetCellIndex(int index, int& i, int& j, int& k)
	{
		k = index / ((GetNi() - 1) * (GetNj() - 1));
		j = (index - k * (GetNi() - 1) * (GetNj() - 1)) / (GetNi() - 1);
		i = index - k * (GetNi() - 1) * (GetNj() - 1) - j * (GetNi() - 1);
	}
	void StructGrid::GetRange(int& iStart, int& iEnd, int& jStart, int& jEnd, int& kStart, int& kEnd)
	{
		iStart = 1;
		iEnd = GetNi() - 1;
		jStart = 1;
		jEnd = GetNj() - 1;
		kStart = 1;
		kEnd = GetNk() - 1;
	}
}