#include "StructGrid.h"
void zaran::StructGrid::SetNi(int ni)
{
	ni_ = ni;
}

void zaran::StructGrid::SetNj(int nj)
{
	nj_ = nj;
}

void zaran::StructGrid::SetNk(int nk)
{
	nk_ = nk;
}

int zaran::StructGrid::GetNi()
{
	return ni_;
}

int zaran::StructGrid::GetNj()
{
	return nj_;
}

int zaran::StructGrid::GetNk()
{
	return nk_;
}

void zaran::StructGrid::GetRange(int& iStart, int& iEnd, int& jStart, int& jEnd, int& kStart, int& kEnd)
{
	iStart = 1;
	iEnd = GetNi();
	jStart = 1;
	jEnd = GetNj();
	kStart = 1;
	kEnd = GetNk();
}

