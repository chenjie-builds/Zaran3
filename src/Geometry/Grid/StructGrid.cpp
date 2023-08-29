#include "StructGrid.h"
using namespace zaran;
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

void StructGrid::GetRange(int& iStart, int& iEnd, int& jStart, int& jEnd, int& kStart, int& kEnd)
{
	iStart = 1;
	iEnd = GetNi();
	jStart = 1;
	jEnd = GetNj();
	kStart = 1;
	kEnd = GetNk();
}

