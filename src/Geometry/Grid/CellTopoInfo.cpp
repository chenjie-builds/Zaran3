#include "CellTopoInfo.h"
using namespace zaran;
CellTopoInfo::CellTopoInfo()
{

}

CellTopoInfo::~CellTopoInfo()
{

}

void CellTopoInfo::SetNode(IArray& node)
{
	node_ = node;

}

void CellTopoInfo::SetFace(IArray& face)
{
	face_ = face;
}

IArray& CellTopoInfo::GetFace()
{
	return face_;
}


IArray& CellTopoInfo::GetNode()
{
	return node_;
}
