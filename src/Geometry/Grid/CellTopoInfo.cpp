#include "CellTopoInfo.h"
namespace zaran
{
	CellTopoInfo::CellTopoInfo()
	{

	}

	CellTopoInfo::~CellTopoInfo()
	{

	}

	void CellTopoInfo::SetNode(std::vector<int>& node)
	{
		node_ = node;

	}

	void CellTopoInfo::SetFace(std::vector<int>& face)
	{
		face_ = face;
	}

	std::vector<int>& CellTopoInfo::GetFace()
	{
		return face_;
	}


	std::vector<int>& CellTopoInfo::GetNode()
	{
		return node_;
	}

}