#include "CellTopoInfo.h"
namespace zaran
{

	CellTopoInfo::CellTopoInfo()
	{

	}

	CellTopoInfo::~CellTopoInfo()
	{

	}


	void CellTopoInfo::SetNodeIndex(Array<IArray>& node_index)
	{
		m_node_index = node_index;
	}
	void CellTopoInfo::SetFaceIndex(Array<IArray>& face_index)
	{
		m_face_index = face_index;
	}
	Array<IArray>& CellTopoInfo::GetNodeIndex()
	{
		return m_node_index;
	}
	Array<IArray>& CellTopoInfo::GetFaceIndex()
	{
		return m_face_index;
	}
}