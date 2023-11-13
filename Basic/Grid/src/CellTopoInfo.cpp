#include "CellTopoInfo.h"
namespace zaran
{

	CellTopoInfo::CellTopoInfo()
	{
		m_center.resize(0);
		m_face_index.resize(0);
		m_node_index.resize(0);
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
	Array<DVector3D>& CellTopoInfo::GetCenterCoord()
	{
		return m_center;
	}
}