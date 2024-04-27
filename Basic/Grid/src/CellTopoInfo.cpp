#include "CellTopoInfo.h"
namespace zaran
{

	CellTopo::CellTopo()
	{
		m_center.resize(0);
		m_face_index.resize(0);
		m_node_index.resize(0);
	}

	CellTopo::~CellTopo()
	{

	}


	void CellTopo::SetNodeIndex(Array<IArray>& node_index)
	{
		m_node_index = node_index;
	}
	void CellTopo::SetFaceIndex(Array<IArray>& face_index)
	{
		m_face_index = face_index;
	}
	Array<IArray>& CellTopo::GetNodeIndex()
	{
		return m_node_index;
	}
	Array<IArray>& CellTopo::GetFaceIndex()
	{
		return m_face_index;
	}
	Array<DVector3D>& CellTopo::GetCenterCoord()
	{
		return m_center;
	}
}