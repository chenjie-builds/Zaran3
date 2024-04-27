#include "FaceTopoInfo.h"
#include"Log.h"
namespace zaran
{
	FaceTopo::FaceTopo()
	{

	}

	FaceTopo::~FaceTopo()
	{
		if (m_face2node)
			delete[] m_face2node;
		if (m_face_node_num)
			delete[] m_face_node_num;
		if (m_node_id)
			delete[] m_node_id;
		if (m_face2cell)
			delete[] m_face2cell;
		if (m_area)
			delete[] m_area;
		if (m_normal)
			delete[] m_normal;
	}

	void FaceTopo::Allocate(int nFace, int* face_node_num)
	{
		m_face_num = nFace;
		m_face_node_num = new int[nFace];
		m_node_id = new int[nFace];
		int nNode = 0;
		for (int i = 0;i < nFace;i++)
		{
			m_node_id[i] = nNode;
			m_face_node_num[i] = face_node_num[i];
			nNode += face_node_num[i];
		}
		m_face2node = new int[nNode];
		m_face2cell = new int[2 * nFace];
		m_area = new double[nFace];
		m_normal = new double[3 * nFace];
	}

	void FaceTopo::SetFace2Node(int iFace, int* node, int nNode)
	{
		if (nNode != m_face_node_num[iFace])
		{
			Log::error("FaceTopoInfo::SetFace2Node: The number of nodes is not equal to the number of nodes in the face");
		}
		for (int i = 0; i < nNode; i++)
		{
			m_face2node[m_node_id[iFace] + i] = node[i];
		}
	}



}