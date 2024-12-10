#include "FaceFNFDM.h"
#include"Log.h"
namespace zaran
{
	FaceFN::FaceFN()
	{

	}

	FaceFN::~FaceFN()
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

	void FaceFN::Allocate(index_type nFace, index_type* face_node_num)
	{
		SetFaceNum(nFace);
		m_face_node_num = new index_type[nFace];
		m_node_id = new index_type[nFace];
		index_type nNode = 0;
		for (index_type i = 0;i < nFace;i++)
		{
			m_node_id[i] = nNode;
			m_face_node_num[i] = face_node_num[i];
			nNode += face_node_num[i];
		}
		m_face2node = new index_type[nNode];
		m_face2cell = new index_type[2 * nFace];
		m_area = new double[nFace];
		m_normal = new double[3 * nFace];
	}

	void FaceFN::SetFace2Node(index_type iFace, index_type* node, index_type nNode)
	{
		if (nNode != m_face_node_num[iFace])
		{
			Log::error("FaceTopoInfo::SetFace2Node: The number of nodes is not equal to the number of nodes in the face");
		}
		for (index_type i = 0; i < nNode; i++)
		{
			m_face2node[m_node_id[iFace] + i] = node[i];
		}
	}



}