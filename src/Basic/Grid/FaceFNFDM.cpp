#include "FaceFNFDM.h"
#include"Log.h"
namespace zaran
{
	FaceFN::FaceFN()
	{
		SetFaceNum(0);
	}

	FaceFN::FaceFN(index_type face_num)
	{
		SetFaceNum(face_num);
	}

	FaceFN::~FaceFN()
	{
	}

	void FaceFN::SetFaceNum(index_type face_num)
	{
		FaceBase::SetFaceNum(face_num);
		m_face_node_num.resize(face_num);
		m_node_id.resize(face_num);
		m_face2cell.resize(face_num * 2);
		m_area.resize(face_num);
		m_normal.resize(face_num * 3);
	}

	void FaceFN::Allocate(index_type face_num, index_type* face_node_num)
	{
		SetFaceNum(face_num);
		index_type nNode = 0;
		for (index_type i = 0; i < face_num; i++)
		{
			m_node_id[i] = nNode;
			m_face_node_num[i] = face_node_num[i];
			nNode += face_node_num[i];
		}
		m_face2node.resize(nNode);
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