//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	FaceTopoInfo.h														||
//*	@brief	面元的拓扑结构,记录左右单元, 包含的节点									||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include"BasicType.h"
#include "FaceBase.h"
namespace zaran
{
	class FaceFN:public FaceBase
	{
	public:
		FaceFN();
		~FaceFN();
		void Allocate(Id nFace, Id* face_node_num);
		void SetLeftCell(Id iFace, Id iCell) { m_face2cell[2 * iFace] = iCell; }
		void SetRightCell(Id iFace, Id iCell) { m_face2cell[2 * iFace + 1] = iCell; }
		void SetArea(Id iFace, double area) { m_area[iFace] = area; }
		void SetNormal(Id iFace, double* normal) { m_normal[3 * iFace] = normal[0]; m_normal[3 * iFace + 1] = normal[1]; m_normal[3 * iFace + 2] = normal[2]; }
		void SetFace2Node(Id iFace, Id* node, Id nNode);
		Id GetFaceNodeNum(Id iFace) { return m_face_node_num[iFace]; }
		Id* GetFace2Node(Id iFace) { return m_face2node + m_node_id[iFace]; }
		Id GetLeftCell(Id iFace) { return m_face2cell[2 * iFace]; }
		Id GetRightCell(Id iFace) { return m_face2cell[2 * iFace + 1]; }
		double* GetNormal(Id iFace) { return m_normal + 3 * iFace; }
		double GetArea(Id iFace) { return m_area[iFace]; }
	private:
		// 面元包含的节点，所有的节点都记录在这里
		Id* m_face2node;
		// 面元包含的节点个数
		Id* m_face_node_num;
		// 面元记录节点起始位置在m_face2node中的索引
		Id* m_node_id;
		// 面元左右单元
		Id* m_face2cell;
		//面元面积
		double* m_area;
		//面元法向
		double* m_normal;

	};
}