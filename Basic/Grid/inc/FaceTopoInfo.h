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
namespace zaran
{
	class FaceTopo
	{
	public:
		FaceTopo();
		~FaceTopo();
		void Allocate(int nFace, int* face_node_num);
		void SetLeftCell(int iFace, int iCell) { m_face2cell[2 * iFace] = iCell; }
		void SetRightCell(int iFace, int iCell) { m_face2cell[2 * iFace + 1] = iCell; }
		void SetArea(int iFace, double area) { m_area[iFace] = area; }
		void SetNormal(int iFace, double* normal) { m_normal[3 * iFace] = normal[0]; m_normal[3 * iFace + 1] = normal[1]; m_normal[3 * iFace + 2] = normal[2]; }
		void SetFace2Node(int iFace, int* node, int nNode);
		int GetFaceNum() { return m_face_num; }
		int GetFaceNodeNum(int iFace) { return m_face_node_num[iFace]; }
		int* GetFace2Node(int iFace) { return m_face2node + m_node_id[iFace]; }
		int GetLeftCell(int iFace) { return m_face2cell[2 * iFace]; }
		int GetRightCell(int iFace) { return m_face2cell[2 * iFace + 1]; }
		double* GetNormal(int iFace) { return m_normal + 3 * iFace; }
		double GetArea(int iFace) { return m_area[iFace]; }
	private:
		// 面元个数
		int m_face_num;
		// 面元包含的节点，所有的节点都记录在这里
		int* m_face2node;
		// 面元包含的节点个数
		int* m_face_node_num;
		// 面元记录节点起始位置在m_face2node中的索引
		int* m_node_id;
		// 面元左右单元
		int* m_face2cell;
		//面元面积
		double* m_area;
		//面元法向
		double* m_normal;

	};
}