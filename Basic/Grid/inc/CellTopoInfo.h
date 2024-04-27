//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	CellTopoInfo.h														||
//*	@brief	Cell Topology setting												||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include"BasicType.h"
namespace zaran
{
	// 单元拓扑类，记录单元包含的节点和面元
	class CellTopo
	{
	public:
		CellTopo();
		~CellTopo();
		void SetNodeIndex(Array<IArray>& node_index);
		void SetFaceIndex(Array<IArray>& face_index);
		void SetCenterCoord(Array<DVector3D>& center);
		Array<IArray>& GetNodeIndex();
		Array<IArray>& GetFaceIndex();
		Array<DVector3D>& GetCenterCoord();
	private:
		// 单元包含的节点下标
		// m_node_index[id] = {node1,node2,node3,...} 表示单元id包含的节点
		Array<IArray> m_node_index;
		// 单元包含的面元
		// m_face_index[id] = {face1,face2,face3,...} 表示单元id包含的面元
		Array<IArray> m_face_index;
		// 单元中心坐标
		// m_center[id] = {x,y,z} 表示单元id的中心坐标
		Array<DVector3D> m_center;
	};
}