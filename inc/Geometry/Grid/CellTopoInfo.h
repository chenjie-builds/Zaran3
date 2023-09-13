//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	CelTopoInfo.h														||
//*	@brief	Cell Topology setting												||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include"BasicType.h"
namespace zaran
{
	// 单元拓扑类，记录单元包含的节点和面元
	class CellTopoInfo
	{
	public:
		CellTopoInfo();
		~CellTopoInfo();
		void SetNodeIndex(Array<IArray>& node_index);
		void SetFaceIndex(Array<IArray>& face_index);
		Array<IArray>& GetNodeIndex();
		Array<IArray>& GetFaceIndex();
	private:
		// 单元包含的节点下标
		// m_node_index[id] = {node1,node2,node3,...} 表示单元id包含的节点
		Array<IArray> m_node_index;
		// 单元包含的面元
		// m_face_index[id] = {face1,face2,face3,...} 表示单元id包含的面元
		Array<IArray> m_face_index;
	};
}