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
	class FaceTopoInfo
	{
	public:
		FaceTopoInfo();
		~FaceTopoInfo();
		void SetNodeIndex(Array<IArray>& node_index);
		void SetLeftCell(IArray& left_cell_index);
		void SetRightCell(IArray& right_cell_index);
		Array<IArray>& GetFace2Node() { return m_node_index; }
		IArray& GetLeftCellIndex() { return m_left_cell_index; }
		IArray& GetRightCellIndex() { return m_right_cell_index; }
	private:
		Array<IArray> m_node_index;
		IArray m_left_cell_index;
		IArray m_right_cell_index;
	};
}