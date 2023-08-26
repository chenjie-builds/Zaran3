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
#include <vector>
namespace zaran
{
	class FaceTopoInfo
	{
	public:
		FaceTopoInfo();
		~FaceTopoInfo();
		void SetFace2Node(std::vector<int>& face2node);
		void SetLeftCell(int& leftCell);
		void SetRightCell(int& rightCell);
		std::vector<int>& GetFace2Node();
		int& GetLeftCell();
		int& GetRightCell();
	private:
		std::vector<int>face2node_;
		int leftCell_;
		int rightCell_;
	};
}