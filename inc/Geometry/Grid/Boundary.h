//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Boundary.h															||
//*	@brief	Describe Boundary													||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include <vector>
#include <Eigen/Dense>
namespace zaran
{
	class Boundary
	{
	public:
		Boundary();
		Boundary(int nodeIndex, int innerNodeIndex, int ghostNodeIndex, Eigen::Vector3d& norm);
		~Boundary() {};
		void SetGhostNodeIndex(const int& ghostNodeIndex);
		void SetInnerNodeIndex(const int& innnerNodeIndex);
		void SetNodeIndex(const int& nodeIndex);
		void SetNorm(const Eigen::Vector3d& boundNorm);
		 int& GetIndex() { return nodeIndex_; }
		 int& GetInnerNodeIndex() { return innerNodeIndex_; }
		 int& GetGhostNodeIndex() { return ghostNodeIndex_; }
		 Eigen::Vector3d& GetNorm() { return norm_; }
		 bool operator==(const Boundary& bound);
	private:
		int nodeIndex_;
		int innerNodeIndex_;
		int ghostNodeIndex_;
		Eigen::Vector3d norm_;
	};
}