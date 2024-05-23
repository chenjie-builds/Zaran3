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
#include"BasicType.h"
namespace zaran
{
	class Boundary
	{
	public:
		Boundary();
		Boundary(int nodeIndex, int innerIndex, int ghostIndex, const double* norm);
		~Boundary() {};
		void SetGhostIndex(const int& ghostIndex);
		void SetInnerIndex(const int& innnerIndex);
		void SetNode(const int& nodeIndex);
		void SetNorm(const DVector3D& boundNorm);
		int& GetIndex() { return m_bound_index; }
		int& GetInnerIndex() { return m_inner_index; }
		int& GetGhostIndex() { return m_ghost_index; }
		const double* GetNorm() { return m_norm; }
		bool operator==(const Boundary& bound);
	private:
		int m_bound_index;
		int m_inner_index;
		int m_ghost_index;
		double m_norm[3];
	};
}