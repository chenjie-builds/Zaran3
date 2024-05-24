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
		Boundary(int idx_bound, int idx_ref, int idx_ghost, const double* norm_bound);
		~Boundary() {};
		void SetIdxGhost(const int& idx_ghost);
		void SetIdxRef(const int& innnerIndex);
		void SetIdxBound(const int& idx_bound);
		void SetNormBound(const DVector3D& boundNorm);
		int& GetIdxBound() { return m_idx_bound; }
		int& GetIdxRef() { return m_idx_ref; }
		int& GetIdxGhost() { return m_idx_ghost; }
		const double* GetNormBound() { return m_norm_bound; }
		bool operator==(const Boundary& bound);
	private:
		int m_idx_bound;
		int m_idx_ref;
		int m_idx_ghost;
		double m_norm_bound[3];
	};
}