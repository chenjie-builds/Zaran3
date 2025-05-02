/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file BoundaryFN.h
 * \brief Boundary Information Manager for Flexible Node FDM.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include"BasicType.h"
namespace zaran
{
	class BoundFN
	{
	public:
		BoundFN();
		BoundFN(int idx_bound, int idx_ref, int idx_ghost, const double* norm_bound);
		~BoundFN() {};
		void SetIdxGhost(const int& idx_ghost);
		void SetIdxRef(const int& innnerIndex);
		void SetIdxBound(const int& idx_bound);
		void SetNorm(const double* boundNorm);
		int& GetIdxBound() { return m_idx_bound; }
		int& GetIdxRef() { return m_idx_ref; }
		int& GetIdxGhost() { return m_idx_ghost; }
		const double* GetNormBound() { return m_norm_bound; }
		bool operator==(const BoundFN& bound);
	private:
		int m_idx_bound;
		int m_idx_ref;
		int m_idx_ghost;
		double m_norm_bound[3];
	};
}