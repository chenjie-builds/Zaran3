#include "BoundFN.h"
using namespace zaran;
BoundFN::BoundFN()
{
	m_idx_bound = -1;
	m_idx_ghost = -1;
	m_idx_ref = -1;
	m_norm_bound[0]=m_norm_bound[1]=m_norm_bound[2]=0.0;
}

BoundFN::BoundFN(int idx_bound, int idx_ref, int idx_ghost,const double* norm_bound)
{
	m_idx_bound = idx_bound;
	m_idx_ref = idx_ref;
	m_idx_ghost = idx_ghost;
	for(int i=0;i<3;i++)
	{
		m_norm_bound[i] = norm_bound[i];
	}
}

void BoundFN::SetIdxGhost(const int& idx_ghost)
{
	m_idx_ghost = idx_ghost;
}

void BoundFN::SetIdxRef(const int& idx_ref)
{
	m_idx_ref = idx_ref;
}

void BoundFN::SetIdxBound(const int& idx_bound)
{
	m_idx_bound = idx_bound;
}
void BoundFN::SetNorm(const double* boundNorm)
{
	for(int i=0;i<3;i++)
	{
		m_norm_bound[i] = boundNorm[i];
	}
}

bool BoundFN::operator==(const BoundFN& bound)
{
	return m_idx_bound == bound.m_idx_bound && m_idx_ghost == bound.m_idx_ghost && m_idx_ref == bound.m_idx_ref && m_norm_bound == bound.m_norm_bound;
}

