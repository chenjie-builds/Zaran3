#include "Boundary.h"
using namespace zaran;
Boundary::Boundary()
{
	m_idx_bound = -1;
	m_idx_ghost = -1;
	m_idx_ref = -1;
	m_norm_bound[0]=m_norm_bound[1]=m_norm_bound[2]=0.0;
}

Boundary::Boundary(int idx_bound, int idx_ref, int idx_ghost,const double* norm_bound)
{
	m_idx_bound = idx_bound;
	m_idx_ref = idx_ref;
	m_idx_ghost = idx_ghost;
	for(int i=0;i<3;i++)
	{
		m_norm_bound[i] = norm_bound[i];
	}
}

void Boundary::SetIdxGhost(const int& idx_ghost)
{
	m_idx_ghost = idx_ghost;
}

void Boundary::SetIdxRef(const int& idx_ref)
{
	m_idx_ref = idx_ref;
}

void Boundary::SetIdxBound(const int& idx_bound)
{
	m_idx_bound = idx_bound;
}
void Boundary::SetNormBound(const DVector3D& boundNorm)
{
	for(int i=0;i<3;i++)
	{
		m_norm_bound[i] = boundNorm[i];
	}
}

bool Boundary::operator==(const Boundary& bound)
{
	return m_idx_bound == bound.m_idx_bound && m_idx_ghost == bound.m_idx_ghost && m_idx_ref == bound.m_idx_ref && m_norm_bound == bound.m_norm_bound;
}

