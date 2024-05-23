#include "Boundary.h"
using namespace zaran;
Boundary::Boundary()
{
	m_bound_index = -1;
	m_ghost_index = -1;
	m_inner_index = -1;
	m_norm[0]=m_norm[1]=m_norm[2]=0.0;
}

Boundary::Boundary(int boundIndex, int innerIndex, int ghostIndex,const double* norm)
{
	m_bound_index = boundIndex;
	m_inner_index = innerIndex;
	m_ghost_index = ghostIndex;
	for(int i=0;i<3;i++)
	{
		m_norm[i] = norm[i];
	}
}

void Boundary::SetGhostIndex(const int& ghostIndex)
{
	m_ghost_index = ghostIndex;
}

void Boundary::SetInnerIndex(const int& innerIndex)
{
	m_inner_index = innerIndex;
}

void Boundary::SetNode(const int& boundIndex)
{
	m_bound_index = boundIndex;
}
void Boundary::SetNorm(const DVector3D& boundNorm)
{
	for(int i=0;i<3;i++)
	{
		m_norm[i] = boundNorm[i];
	}
}

bool Boundary::operator==(const Boundary& bound)
{
	return m_bound_index == bound.m_bound_index && m_ghost_index == bound.m_ghost_index && m_inner_index == bound.m_inner_index && m_norm == bound.m_norm;
}

