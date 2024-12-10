#include "Residual.h"
#include"Log.h"
namespace zaran
{

	ResInfo::ResInfo(index_type equ_num) : m_equ_num(equ_num)
	{
		m_res_inf.resize(equ_num);
		m_ref_L2.resize(equ_num);
		m_res_inf_coord.resize(equ_num * 3);
		m_res_inf_idx.resize(equ_num);
	}
	ResInfo::~ResInfo()
	{
	}

	double ResInfo::GetInfNorm(index_type iVar)
	{
		return m_res_inf[iVar];
	}

	double ResInfo::GetL2Norm(index_type iVar)
	{
		return m_ref_L2[iVar];
	}

	const double* ResInfo::GetInfNormCoord(index_type iVar) const
	{
		return  &m_res_inf_coord[3 * iVar];
	}

	int ResInfo::GetInfNormIdx(index_type iVar) const
	{
		return m_res_inf_idx[iVar];
	}

	void ResInfo::SetInfNorm(index_type iVar, double value)
	{
		m_res_inf[iVar] = value;
	}

	void ResInfo::SetL2Norm(index_type iVar, double value)
	{
		m_ref_L2[iVar] = value;
	}

	void ResInfo::SetInfNormCoord(index_type iVar, const double* coord)
	{
		for (size_t i_dim = 0; i_dim < 3; i_dim++)
		{
			m_res_inf_coord[3 * iVar + i_dim] = coord[i_dim];
		}
	}

	void ResInfo::SetInfNormIdx(index_type iVar, index_type idx)
	{
		m_res_inf_idx[iVar] = idx;
	}

}