#include "BoundStruct.h"
namespace zaran
{

	BoundStruct::BoundStruct(int i_bnd, int j_bnd, int k_bnd, const int* direction, const double* norm, int i_bnd_tgt, int j_bnd_tgt, int k_bnd_tgt, const int* direction_tgt, int target_block) : m_i_bnd_src(i_bnd), m_j_bnd_src(j_bnd), m_k_bnd_src(k_bnd), m_i_bnd_tgt(i_bnd_tgt), m_j_bnd_tgt(j_bnd_tgt), m_k_bnd_tgt(k_bnd_tgt), m_target_block(target_block)
	{
		for (int i = 0; i < 3; i++)
		{
			m_dir_src[i] = direction[i];
			m_norm[i] = norm[i];
			if (direction_tgt != nullptr)
			{
				m_dir_tgt[i] = direction_tgt[i];
			}
		}
	}
	BoundStruct::BoundStruct()
	{
		for (int i = 0; i < 3; i++)
		{
			m_i_bnd_src = 0;
			m_j_bnd_src = 0;
			m_k_bnd_src = 0;
			m_i_bnd_tgt = 0;
			m_j_bnd_tgt = 0;
			m_k_bnd_tgt = 0;
			m_dir_src[i] = 0;
			m_norm[i] = 0;
			m_dir_tgt[i] = 0;


		}
	}
	BoundStruct::~BoundStruct()
	{
	}

	void BoundStruct::SetIdx(const int& i_bound, const int& j_bound, const int& k_bound)
	{
		m_i_bnd_src = i_bound;
		m_j_bnd_src = j_bound;
		m_k_bnd_src = k_bound;
	}

	void BoundStruct::SetDirection(const int* direction)
	{
		for (int i = 0; i < 3; i++)
		{
			m_dir_src[i] = direction[i];
		}
	}

	void BoundStruct::SetNorm(const double* norm_bound)
	{
		for (int i = 0; i < 3; i++)
		{
			m_norm[i] = norm_bound[i];
		}
	}

	void BoundStruct::GetIdx(int& i, int& j, int& k)
	{
		i = m_i_bnd_src;
		j = m_j_bnd_src;
		k = m_k_bnd_src;
	}
	const int* BoundStruct::GetDirection()
	{
		return m_dir_src;
	}
	const double* BoundStruct::GetNorm()
	{
		return m_norm;
	}
	bool BoundStruct::operator==(const BoundStruct& bound)
	{
		return m_i_bnd_src == bound.m_i_bnd_src && m_j_bnd_src == bound.m_j_bnd_src && m_k_bnd_src == bound.m_k_bnd_src && m_dir_src == bound.m_dir_src && m_norm == bound.m_norm && m_i_bnd_tgt == bound.m_i_bnd_tgt && m_j_bnd_tgt == bound.m_j_bnd_tgt && m_k_bnd_tgt == bound.m_k_bnd_tgt && m_dir_tgt == bound.m_dir_tgt;
	}
	void BoundStruct::SetIdxTgt(const int& i_bnd_tgt, const int& j_bnd_tgt, const int& k_bnd_tgt)
	{
		m_i_bnd_tgt = i_bnd_tgt;
		m_j_bnd_tgt = j_bnd_tgt;
		m_k_bnd_tgt = k_bnd_tgt;
	}
	void BoundStruct::SetDirectionTgt(const int* direction_tgt)
	{
		for (int i = 0; i < 3; i++)
		{
			m_dir_tgt[i] = direction_tgt[i];
		}
	}
	void BoundStruct::GetIdxTgt(int& i_bnd_tgt, int& j_bnd_tgt, int& k_bnd_tgt)
	{
		i_bnd_tgt = m_i_bnd_tgt;
		j_bnd_tgt = m_j_bnd_tgt;
		k_bnd_tgt = m_k_bnd_tgt;
	}
	const int* BoundStruct::GetDirectionTgt()
	{
		return m_dir_tgt;
	}
	void BoundStruct::SetTargetBlock(const int& target_block)
	{
		m_target_block = target_block;
	}
	int BoundStruct::GetTargetBlock()
	{
		return m_target_block;
	}
}