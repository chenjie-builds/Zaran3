#include "IdxProxy.h"
namespace zaran
{
	IdProxy::IdProxy(index_type idx)
	{
		m_idx = idx;
	}

	IdProxy::~IdProxy()
	{
	}
	index_type IdProxy::GetIdx() const
	{
		return m_idx;
	}
	void IdProxy::SetIdx(index_type idx)
	{
		m_idx = idx;
	}
}