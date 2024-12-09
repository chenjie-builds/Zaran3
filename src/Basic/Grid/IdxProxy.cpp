#include "IdxProxy.h"
namespace zaran
{
	IdProxy::IdProxy(Id idx)
	{
		m_idx = idx;
	}

	IdProxy::~IdProxy()
	{
	}
	Id IdProxy::GetIdx() const
	{
		return m_idx;
	}
	void IdProxy::SetIdx(Id idx)
	{
		m_idx = idx;
	}
}