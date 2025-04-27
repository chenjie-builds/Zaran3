#include "IdxProxy.h"
namespace zaran
{
	IdProxy::IdProxy(index_type max_idx)
	{
		m_max_idx = max_idx;
	}

	index_type IdProxy::operator()(index_type idx) const
	{
#if DEBUG_MODE
		if (idx < 0 || idx >= m_max_idx)
		{
			Log::error("Invalid index:{}", idx);
			Log::error("Max index:{}", m_max_idx);
		}
#endif
		return idx;
	}


}