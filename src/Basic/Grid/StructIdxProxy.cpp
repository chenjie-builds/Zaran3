#include "StructIdxProxy.h"  
namespace zaran
{
	IdProxyStruct::IdProxyStruct(index_type max_i, index_type max_j, index_type max_k)
		: IdProxy(max_i* max_j* max_k) // Explicitly call the base class constructor  
	{
		m_ni = max_i;
		m_nj = max_j;
		m_nk = max_k;
		m_ij_plane_size = max_i * max_j;
	}

	index_type IdProxyStruct::operator()(index_type i, index_type j, index_type k)
	{
		return i + j * m_ni + k * m_ij_plane_size;
	}

}