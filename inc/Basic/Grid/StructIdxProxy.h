/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file IdProxyStruct.h
 * \brief Index proxy class, used for index conversion of different grid types.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "BasicType.h"
#include "IdxProxy.h"
namespace zaran
{
    // 三维结构网格索引代理类
    // 用于计算三维结构网格的索引，便于和其他网格类型统一
    class IdProxyStruct : public IdProxy
    {
    public:
        IdProxyStruct(index_type max_i, index_type max_j, index_type max_k);
        index_type operator()(index_type i, index_type j, index_type k);
		void GetIdxStruct(index_type idx, index_type& i, index_type& j, index_type& k) const
		{
			k = idx / m_ij_plane_size;
			j = (idx - k * m_ij_plane_size) / m_ni;
			i = idx - k * m_ij_plane_size - j * m_ni;
		}
    private:
		// 网格的节点个数
        index_type m_ni, m_nj, m_nk;
        index_type m_ij_plane_size;
    };
}