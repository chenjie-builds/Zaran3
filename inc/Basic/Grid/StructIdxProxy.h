//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	.h													||
//*	@brief							||
//*	@author	Chen Jie.															||
//==============================================================================||
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
        IdProxyStruct() = default;
        IdProxyStruct(index_type ni, index_type nj, index_type nk);
        ~IdProxyStruct();
        void SetIdx(index_type i, index_type j, index_type k);
        void SetIdx(index_type idx);
        void GetIdxStruct(index_type&i, index_type&j, index_type&k) const ;
        using IdProxy::GetIdx;
        index_type GetIdx(index_type i, index_type j, index_type k) const
        {
            return i + j * m_ni + k * m_ni * m_nj;
        }
    private:
        index_type m_idx_i, m_idx_j, m_idx_k;
		// 网格的节点个数
        index_type m_ni, m_nj, m_nk;
    };
}