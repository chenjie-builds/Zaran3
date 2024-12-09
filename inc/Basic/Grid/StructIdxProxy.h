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
        IdProxyStruct(Id ni, Id nj, Id nk);
        ~IdProxyStruct();
        void SetIdx(Id i, Id j, Id k);
        void SetIdx(Id idx);
        void GetIdxStruct(Id&i, Id&j, Id&k) const ;
        using IdProxy::GetIdx;
        Id GetIdx(Id i, Id j, Id k) const
        {
            return i + j * m_ni + k * m_ni * m_nj;
        }
    private:
        Id m_idx_i, m_idx_j, m_idx_k;
		// 网格的节点个数
        Id m_ni, m_nj, m_nk;
    };
}