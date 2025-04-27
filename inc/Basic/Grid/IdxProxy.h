//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	IdProxy.h													||
//*	@brief							||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "BasicType.h"
namespace zaran
{
    /// @brief 索引代理类
    /// 用于不同网格类型的索引转换
    class IdProxy
    {
    public:
        IdProxy(index_type max_idx);
        index_type operator()(index_type idx) const;
    private:
        index_type m_max_idx;
    };
}