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
        IdProxy() = default;
        IdProxy(Id idx);
        ~IdProxy();
        Id GetIdx() const;
        void SetIdx(Id idx);
    private:
        Id m_idx; // 索引
    };
}