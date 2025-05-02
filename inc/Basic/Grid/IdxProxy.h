/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file IdProxy.h
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