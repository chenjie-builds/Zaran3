//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Idx.h													||
//*	@brief							||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
namespace zaran
{
    /// @brief 索引代理类
    /// 用于不同网格类型的索引转换
    class Idx
    {
    public:
        Idx() = default;
        Idx(int idx);
        ~Idx();
        int GetIdx() const;
        void SetIdx(int idx);
    private:
        int m_idx; // 索引
    };
}