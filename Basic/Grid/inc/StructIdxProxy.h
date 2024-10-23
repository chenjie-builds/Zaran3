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
#include "IdxProxy.h"
namespace zaran
{
    // 三维结构网格索引代理类
    // 用于计算三维结构网格的索引，便于和其他网格类型统一
    class StructIdxProxy : public IdxProxy
    {
    public:
        StructIdxProxy() = default;
        StructIdxProxy(int ni, int nj, int nk);
        ~StructIdxProxy();
        void SetIdx(int i, int j, int k);
        void SetI(int i);
        void SetJ(int j);
        void SetK(int k);
        int GetNi() const;
        int GetNj() const;
        int GetNk() const;
        void GetIdxStruct(int &i, int &j, int &k) ;
        using IdxProxy::SetIdx;
        using IdxProxy::GetIdx;
        int GetIdx(int i, int j, int k) const
        {
            return i + j * m_ni + k * m_ni * m_nj;
        }
    private:
        GridStruct *m_grid;
        int m_ni, m_nj, m_nk;
        int m_idx_i, m_idx_j, m_idx_k;
    };
}