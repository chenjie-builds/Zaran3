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
#include "GridStruct.h"
namespace zaran
{
    // 三维结构网格索引代理类
    // 用于计算三维结构网格的索引，便于和其他网格类型统一
    class StructIdxProxy
    {
    public:
        StructIdxProxy(GridStruct* grid);
        ~StructIdxProxy();
        int GetNi()const;
        int GetNj()const;
        int GetNk()const;
        int GetIdx(int i, int j, int k)const;
        void GetIdxStruct(int idx, int& i, int& j, int& k)const;
    private:
        GridStruct* m_grid;
    };
}