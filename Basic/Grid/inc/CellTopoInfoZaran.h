//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	CellTopoInfoZaran.h														||
//*	@brief	Zaran Cell Topology setting												||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include"CellTopoInfo.h"
namespace zaran
{
    enum class CellType
    {
        Fluid=1,//流体单元
        Solid=2,//固体单元
        FluidSolid=3,//流体单元，邻居单元有固体单元
        SolidFluid=4,//固体单元，邻居单元有流体单元
        Unknown=5,//未知单元

    };
    class CellTopoZaran : public CellTopo
    {
    public:
        CellTopoZaran();
        ~CellTopoZaran();
        void SetType(Array<CellType>& cellType);
        Array<CellType>& GetType();
    private:
        Array<CellType> m_cellType;
    };
} // namespace zaran
