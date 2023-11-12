//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Create_Zaran_3D.h													||
//*	@brief	生成扎染三维网格           											 ||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "GridListFactory.h"
#include"ReadSTL.h"
#include"PolyData.h"
#include <vtkImplicitPolyDataDistance.h>
#include<vtkSelectEnclosedPoints.h>
namespace zaran
{
    class GridListFactoryZaran3D :public GridListFactory
    {
    public:
        GridListFactoryZaran3D() {};
        void Create(Ptr<GridList>& gridList) override;
    private:
        void CreateStructPart(Ptr<GridList>& gridList);
        void TagCell(Ptr<GridList>& gridList);
        void CrateBoundPatch(Ptr<GridList>& gridList);
        void ReadModel();
    private:
        Ptr<PolyDataModel> m_polyDataModel;
    };
}