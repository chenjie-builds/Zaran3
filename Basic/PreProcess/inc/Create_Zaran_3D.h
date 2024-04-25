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
#include "Grid_Zaran_3D.h"
#include"ReadSTL.h"
#include"PolyData.h"
#include <vtkImplicitPolyDataDistance.h>
#include<vtkSelectEnclosedPoints.h>
namespace zaran
{
    class GridFactoryZaran3D :public GridListFactory
    {
    public:
        GridFactoryZaran3D() {};
        void Create(Grid*& grid) override;
    private:
        void CreateStructPart(Grid_Zaran_3D* grid);
        void TagCell(Grid_Zaran_3D* grid);
        void CrateBoundPatch(Grid_Zaran_3D* grid);
        void ReadModel();
    private:
        Ptr<PolyDataModel> m_polyDataModel;
    };
}