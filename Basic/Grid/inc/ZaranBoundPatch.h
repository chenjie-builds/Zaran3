//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	ZaranBoundPatch.h														||
//*	@brief	Zaran 边界补丁														||
//*			过渡节点编号，边界坐标位置                                  			||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "BasicType.h"
namespace zaran
{
    class ZaranBoundPatch
    {
    public:
        ZaranBoundPatch() :m_index(), m_coordinate(), m_normal() {}
        ~ZaranBoundPatch() {}
        void SetIndex(const Array<IArray>& index) { m_index = index; }
        void SetCoordinate(const Array<DVector3D>& coordinate) { m_coordinate = coordinate; }
        void SetNormal(const Array<DVector3D>& normal) { m_normal = normal; }
        Array<IArray>& GetIndex() { return m_index; }
        Array<DVector3D>& GetCoordinate() { return m_coordinate; }
         Array<DVector3D>& GetNormal()  { return m_normal; }
         int GetPatchNum() { return m_index.size(); }
    private:
        // 过渡节点编号
        Array<IArray> m_index;
        // 边界坐标位置
        Array<DVector3D> m_coordinate;
        // 边界法向量
        Array<DVector3D> m_normal;
    };

} // namespace zaran
