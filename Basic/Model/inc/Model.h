//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Model.h																||
//*	@brief	模型基类, 描述几何外形												||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "BasicType.h"
#include "CommonPara.h"
namespace zaran
{
	class Model
	{
	public:
		virtual bool InModel(const double *point_input) const = 0;		  // 判断点pt是否在
		virtual void GetClosestPoint(const double *point_input,double * point_find) const = 0; // 求出点pt在模型离散点数组中最近的点
		virtual double GetClosestDistance(const double *point_input) const = 0;
		const Box &GetBox() const { return m_box; }
		void SetBox(const Box &box) { m_box = box; }
	private:
		// 模型的包围盒
		Box m_box;
	};

}