/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file Model.h
 * \brief Model Base class, Define the Geometry Model.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
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