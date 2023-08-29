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
#include"BasicType.h"
namespace zaran
{
	class Model
	{
	public:
		virtual bool InModel(const DVector3D& pt)const = 0;//判断点pt是否在
		virtual void GenModelPoint(const double delta) = 0;//生成模型离散点数组
		virtual DVector3D GetClosestPoint(const DVector3D& pt)const = 0;//求出点pt在模型离散点数组中最近的点
		virtual double NearestDistance(const DVector3D& pt)const = 0;
		void SetBoxMin(const DVector3D& box_min);
		void SetBoxMax(const DVector3D& box_max);
		const DVector3D& GetBoxMax()const { return box_max_; }
		const DVector3D& GetBoxMin()const { return box_min_; }

		//	virtual vector<Point> ModelPoint() { return modelPoint_; };
	private:
		DVector3D box_min_, box_max_;
	};

}