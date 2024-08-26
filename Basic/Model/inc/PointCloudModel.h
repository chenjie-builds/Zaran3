//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	PointCloudModel.h													||
//*	@brief	点云模型类															||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include"Model.h"
//#include "Utility/include/KDTree.h"
#include <vtkIdList.h>
#include <vtkNew.h>
#include <vtkPointSource.h>
#include <vtkPoints.h>
#include <vtkKdTreePointLocator.h>
/*
点云类，模型是由点云组成
*/
namespace zaran
{
	class PointCloudModel :public Model
	{
	private:
		vtkSmartPointer< vtkKdTreePointLocator> pointTree_;
	public:
		PointCloudModel() {};
		PointCloudModel(const Array<DVector3D>& point_list);
		bool InModel(const DVector3D& pt)const  override;
		DVector3D GetClosestPoint(const DVector3D& pt)const override;
		void GenModelPoint(const double delta) override;//生成模型离散点数组
		double NearestDistance(const DVector3D& pt)const override;//求出pt离模型最近的点

	};
}