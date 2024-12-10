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
		vtkSmartPointer< vtkKdTreePointLocator> m_point_cloud;
	public:
		PointCloudModel() {};
		PointCloudModel(const dynamic_array<Eigen::Vector3d>& point_list);
		bool InModel(const double* point_input)const  override;
		void GetClosestPoint(const double*point_input,double*point_find)const override;
		double GetClosestDistance(const double* point_input)const override;//求出pt离模型最近的点

	};
}