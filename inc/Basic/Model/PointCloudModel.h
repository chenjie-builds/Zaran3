/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file PointCloudModel.h
 * \brief Point cloud model class.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
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
		PointCloudModel(const dynamic_array<dynamic_array<double>>& point_list);
		bool InModel(const double* point_input)const  override;
		void GetClosestPoint(const double*point_input,double*point_find)const override;
		double GetClosestDistance(const double* point_input)const override;//求出pt离模型最近的点

	};
}