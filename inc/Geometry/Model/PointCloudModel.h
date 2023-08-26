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
#include <vector>
//#include "Utility/include/KDTree.h"
#include <vtkIdList.h>
#include <vtkNew.h>
#include <vtkPointSource.h>
#include <vtkPoints.h>
#include <vtkKdTreePointLocator.h>
/*
点云类，模型是由点云组成
*/
class PointCloudModel :public Model
{
private:
	vtkSmartPointer< vtkKdTreePointLocator> pointTree_;
public:
	PointCloudModel() {};
	PointCloudModel(const std::vector<Eigen::Vector3d>& point_list);
	bool InModel(const Eigen::Vector3d& pt)const  override;
	Eigen::Vector3d GetClosestPoint(const Eigen::Vector3d& pt)const override;
	void GenModelPoint(const double delta) override;//生成模型离散点数组
	double NearestDistance(const Eigen::Vector3d& pt)const override;//求出pt离模型最近的点

};