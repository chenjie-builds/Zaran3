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
#include<string>
#include<Eigen\Dense>
class Model
{
public:
	virtual bool InModel(const Eigen::Vector3d & pt)const = 0;//判断点pt是否在
	virtual void GenModelPoint(const double delta) = 0;//生成模型离散点数组
	virtual Eigen::Vector3d GetClosestPoint(const Eigen::Vector3d& pt)const = 0;//求出点pt在模型离散点数组中最近的点
	virtual double NearestDistance(const Eigen::Vector3d& pt)const = 0;
	void SetBoxMin(const Eigen::Vector3d& box_min);
	void SetBoxMax(const Eigen::Vector3d& box_max);
	const Eigen::Vector3d& GetBoxMax()const { return box_max_; }
	const Eigen::Vector3d& GetBoxMin()const { return box_min_; }

	//	virtual vector<Point> ModelPoint() { return modelPoint_; };
private:
	Eigen::Vector3d box_min_, box_max_;
};
