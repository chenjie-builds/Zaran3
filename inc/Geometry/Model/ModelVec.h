//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	ModelVec.h															||
//*	@brief	模型数组类, 存储一系列模型												||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include"PointCloudModel.h"
#include <vector>
/*
外形几何类，所有的外形都存在这里
modNum：		外形的个数
geo：		存储外形的指针
filename：	存储外形信息文件的文件名
目前只有一个多边形外形2021年1月20日11:30:47
*/
class ModelVec
{
private:
	std::vector<Model*> modelVec;
public:
	ModelVec() { modelVec.resize(0); };
	ModelVec(std::vector<Model*> mod) :modelVec(mod) {};
	ModelVec(const ModelVec& mod) { modelVec = mod.modelVec; };
	void AddMod(Model* mod) { modelVec.push_back(mod); }
	void AddMod(std::vector<Model*> mod);
	void AddMod(ModelVec* mod);
	const Model& GetModel(size_t iModel)const { return *modelVec[iModel]; }
	Eigen::Vector3d GetBoxMax(size_t imodel) { return modelVec[imodel]->GetBoxMax(); }
	Eigen::Vector3d GetBoxMin(size_t imodel) { return modelVec[imodel]->GetBoxMin(); }
	bool InModel(const Eigen::Vector3d& pt)const;
	void GenPointCloud(const double delta);
	Eigen::Vector3d GetClosestPoint(const Eigen::Vector3d& pt)const;
	double NearestDistance(const Eigen::Vector3d& pt)const;
	const size_t GetModelNum()const { return modelVec.size(); }
};