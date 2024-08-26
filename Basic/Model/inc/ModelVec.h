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
#include"BasicType.h"
/*
外形几何类，所有的外形都存在这里
modNum：		外形的个数
geo：		存储外形的指针
filename：	存储外形信息文件的文件名
目前只有一个多边形外形2021年1月20日11:30:47
*/
namespace zaran
{
	class ModelVec
	{
	private:
		Array<Model*> modelVec;
	public:
		ModelVec() { modelVec.resize(0); };
		ModelVec(Array<Model*> mod) :modelVec(mod) {};
		ModelVec(const ModelVec& mod) { modelVec = mod.modelVec; };
		void AddMod(Model* mod) { modelVec.push_back(mod); }
		void AddMod(Array<Model*> mod);
		void AddMod(ModelVec* mod);
		const Model& GetModel(size_t iModel)const { return *modelVec[iModel]; }
		DVector3D GetBoxMax(size_t imodel) { return modelVec[imodel]->GetBoxMax(); }
		DVector3D GetBoxMin(size_t imodel) { return modelVec[imodel]->GetBoxMin(); }
		bool InModel(const DVector3D& pt)const;
		void GenPointCloud(const double delta);
		DVector3D GetClosestPoint(const DVector3D& pt)const;
		double NearestDistance(const DVector3D& pt)const;
		const size_t GetModelNum()const { return modelVec.size(); }
	};
}