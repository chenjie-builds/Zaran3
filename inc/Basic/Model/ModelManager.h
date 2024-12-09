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
#include "PointCloudModel.h"
#include "BasicType.h"
/*
外形几何类，所有的外形都存在这里
modNum：		外形的个数
geo：		存储外形的指针
filename：	存储外形信息文件的文件名
目前只有一个多边形外形2021年1月20日11:30:47
*/
namespace zaran
{
	class ModelManager
	{
	public:
		ModelManager();
		ModelManager(std::vector<std::shared_ptr<Model >> model) : m_model(model) {};
		ModelManager(const ModelManager &model) { m_model = model.m_model; };
		void AddModel(std::shared_ptr<Model> model);
		void AddModel(std::vector<std::shared_ptr<Model>> model);
		void AddModel(ModelManager *model);
		const std::shared_ptr<Model> GetModel(size_t iModel) const { return m_model[iModel]; }
		bool InModel(const double *pt) const;
		void GetClosestPoint(const double *pt, double *point_find) const;
		double GetClosestDistance(const double *pt) const;
		const size_t GetModelNum() const { return m_model.size(); }
		Box& GetBox() { return m_box; }

	private:
		std::vector<std::shared_ptr<Model>> m_model;
		Box m_box;
	};
}