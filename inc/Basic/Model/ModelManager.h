/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file ModelManager.h
 * \brief Model Manager class, used to manage multiple models.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
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
		ModelManager(dynamic_array<shared_ptr<Model >> model) : m_model(model) {};
		ModelManager(const ModelManager &model) { m_model = model.m_model; };
		void AddModel(shared_ptr<Model> model);
		void AddModel(dynamic_array<shared_ptr<Model>> model);
		void AddModel(ModelManager *model);
		const shared_ptr<Model> GetModel(size_t iModel) const { return m_model[iModel]; }
		bool InModel(const double *pt) const;
		void GetClosestPoint(const double *pt, double *point_find) const;
		double GetClosestDistance(const double *pt) const;
		const size_t GetModelNum() const { return m_model.size(); }
		Box& GetBox() { return m_box; }

	private:
		dynamic_array<shared_ptr<Model>> m_model;
		Box m_box;
	};
}