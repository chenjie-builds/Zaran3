#include "ModelManager.h"
#include "MathBasic.h"
using namespace zaran;
zaran::ModelManager::ModelManager()
{
	m_box.x_min = LARGE_NUMBER;
	m_box.x_max = -LARGE_NUMBER;
	m_box.y_min = LARGE_NUMBER;
	m_box.y_max = -LARGE_NUMBER;
	m_box.z_min = LARGE_NUMBER;
	m_box.z_max = -LARGE_NUMBER;
	m_model.clear();
}
void zaran::ModelManager::AddModel(Model *model)
{
	m_model.push_back(model);
	auto box = model->GetBox();
	m_box.x_min = Min(m_box.x_min, box.x_min);
	m_box.x_max = Max(m_box.x_max, box.x_max);
	m_box.y_min = Min(m_box.y_min, box.y_min);
	m_box.y_max = Max(m_box.y_max, box.y_max);
	m_box.z_min = Min(m_box.z_min, box.z_min);
	m_box.z_max = Max(m_box.z_max, box.z_max);
}
void ModelManager::AddModel(std::vector<Model *> modVec)
{
	for (int iModel = 0; iModel < modVec.size(); ++iModel)
		AddModel(modVec[iModel]);
}
void ModelManager::AddModel(ModelManager *mod)
{
	for (int iModel = 0; iModel < mod->m_model.size(); ++iModel)
		AddModel(mod->m_model[iModel]);
}

bool ModelManager::InModel(const double *pt) const
{
	bool inMod = false;
	for (int iModel = 0; iModel < m_model.size(); ++iModel)
	{
		inMod = inMod || m_model[iModel]->InModel(pt);
	}
	return inMod;
}

void zaran::ModelManager::GetClosestPoint(const double *pt, double *point_find) const
{
	if (m_model.size() == 0)
	{
		point_find[0] = 0;
		point_find[1] = 0;
		point_find[2] = 0;
		return;
	}
	double min_dist = LARGE_NUMBER;
	double dist;
	for (int iModel = 0; iModel < m_model.size(); ++iModel)
	{
		dist = m_model[iModel]->GetClosestDistance(pt);
		if (dist < min_dist)
		{
			min_dist = dist;
			m_model[iModel]->GetClosestPoint(pt, point_find);
		}
	}
}

double ModelManager::GetClosestDistance(const double *pt) const
{
	double min_dis = LARGE_NUMBER;
	for (int iModel = 0; iModel < m_model.size(); ++iModel)
	{
		min_dis = Min(min_dis, m_model[iModel]->GetClosestDistance(pt));
	}
	return min_dis;
}
