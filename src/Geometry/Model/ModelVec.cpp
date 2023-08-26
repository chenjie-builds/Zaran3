#include"ModelVec.h"
#include "MathBasic.h"
void ModelVec::AddMod(std::vector<Model*> modVec)
{
	for (int iModel = 0; iModel < modVec.size(); ++iModel)
		AddMod(modVec[iModel]);
}
void ModelVec::AddMod(ModelVec* mod)
{
	for (int iModel = 0; iModel < mod->modelVec.size(); ++iModel)
		AddMod(mod->modelVec[iModel]);
}

bool ModelVec::InModel(const Eigen::Vector3d& pt)const
{
	bool inMod = false;
	for (int iModel = 0; iModel < modelVec.size(); ++iModel)
	{
		inMod = inMod || modelVec[iModel]->InModel(pt);
	}
	return inMod;
}

void ModelVec::GenPointCloud(const double delta)
{
	for (int iModel = 0; iModel < modelVec.size(); ++iModel)
	{
		modelVec[iModel]->GenModelPoint(delta);
	}
}

Eigen::Vector3d ModelVec::GetClosestPoint(const Eigen::Vector3d& pt)const
{
	Eigen::Vector3d NearPt = modelVec[0]->GetClosestPoint(pt);
	double min_dist = Distance(NearPt, pt);
	Eigen::Vector3d TmpPt;
	if (modelVec.size() > 0)
	{
		for (int iModel = 1; iModel < modelVec.size(); ++iModel)
		{
			TmpPt = modelVec[iModel]->GetClosestPoint(pt);
			if (Distance(TmpPt, pt) < min_dist)
			{
				min_dist = Distance(TmpPt, pt);
				NearPt = TmpPt;
			}
		}
	}
	return NearPt;
}

double ModelVec::NearestDistance(const Eigen::Vector3d& pt)const
{
	double minDistance = modelVec[0]->NearestDistance(pt);
	double tmpDis;
	if (modelVec.size() > 0)
	{
		for (int iModel = 1; iModel < modelVec.size(); ++iModel)
		{
			minDistance = Min(minDistance, modelVec[iModel]->NearestDistance(pt));
		}
	}
	return minDistance;
}
