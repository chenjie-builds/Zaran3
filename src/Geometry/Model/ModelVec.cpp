#include"ModelVec.h"
#include "MathBasic.h"
using namespace zaran;
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

bool ModelVec::InModel(const DVector3D& pt)const
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

DVector3D ModelVec::GetClosestPoint(const DVector3D& pt)const
{
	DVector3D NearPt = modelVec[0]->GetClosestPoint(pt);
	double min_dist = DistanceOfTwoPoints(NearPt.data(), pt.data());
	DVector3D TmpPt;
	if (modelVec.size() > 0)
	{
		for (int iModel = 1; iModel < modelVec.size(); ++iModel)
		{
			TmpPt = modelVec[iModel]->GetClosestPoint(pt);
			if (DistanceOfTwoPoints(TmpPt.data(), pt.data()) < min_dist)
			{
				min_dist = DistanceOfTwoPoints(TmpPt.data(), pt.data());
				NearPt = TmpPt;
			}
		}
	}
	return NearPt;
}

double ModelVec::NearestDistance(const DVector3D& pt)const
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
