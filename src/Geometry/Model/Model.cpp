#include"Model.h"
#include"MathBasic.h"
#include<iostream>
#include<fstream>
#include<algorithm>
void Model::SetBoxMin(const Eigen::Vector3d& box_min)
{
	box_min_ = box_min;
}
void Model::SetBoxMax(const Eigen::Vector3d& box_max)
{
	box_max_ = box_max;
}
//Point ModBase::GetClosestPoint(const Point& pt)
//{
	//int i, id = -1;
	//double mindis = DBL_MAX;
	//double dis;
	//for (i = 0; i < modelPoint_.size(); i++)
	//{
	//	dis = distance(modelPoint_[i], pt);
	//	mindis = std::min(mindis, dis);
	//	if (abs(mindis - dis) < delta)
	//		id = i;
	//}
	//return modelPoint_[id];
//}

//double ModBase::NearestDistance(const Point& pt)
//{
	//int i, j;
	//double min_dist = DBL_MAX;
	//for (i = 0; i < modelPoint_.size(); j = i++)
	//{
	//	min_dist = std::min(min_dist, distance(modelPoint_[i], pt));
	//}
	//return min_dist;
//}