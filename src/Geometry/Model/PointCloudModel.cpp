#include"PointCloudModel.h"
#include"CommonPara.h"
#include "MathBasic.h"
PointCloudModel::PointCloudModel(const std::vector<Eigen::Vector3d>& point_list)
{
	vtkNew<vtkPoints>points;
	Eigen::Vector3d max, min;
	max = Eigen::Vector3d(-LARGE_NUMBER, -LARGE_NUMBER, -LARGE_NUMBER);
	min = Eigen::Vector3d(LARGE_NUMBER, LARGE_NUMBER, LARGE_NUMBER);
	for (size_t iPoint = 0; iPoint < point_list.size(); ++iPoint)
	{
		auto& x = point_list[iPoint].x();
		auto& y = point_list[iPoint].y();
		auto& z = point_list[iPoint].z();
		points->InsertPoint(iPoint, x, y, z);
		max = { std::max(max.x(), x), std::max(max.y(), y), std::max(max.z(), z) };
		min = { std::min(min.x(), x), std::min(min.y(), y), std::min(min.z(), z) };

	}
	this->SetBoxMax(max);
	this->SetBoxMin(min);
	vtkNew<vtkPolyData> polydata;
	polydata->SetPoints(points);
	pointTree_ = vtkNew<vtkKdTreePointLocator>();
	pointTree_->SetDataSet(polydata);
	pointTree_->BuildLocator();
}

bool PointCloudModel::InModel(const Eigen::Vector3d& pt)const
{
	//点云模型，所有输入的点都在外部
	return false;
}

void PointCloudModel::GenModelPoint(const double delta)
{
	// Model is construct by point cloud
	// do nothing!
}
Eigen::Vector3d PointCloudModel::GetClosestPoint(const Eigen::Vector3d& pt)const
{
	double coord[3] = { pt.x(), pt.y(), pt.z() };
	vtkIdType id = pointTree_->FindClosestPoint(coord);
	auto closetPt = pointTree_->GetDataSet()->GetPoint(id);
	return Eigen::Vector3d(closetPt[0], closetPt[1], closetPt[2]);
}
double PointCloudModel::NearestDistance(const Eigen::Vector3d& pt)const
{
	Eigen::Vector3d near_point = GetClosestPoint(pt);
	return Distance(pt, near_point);
}
