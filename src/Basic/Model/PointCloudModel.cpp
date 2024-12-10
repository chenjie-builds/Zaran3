#include "PointCloudModel.h"
#include "CommonPara.h"
#include "MathBasic.h"
using namespace zaran;
PointCloudModel::PointCloudModel(const dynamic_array<Eigen::Vector3d> &point_list)
{
	vtkNew<vtkPoints> points;
	Eigen::Vector3d max, min;
	max = Eigen::Vector3d(-LARGE_NUMBER, -LARGE_NUMBER, -LARGE_NUMBER);
	min = Eigen::Vector3d(LARGE_NUMBER, LARGE_NUMBER, LARGE_NUMBER);
	Box box;
	for (size_t iPoint = 0; iPoint < point_list.size(); ++iPoint)
	{
		auto &x = point_list[iPoint].x();
		auto &y = point_list[iPoint].y();
		auto &z = point_list[iPoint].z();
		points->InsertPoint(iPoint, x, y, z);
		box.x_min = std::min(box.x_min, x);
		box.x_max = std::max(box.x_max, x);
		box.y_min = std::min(box.y_min, y);
		box.y_max = std::max(box.y_max, y);
		box.z_min = std::min(box.z_min, z);
		box.z_max = std::max(box.z_max, z);
	}
	SetBox(box);
	vtkNew<vtkPolyData> polydata;
	polydata->SetPoints(points);
	m_point_cloud = vtkNew<vtkKdTreePointLocator>();
	m_point_cloud->SetDataSet(polydata);
	m_point_cloud->BuildLocator();
}

bool PointCloudModel::InModel(const double *point_input) const
{
	// 点云模型，所有输入的点都在外部
	return false;
}

void PointCloudModel::GetClosestPoint(const double *point_input, double *point_find) const
{
	double coord[3] = {point_input[0], point_input[1], point_input[2]};
	vtkIdType id = m_point_cloud->FindClosestPoint(coord);
	auto closetPt = m_point_cloud->GetDataSet()->GetPoint(id);
	point_find[0] = closetPt[0];
	point_find[1] = closetPt[1];
	point_find[2] = closetPt[2];
}
double PointCloudModel::GetClosestDistance(const double *point_input) const
{
	double coord[3] = {point_input[0], point_input[1], point_input[2]};
	vtkIdType id = m_point_cloud->FindClosestPoint(coord);
	auto closetPt = m_point_cloud->GetDataSet()->GetPoint(id);
	return DistanceOfTwoPoints(point_input, closetPt);
}
