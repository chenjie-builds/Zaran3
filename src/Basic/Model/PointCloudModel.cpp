#include "PointCloudModel.h"
#include "CommonPara.h"
#include "MathBasic.h"
using namespace zaran;
PointCloudModel::PointCloudModel(const dynamic_array<Eigen::Vector3d> &point_list)
{
	vtkNew<vtkPoints> points;
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
PointCloudModel::PointCloudModel(const dynamic_array<dynamic_array<double>>& point_list)
{
	vtkNew<vtkPoints> points;
    Box box;
	for (size_t iPoint = 0; iPoint < point_list.size(); ++iPoint)
	{
		auto& x = point_list[iPoint][0];
		auto& y = point_list[iPoint][1];
		auto& z = point_list[iPoint][2];
		points->InsertPoint(iPoint, x, y, z);
		box.x_min = std::min(box.x_min, x);
		box.x_max = std::max(box.x_max, x);
		box.y_min = std::min(box.y_min, y);
		box.y_max = std::max(box.y_max, y);
		box.z_min = std::min(box.z_min, z);
		box.z_max = std::max(box.z_max, z);
	}
	SetBox(box);
	vtkNew<vtkPolyData> point_cloud_poly_data;
	point_cloud_poly_data->SetPoints(points);
	m_point_cloud = vtkNew<vtkKdTreePointLocator>();
	m_point_cloud->SetDataSet(point_cloud_poly_data);
	m_point_cloud->BuildLocator();

}

bool PointCloudModel::InModel(const double *point_input) const
{
	// 点云模型，所有输入的点都在外部
	return false;
}

void PointCloudModel::GetClosestPoint(const double *point_input, double *point_find) const
{
	double input_coord[3] = {point_input[0], point_input[1], point_input[2]};
	vtkIdType closest_point_id = m_point_cloud->FindClosestPoint(input_coord);
	auto closest_point = m_point_cloud->GetDataSet()->GetPoint(closest_point_id);
	point_find[0] = closest_point[0];
	point_find[1] = closest_point[1];
	point_find[2] = closest_point[2];
}
double PointCloudModel::GetClosestDistance(const double *point_input) const
{
	double input_coord[3] = {point_input[0], point_input[1], point_input[2]};
	vtkIdType closest_point_id = m_point_cloud->FindClosestPoint(input_coord);
	auto closest_point = m_point_cloud->GetDataSet()->GetPoint(closest_point_id);
	return DistanceOfTwoPoints(point_input, closest_point);
}
