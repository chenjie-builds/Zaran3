#include"PolygonModel.h"
#include"CommonPara.h"
#include "MathBasic.h"
#include <Eigen/Dense>
PloygonModel::PloygonModel(std::vector<Eigen::Vector3d>& vertex_)
{
	vertex_ = vertex_;
	vertexIndex_.resize(vertex_.size());
}
bool PloygonModel::InModel(const Eigen::Vector3d& pt)const
{
	int iVertex, jVertex;
	bool isInModel = false;
	for (iVertex = 0, jVertex = vertex_.size() - 1; iVertex < vertex_.size(); jVertex = iVertex++)
	{
		if (((vertex_[iVertex].y() > pt.y()) != (vertex_[jVertex].y() > pt.y())) &&
			(pt.x() < (vertex_[jVertex].x() - vertex_[iVertex].x()) * (pt.y() - vertex_[iVertex].y()) / (vertex_[jVertex].y() - vertex_[iVertex].y()) + vertex_[iVertex].x()))
			isInModel = !isInModel;
	}
	return isInModel;
}
void PloygonModel::GenModelPoint(const double delta)
{
	int iVertex;
	double dis, dx, dy;
	modelPoint_.push_back(vertex_[0]);
	vertexIndex_[0] = 0;
	for (iVertex = 0; iVertex < vertex_.size() - 1; ++iVertex)
	{
		dis = Distance(vertex_[iVertex].data(), vertex_[iVertex + 1].data());
		dx = (vertex_[iVertex + 1].x() - vertex_[iVertex].x()) * delta / dis;
		dy = (vertex_[iVertex + 1].y() - vertex_[iVertex].y()) * delta / dis;
		while (Distance(modelPoint_[modelPoint_.size() - 1].data(), vertex_[iVertex + 1].data()) > delta)
			modelPoint_.push_back({ modelPoint_[modelPoint_.size() - 1].x() + dx,modelPoint_[modelPoint_.size() - 1].y() + dy,0 });
		modelPoint_.push_back(vertex_[iVertex + 1]);
		vertexIndex_[iVertex + 1] = modelPoint_.size() - 1;
	}
	dis = Distance(vertex_[vertex_.size() - 1].data(), vertex_[0].data());
	dx = (vertex_[0].x() - vertex_[vertex_.size() - 1].x()) * delta / dis;
	dy = (vertex_[0].y() - vertex_[vertex_.size() - 1].y()) * delta / dis;
	while (Distance(modelPoint_[modelPoint_.size() - 1].data(), vertex_[0].data()) > delta)
		modelPoint_.push_back({ modelPoint_[modelPoint_.size() - 1].x() + dx,modelPoint_[modelPoint_.size() - 1].y() + dy ,0});
}
Eigen::Vector3d PloygonModel::GetClosestPoint(const Eigen::Vector3d& pt)const
{
	int iPoint, id = -1;
	double mindis = DBL_MAX;
	double dis;
	for (iPoint = 0; iPoint < modelPoint_.size(); ++iPoint)
	{
		dis = Distance(modelPoint_[iPoint], pt);
		mindis = std::min(mindis, dis);
		if (abs(mindis - dis) < SMALL_NUMBER)
			id = iPoint;
	}
	return modelPoint_[id];
}
double PloygonModel::NearestDistance(const Eigen::Vector3d& pt)const
{
	int iVertex, jVertex;
	double minDistance = DBL_MAX;
	for (iVertex = 0, jVertex = vertex_.size() - 1; iVertex < vertex_.size(); jVertex = iVertex++)
	{
		double x1 = vertex_[iVertex].x();
		double y1 = vertex_[iVertex].y();
		double x2 = vertex_[jVertex].x();
		double y2 = vertex_[jVertex].y();
		//直线方程Ax+By+C=0
		double A = y2 - y1;
		double B = x1 - x2;
		double C = x2 * y1 - x1 * y2;

		Eigen::Vector3d ptCross((B * B * pt.x() - A * B * pt.y() - A * C) / (A * A + B * B),
			(A * A * pt.y() - A * B * pt.x() - B * C) / (A * A + B * B),0);
		double d1 = Distance(vertex_[iVertex], ptCross);
		double d2 = Distance(vertex_[jVertex], ptCross);
		double d3 = Distance(vertex_[iVertex], vertex_[jVertex]);
		double d4 = Distance(ptCross, pt);
		double d5 = d3 - d2 - d1;
		double d6 = abs(A * pt.x() + B * pt.y() + C) / sqrt(A * A + B * B);
		double d7 = A * x1 + B * y1 + C;
		double d8 = A * x2 + B * y2 + C;
		if (abs(Distance(vertex_[iVertex], ptCross) + Distance(vertex_[jVertex], ptCross) - Distance(vertex_[iVertex], vertex_[jVertex])) > 1e-15)
			Min(minDistance, Min(Distance(vertex_[iVertex], pt), Distance(vertex_[jVertex], pt)));
		else
			Min(minDistance, Distance(ptCross, pt));

	}
	return minDistance;
}