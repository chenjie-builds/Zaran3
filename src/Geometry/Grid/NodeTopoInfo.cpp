#include "NodeTopoInfo.h"
using namespace zaran;
NodeTopoInfo::NodeTopoInfo()
{

}

NodeTopoInfo::~NodeTopoInfo()
{

}

void NodeTopoInfo::SetType(const NodeType& type)
{
	type_ = type;
}

void NodeTopoInfo::SetTag(const int& tag)
{
	tag_ = tag;
}


const NodeType& NodeTopoInfo::GetType() const
{
	return type_;
}

const int& NodeTopoInfo::GetTag() const
{
	return tag_;
}




std::vector<int>& NodeTopoInfo::GetNeighborCloud()
{
	return neighborCloud_;
}

Eigen::Vector3d& NodeTopoInfo::GetCoordinate()
{
	return coordinate_;
}

void NodeTopoInfo::SetCoordinate(const double& x, const double& y, const double& z/*=0*/)
{
	coordinate_[0] = x;
	coordinate_[1] = y;
	coordinate_[2] = z;
}

void NodeTopoInfo::SetNeighborTemplateI(std::vector<int>& neighborTemplateI)
{
	neighborTemplateI_ = neighborTemplateI;
}

void NodeTopoInfo::SetNeighborTemplateJ(std::vector<int>& neighborTemplateJ)
{
	neighborTemplateJ_ = neighborTemplateJ;
}

void NodeTopoInfo::SetNeighborTemplateK(std::vector<int>& neighborTemplateK)
{
	neighborTemplateK_ = neighborTemplateK;
}

std::vector<int>& NodeTopoInfo::GetNeighborTemplateI()
{
	return neighborTemplateI_;
}

std::vector<int>& NodeTopoInfo::GetNeighborTemplateJ()
{
	return neighborTemplateJ_;
}

std::vector<int>& NodeTopoInfo::GetNeighborTemplateK()
{
	return neighborTemplateK_;
}

void NodeTopoInfo::SetNeighborCloud(const std::vector<int>& neighborCloud)
{
	neighborCloud_ = neighborCloud;
}
