#include "NodeTopoInfo.h"
using namespace zaran;
NodeTopoInfo::NodeTopoInfo()
{
	type_ = NodeType::undefined;
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




IArray& NodeTopoInfo::GetNeighborCloud()
{
	return neighborCloud_;
}

DVector3D& NodeTopoInfo::GetCoordinate()
{
	return coordinate_;
}

void NodeTopoInfo::SetCoordinate(const double& x, const double& y, const double& z/*=0*/)
{
	coordinate_[0] = x;
	coordinate_[1] = y;
	coordinate_[2] = z;
}

void NodeTopoInfo::SetNeighborTemplateI(IArray& neighborTemplateI)
{
	neighborTemplateI_ = neighborTemplateI;
}

void NodeTopoInfo::SetNeighborTemplateJ(IArray& neighborTemplateJ)
{
	neighborTemplateJ_ = neighborTemplateJ;
}

void NodeTopoInfo::SetNeighborTemplateK(IArray& neighborTemplateK)
{
	neighborTemplateK_ = neighborTemplateK;
}

IArray& NodeTopoInfo::GetNeighborTemplateI()
{
	return neighborTemplateI_;
}

IArray& NodeTopoInfo::GetNeighborTemplateJ()
{
	return neighborTemplateJ_;
}

IArray& NodeTopoInfo::GetNeighborTemplateK()
{
	return neighborTemplateK_;
}

void NodeTopoInfo::SetNeighborCloud(const IArray& neighborCloud)
{
	neighborCloud_ = neighborCloud;
}
