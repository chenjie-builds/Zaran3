#include"STL.h"
namespace zaran
{
	bool PolyDataModel::InModel(const DVector3D& pt) const
	{
		//int inside = m_cell_locator->
		//if (inside == 1)
		//	return true;
		//else
			return false;
	}
	bool PolyDataModel::IsClosed()
	{
		return false;
	}
	void PolyDataModel::SetTolerance(double tol)
	{
		m_tolerance = tol;
	}
	//bool PolyDataModel::IsInside(double* p)
	//{
	//	int inside = m_cell_locator->InsideOrOutside(p);
	//	if (inside == 1)
	//		return true;
	//	else
	//		return false;
	//}

	//void PolyDataModel::FindClosestPoint(double* p, double* closestPoint)
	//{
	//	vtkIdType cellId;
	//	int subId;
	//	double dist2;
	//	m_cell_locator->FindClosestPoint(p, closestPoint, cellId, subId, dist2);
	//}
}
