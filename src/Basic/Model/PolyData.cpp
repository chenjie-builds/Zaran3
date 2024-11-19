#include "PolyData.h"
#include <vtkFeatureEdges.h>
#include <vtkSTLReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkVertexGlyphFilter.h>
namespace zaran
{

	bool PolyDataModel::IsClosed()
	{
		vtkNew<vtkFeatureEdges> featureEdges;
		featureEdges->FeatureEdgesOff();
		featureEdges->BoundaryEdgesOn();
		featureEdges->NonManifoldEdgesOn();
		featureEdges->ManifoldEdgesOff();
		featureEdges->SetInputData(m_polyData);
		featureEdges->Update();
		int numberOfOpenEdges = featureEdges->GetOutput()->GetNumberOfCells();
		if (numberOfOpenEdges == 0)
			return true;
		else
			return false;
	}
	void PolyDataModel::SetTolerance(double tol)
	{
		m_tolerance = tol;
	}
	void PolyDataModel::SetPolyData(vtkSmartPointer<vtkPolyData> &polyData, double tolerance)
	{
		m_polyData = polyData;
		m_tolerance = tolerance;
		m_cell_locator = vtkSmartPointer<vtkOBBTree>::New();
		m_cell_locator->SetDataSet(m_polyData);
		m_cell_locator->BuildLocator();
		m_distance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
		m_distance->SetInput(m_polyData);
		m_cellLocator = vtkSmartPointer<vtkCellLocator>::New();
		m_cellLocator->SetDataSet(m_polyData);
		m_cellLocator->BuildLocator();
		double bounds[6];
		m_polyData->GetBounds(bounds);
		Box box;
		box.x_min = bounds[0];
		box.x_max = bounds[1];
		box.y_min = bounds[2];
		box.y_max = bounds[3];
		box.z_min = bounds[4];
		box.z_max = bounds[5];
		SetBox(box);
	}
	bool PolyDataModel::InModel(const double *pt) const
	{
		int inside = m_cell_locator->InsideOrOutside(pt);
		if (inside == -1)
			return true;
		else
			return false;
	}
	void PolyDataModel::GetClosestPoint(const double *point_input, double *point_find) const
	{
		vtkIdType cellId;
		int subId;
		double dist2;
		m_cellLocator->FindClosestPoint(point_input, point_find, cellId, subId, dist2);
	}
	double PolyDataModel::GetClosestDistance(const double *pt) const
	{
		double point_find[3];
		vtkIdType cellId;
		int subId;
		double dist2;
		return abs(m_distance->FunctionValue(pt));
		// m_cellLocator->FindClosestPoint(pt, point_find, cellId, subId, dist2);
		// return dist2;
	}
}
