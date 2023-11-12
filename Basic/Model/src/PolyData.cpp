#include"PolyData.h"
#include <vtkFeatureEdges.h>
#include <vtkSTLReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include<vtkVertexGlyphFilter.h>
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
	void PolyDataModel::SetPolyData(vtkSmartPointer<vtkPolyData>& polyData, double tolerance)
	{
		m_polyData = polyData;
		m_tolerance = tolerance;
		m_cell_locator = vtkSmartPointer<vtkOBBTree>::New();
		m_cell_locator->SetDataSet(m_polyData);
		m_cell_locator->BuildLocator();



	}
	bool PolyDataModel::InModel(const DVector3D& pt)const
	{
		int inside = m_cell_locator->InsideOrOutside(pt.data());
		if (inside == -1)
			return true;
		else
			return false;
	}
	DVector3D PolyDataModel::GetClosestPoint(const DVector3D& pt) const
	{
		DVector3D closestPoint;
		vtkIdType cellId;
		int subId;
		double dist2;
		m_cell_locator->FindClosestPoint(pt.data(), closestPoint.data(), cellId, subId, dist2);
		return closestPoint;
	}
	double PolyDataModel::NearestDistance(const DVector3D& pt) const
	{
		DVector3D closestPoint;
		vtkIdType cellId;
		int subId;
		double dist2;
		m_cell_locator->FindClosestPoint(pt.data(), closestPoint.data(), cellId, subId, dist2);
		return sqrt(dist2);
	}
	void PolyDataModel::ShowModel()
	{
		vtkNew<vtkVertexGlyphFilter> glyphFilter;
		glyphFilter->SetInputData(m_polyData);
		glyphFilter->Update();

		// Visualize

		vtkNew<vtkPolyDataMapper> mapper;
		mapper->SetInputConnection(glyphFilter->GetOutputPort());
		vtkNew<vtkNamedColors> colors;
		vtkNew<vtkActor> actor;
		actor->SetMapper(mapper);

		vtkNew<vtkRenderer> renderer;
		renderer->AddActor(actor);
		renderer->SetBackground(colors->GetColor3d("Gainsboro").GetData());

		vtkNew<vtkRenderWindow> renderWindow;
		renderWindow->AddRenderer(renderer);
		renderWindow->SetWindowName("ReadTextFile");

		vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
		renderWindowInteractor->SetRenderWindow(renderWindow);

		renderWindow->Render();
		renderWindowInteractor->Start();

	}
}
