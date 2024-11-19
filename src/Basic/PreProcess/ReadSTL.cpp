#include "ReadSTL.h"
#include <vtkSTLReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkNamedColors.h>
#include <vtkProperty.h>
#include "Log.h"
namespace zaran
{

	void STLReader::ReadSTLFile(const char *filename)
	{
		Log::info("Read STL file: {}", filename);
		vtkSmartPointer<vtkSTLReader> reader =
			vtkSmartPointer<vtkSTLReader>::New();
		reader->SetFileName(filename);
		reader->Update();
		m_poly_data = vtkSmartPointer<vtkPolyData>::New();
		m_poly_data = reader->GetOutput();
		vtkIdType numberOfCells = m_poly_data->GetNumberOfCells();
		auto numberOfPoints = m_poly_data->GetNumberOfPoints();
		auto numberOfFaces = m_poly_data->GetNumberOfLines();
		auto numberOfEdges = m_poly_data->GetNumberOfStrips();
		Log::info("Read STL file done");
		Log::info("Number of cells: {}", numberOfCells);
		Log::info("Number of points: {}", numberOfPoints);
	}
}