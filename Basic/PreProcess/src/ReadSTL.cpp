#include"ReadSTL.h"
#include <vtkSTLReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include<vtkVertexGlyphFilter.h>
#include<vtkNamedColors.h>
#include<vtkProperty.h>
namespace zaran
{

	void STLReader::ReadSTLFile(const char* filename)
	{
		vtkSmartPointer<vtkSTLReader> reader =
			vtkSmartPointer<vtkSTLReader>::New();
		reader->SetFileName(filename);
		reader->Update();
		m_mesh = vtkSmartPointer<vtkPolyData>::New();
		m_mesh = reader->GetOutput();
		vtkIdType numberOfCells = m_mesh->GetNumberOfCells();
		auto numberOfPoints = m_mesh->GetNumberOfPoints();
		auto numberOfFaces = m_mesh->GetNumberOfLines();
		auto numberOfEdges = m_mesh->GetNumberOfStrips();
		std::cout << "Number of points: " << numberOfPoints << std::endl;
		std::cout << "Number of cells: " << numberOfCells << std::endl;
	}
}