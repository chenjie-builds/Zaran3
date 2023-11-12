#pragma once
// read stl file by using vtk
#include <vtkSTLReader.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include"BasicType.h"
namespace zaran
{

	class STLReader
	{
	public:
		STLReader() {};
		~STLReader() {};
		void ReadSTLFile(const char* filename);
		vtkSmartPointer<vtkPolyData>&  GetMesh() { return m_mesh; }
	private:
		vtkSmartPointer<vtkPolyData> m_mesh;
	};
}