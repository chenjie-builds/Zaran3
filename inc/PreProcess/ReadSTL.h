#pragma once
// read stl file by using vtk
#include <vtkSTLReader.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
class STLReader
{
public:
	STLReader() {};
	~STLReader() {};
	void ReadSTLFile(const char* filename);
private:
	vtkSmartPointer<vtkSTLReader> reader;
};