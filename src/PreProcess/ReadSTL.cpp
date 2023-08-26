#include"ReadSTL.h"
#include<vtkFeatureEdges.h>
void STLReader::ReadSTLFile(const char* filename)
{

	vtkSmartPointer<vtkSTLReader> reader =
		vtkSmartPointer<vtkSTLReader>::New();

	reader->SetFileName(filename);

	reader->Update();

	vtkSmartPointer<vtkPolyData> mesh = reader->GetOutput();

	vtkSmartPointer<vtkPoints> points = mesh->GetPoints();

	vtkSmartPointer<vtkDataArray> dataArray = points->GetData();

	vtkCellArray* cells = mesh->GetPolys();
	cells->InsertNextCell({ 0,1,3 });

	vtkIdType numberOfCells = mesh->GetNumberOfCells();
	auto numberOfPoints = mesh->GetNumberOfPoints();
	auto numberOfFaces = mesh->GetNumberOfLines();
	auto numberOfEdges = mesh->GetNumberOfStrips();
	std::cout << "Number of points: " << numberOfPoints << std::endl;
	std::cout << "Number of faces: " << numberOfFaces << std::endl;
	std::cout << "Number of edges: " << numberOfEdges << std::endl;
	vtkSmartPointer<vtkIdList> faceIndex = vtkSmartPointer<vtkIdList>::New();
	vtkIdType vertexIndex = 0;
	double vertexArray[3] = { 0.0, 0.0, 0.0 };
	for (size_t i = 0; i < numberOfCells; i++)
	{

		mesh->GetCellPoints(i, faceIndex);
		faceIndex->GetNumberOfIds();
		vertexIndex = faceIndex->GetId(0);

		vertexArray[0] = dataArray->GetComponent(vertexIndex, 0);

		vertexArray[1] = dataArray->GetComponent(vertexIndex, 1);

		vertexArray[2] = dataArray->GetComponent(vertexIndex, 2);
		//std::cout << vertexArray[0] << " " << vertexArray[1] << " " << vertexArray[2] << std::endl;
	}
	vtkSmartPointer<vtkFeatureEdges> featureEdges =
		vtkSmartPointer<vtkFeatureEdges>::New();
	featureEdges->SetInputData(mesh);
	featureEdges->BoundaryEdgesOn();
	featureEdges->FeatureEdgesOff();
	featureEdges->ManifoldEdgesOff();
	featureEdges->NonManifoldEdgesOff();
	featureEdges->Update();
	int numberOfOpenEdges = featureEdges->GetOutput()->GetNumberOfCells();
	if (numberOfOpenEdges > 0)
	{
		std::cout << "Surface is not closed" << std::endl;
	}
	else
	{
		std::cout << "Surface is closed" << std::endl;
	}

}