#include "Create_Zaran_3D.h"
#include "Grid_Zaran_3D.h"
#include"ReadSTL.h"
#include"PolyData.h"
#include "Log.h"
#include"CellTopoInfoZaran.h"
#include <vtkImplicitPolyDataDistance.h>
#include<vtkCellLocator.h>
#include "MathBasic.h"
#include<fstream>
void zaran::GridListFactoryZaran3D::Create(Ptr<GridList>& gridList)
{
	if (gridList.get() == nullptr)
	{
		gridList = std::make_shared<GridList>();
	}
	ZaranLog::info("Start create grid");
	CreateStructPart(gridList);
	ZaranLog::info("Create grid finished");
	ZaranLog::info("Start read model");
	ReadModel();
	ZaranLog::info("Start tag cell type");
	TagCell(gridList);
	ZaranLog::info("Tag cell type finished");
	CrateBoundPatch(gridList);
}

void zaran::GridListFactoryZaran3D::CreateStructPart(Ptr<GridList>& gridList)
{
	if (gridList.get() == nullptr)
	{
		gridList = std::make_shared<GridList>();
	}

	Ptr<Grid_Zaran_3D>grid = std::make_shared<Grid_Zaran_3D>();
	grid->SetDimension(Dimension::three);
	grid->SetIndex(0);
	grid->SetLevel(0);
	grid->SetName("noname");
	grid->SetType(GridType::Zaran_3D);
	auto& nodeTopo = grid->GetNodeTopo();
	auto& nodeCoord = nodeTopo->GetCoordinate();
	int xNodeNum = GlobalData::GetInt("nx");
	int yNodeNum = GlobalData::GetInt("ny");
	int zNodeNum = GlobalData::GetInt("nz");
	xNodeNum++;
	yNodeNum++;
	zNodeNum++;
	int nodeNum = 0;
	double xMin = GlobalData::GetDouble("xMin");
	double yMin = GlobalData::GetDouble("yMin");
	double zMin = GlobalData::GetDouble("zMin");
	double xMax = GlobalData::GetDouble("xMax");
	double yMax = GlobalData::GetDouble("yMax");
	double zMax = GlobalData::GetDouble("zMax");
	grid->SetBox(xMin, xMax, yMin, yMax, zMin, zMax);
	double dx = (xMax - xMin) / (xNodeNum - 1);
	double dy = (yMax - yMin) / (yNodeNum - 1);
	double dz = (zMax - zMin) / (zNodeNum - 1);
	int i, j, k;
	double x, y, z;
	xNodeNum += 2;
	yNodeNum += 2;
	zNodeNum += 2;
	grid->SetNodeNum((xNodeNum), (yNodeNum), (zNodeNum));
	nodeCoord.resize((xNodeNum) * (yNodeNum) * (zNodeNum));

	auto& nodeType = nodeTopo->GetType();
	nodeType.resize((xNodeNum) * (yNodeNum) * (zNodeNum));
	int iNode = 0;
	for (k = 0; k < zNodeNum; ++k)
	{
		z = k * dz + zMin - dz;
		for (j = 0; j < yNodeNum; ++j)
		{
			y = j * dy + yMin - dy;
			for (i = 0; i < xNodeNum; ++i)
			{
				x = i * dx + xMin - dx;
				iNode = grid->GetNodeIndex(i, j, k);
				nodeCoord[iNode] = { x,y,z };
				if (i == 0 || j == 0 || k == 0 || i == xNodeNum + 1 || j == yNodeNum + 1 || k == zNodeNum + 1)
				{
					nodeType[iNode] = NodeType::ghost;
					continue;
				}
				if (i == 1 || j == 1 || k == 1 || i == xNodeNum || j == yNodeNum || k == zNodeNum)
				{
					nodeType[iNode] = NodeType::inlet;
				}
			}
		}
	}
	ZaranLog::info("Start create cell topo1");
	grid->SetInnerNodeNum((xNodeNum - 4) * (yNodeNum - 4) * (zNodeNum - 4));
	grid->SetTotalNodeNum((xNodeNum) * (yNodeNum) * (zNodeNum));
	grid->SetBoundNodeNum(xNodeNum * yNodeNum * zNodeNum - grid->GetInnerNodeNum());
	grid->SetInterNodeInfo(std::make_shared<InterNodeInfo>());
	grid->SetBoundaryMap(std::make_shared<BoundaryMap>());
	ZaranLog::info("Start create cell topo2");
	Ptr<CellTopoInfoZaran>& cellTopo = grid->GetCellTopo();
	ZaranLog::info("Start create cell topo2");
	auto& cell2node = cellTopo->GetNodeIndex();
	auto& cell_center = cellTopo->GetCenterCoord();
	ZaranLog::info("Start create cell topo");
	ZaranLog::info("cell2node size: {}", cell2node.size());
	ZaranLog::info("cell_center size: {}", cell_center.size());
	ZaranLog::info("xNodeNum: {}, yNodeNum: {}, zNodeNum: {}", xNodeNum, yNodeNum, zNodeNum);
	cell2node.resize((xNodeNum - 1) * (yNodeNum - 1) * (zNodeNum - 1));
	cell_center.resize((xNodeNum - 1) * (yNodeNum - 1) * (zNodeNum - 1));
	int iCell;
	ZaranLog::info("Start create cell topo");
	for (k = 0; k < zNodeNum - 1; ++k)
	{
		for (j = 0; j < yNodeNum - 1; ++j)
		{
			for (i = 0; i < xNodeNum - 1; ++i)
			{
				iCell = grid->GetCellIndex(i, j, k);
				cell2node[iCell] = (IArray{ grid->GetNodeIndex(i, j, k), grid->GetNodeIndex(i + 1, j, k), grid->GetNodeIndex(i + 1, j + 1, k), grid->GetNodeIndex(i, j + 1, k), grid->GetNodeIndex(i, j, k + 1), grid->GetNodeIndex(i + 1, j, k + 1), grid->GetNodeIndex(i + 1, j + 1, k + 1), grid->GetNodeIndex(i, j + 1, k + 1) });
				cell_center[iCell][0] = (i - 0.5) * dx + xMin;
				cell_center[iCell][1] = (j - 0.5) * dy + yMin;
				cell_center[iCell][2] = (k - 0.5) * dz + zMin;
			}
		}
	}
	ZaranLog::info("Create cell topo finished");
	auto boundMap = grid->GetBoundaryMap();
	int cellIndex, innerNodeIndex, ghostNodeIndex;
	DVector3D boundNorm;
	Boundary bound;
	//i方向两个面分别为入口和出口
	for (k = 1; k < zNodeNum - 1; ++k)
	{

		for (j = 1; j < yNodeNum - 1; ++j)
		{
			i = 1;
			cellIndex = grid->GetCellIndex(i, j, k);
			bound.SetNodeIndex(cellIndex);
			innerNodeIndex = grid->GetCellIndex(i + 1, j, k);
			bound.SetInnerIndex(innerNodeIndex);
			ghostNodeIndex = grid->GetCellIndex(i - 1, j, k);
			bound.SetGhostIndex(ghostNodeIndex);
			boundNorm = cell_center[ghostNodeIndex] - cell_center[cellIndex];
			bound.SetNorm(boundNorm);
			boundMap->AddBoundary("inlet", bound);

			i = xNodeNum - 3;
			cellIndex = grid->GetCellIndex(i, j, k);
			bound.SetNodeIndex(cellIndex);
			innerNodeIndex = grid->GetCellIndex(i - 1, j, k);
			bound.SetInnerIndex(innerNodeIndex);
			ghostNodeIndex = grid->GetCellIndex(i + 1, j, k);
			bound.SetGhostIndex(ghostNodeIndex);
			boundNorm = cell_center[ghostNodeIndex] - cell_center[cellIndex];
			bound.SetNorm(boundNorm);
			boundMap->AddBoundary("outlet", bound);
		}
	}

	//其余四个面均为壁面
	for (i = 1; i < xNodeNum - 1; ++i)
	{
		//j方向
		for (k = 1; k < zNodeNum - 1; ++k)
		{
			j = 1;
			cellIndex = grid->GetCellIndex(i, j, k);
			bound.SetNodeIndex(cellIndex);
			innerNodeIndex = grid->GetCellIndex(i, j + 1, k);
			bound.SetInnerIndex(innerNodeIndex);
			ghostNodeIndex = grid->GetCellIndex(i, j - 1, k);
			bound.SetGhostIndex(ghostNodeIndex);
			boundNorm = cell_center[ghostNodeIndex] - cell_center[cellIndex];
			bound.SetNorm(boundNorm);
			boundMap->AddBoundary("outlet", bound);
			j = yNodeNum - 3;
			cellIndex = grid->GetCellIndex(i, j, k);
			bound.SetNodeIndex(cellIndex);
			innerNodeIndex = grid->GetCellIndex(i, j - 1, k);
			bound.SetInnerIndex(innerNodeIndex);
			ghostNodeIndex = grid->GetCellIndex(i, j + 1, k);
			bound.SetGhostIndex(ghostNodeIndex);
			boundNorm = cell_center[ghostNodeIndex] - cell_center[cellIndex];
			bound.SetNorm(boundNorm);
			boundMap->AddBoundary("outlet", bound);
		}

		//j方向
		for (j = 1; j < yNodeNum - 1; ++j)
		{
			k = 1;
			cellIndex = grid->GetCellIndex(i, j, k);
			bound.SetNodeIndex(cellIndex);
			innerNodeIndex = grid->GetCellIndex(i, j, k + 1);
			bound.SetInnerIndex(innerNodeIndex);
			ghostNodeIndex = grid->GetCellIndex(i, j, k - 1);
			bound.SetGhostIndex(ghostNodeIndex);
			boundNorm = cell_center[ghostNodeIndex] - cell_center[cellIndex];
			bound.SetNorm(boundNorm);
			boundMap->AddBoundary("outlet", bound);
			k = zNodeNum - 3;
			cellIndex = grid->GetCellIndex(i, j, k);
			bound.SetNodeIndex(cellIndex);
			innerNodeIndex = grid->GetCellIndex(i, j, k - 1);
			bound.SetInnerIndex(innerNodeIndex);
			ghostNodeIndex = grid->GetCellIndex(i, j, k + 1);
			bound.SetGhostIndex(ghostNodeIndex);
			boundNorm = cell_center[ghostNodeIndex] - cell_center[cellIndex];
			bound.SetNorm(boundNorm);
			boundMap->AddBoundary("outlet", bound);
		}
	}
	gridList->AddGrid(grid);
}

void zaran::GridListFactoryZaran3D::TagCell(Ptr<GridList>& gridList)
{
	Ptr<zaran::Grid_Zaran_3D> grid = std::static_pointer_cast<Grid_Zaran_3D>(gridList->GetGrid(0));
	auto cellTopo = std::static_pointer_cast<CellTopoInfoZaran>(grid->GetCellTopo());
	auto& cell_type = cellTopo->GetType();
	auto& cell_center = cellTopo->GetCenterCoord();
	int ni, nj, nk;
	grid->GetNodeNum(ni, nj, nk);
	cell_type.resize((ni - 1) * (nj - 1) * (nk - 1));
	int i, j, k;
	double x_min, x_max, y_min, y_max, z_min, z_max;
	grid->GetBox(x_min, x_max, y_min, y_max, z_min, z_max);
	double dx = (x_max - x_min) / (ni - 1);
	double dy = (y_max - y_min) / (nj - 1);
	double dz = (z_max - z_min) / (nk - 1);
	DArray dist_cell_to_model;//单元中心到模型的距离
	dist_cell_to_model.resize((cell_type.size()));
	IArray cell_in_model;//单元中心是否在模型内
	cell_in_model.resize((cell_type.size()));
	ZaranLog::info("Start tag cell type");
	auto poly_data = m_polyDataModel->GetPolyData();
	vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
	implicitPolyDataDistance->SetInput(poly_data);
#pragma omp parallel for
	for (int iCell = 0;iCell < cell_type.size();++iCell)
	{
		dist_cell_to_model[iCell] = abs(implicitPolyDataDistance->FunctionValue(cell_center[iCell].data()));
	}
	double tol = 0.5 * sqrt(dx * dx + dy * dy + dz * dz);
	for (int iCell = 0;iCell < cell_type.size();++iCell)
	{
		if (dist_cell_to_model[iCell] < tol)
		{
			cell_type[iCell] = CellType::Solid;
		}
		else
		{
			cell_type[iCell] = CellType::Unknown;
		}
	}
	int nFluid = 1;
	cell_type[0] = CellType::Fluid;
	while (nFluid != 0)
	{
		nFluid = 0;
		for (int iCell = 0;iCell < cell_type.size();++iCell)
		{
			if (cell_type[iCell] == CellType::Fluid)
			{
				grid->GetCellIndex(iCell, i, j, k);
				if (i > 0 && cell_type[grid->GetCellIndex(i - 1, j, k)] == CellType::Unknown)
				{
					cell_type[grid->GetCellIndex(i - 1, j, k)] = CellType::Fluid;
					nFluid++;
				}
				if (i < ni - 2 && cell_type[grid->GetCellIndex(i + 1, j, k)] == CellType::Unknown)
				{
					cell_type[grid->GetCellIndex(i + 1, j, k)] = CellType::Fluid;
					nFluid++;
				}
				if (j > 0 && cell_type[grid->GetCellIndex(i, j - 1, k)] == CellType::Unknown)
				{
					cell_type[grid->GetCellIndex(i, j - 1, k)] = CellType::Fluid;
					nFluid++;
				}
				if (j < nj - 2 && cell_type[grid->GetCellIndex(i, j + 1, k)] == CellType::Unknown)
				{
					cell_type[grid->GetCellIndex(i, j + 1, k)] = CellType::Fluid;
					nFluid++;
				}
				if (k > 0 && cell_type[grid->GetCellIndex(i, j, k - 1)] == CellType::Unknown)
				{
					cell_type[grid->GetCellIndex(i, j, k - 1)] = CellType::Fluid;
					nFluid++;
				}
				if (k < nk - 2 && cell_type[grid->GetCellIndex(i, j, k + 1)] == CellType::Unknown)
				{
					cell_type[grid->GetCellIndex(i, j, k + 1)] = CellType::Fluid;
					nFluid++;
				}
			}
		}
		ZaranLog::info("nFluid: {}", nFluid);
	}
	for (int iCell = 0;iCell < cell_type.size();++iCell)
	{
		if (cell_type[iCell] == CellType::Unknown)
			cell_type[iCell] = CellType::Solid;
	}
	int n_new_solid = 1;
	int iCell = 0;
	while (n_new_solid != 0)
	{
		n_new_solid = 0;
		for (k = 1; k < nk - 1; ++k)
		{
			for (j = 1; j < nj - 1; ++j)
			{
				for (i = 1; i < ni - 1; ++i)
				{
					iCell = grid->GetCellIndex(i, j, k);
					if (cell_type[iCell] != CellType::Fluid)
						continue;
					if (cell_type[grid->GetCellIndex(i - 1, j, k)] == CellType::Solid &&
						cell_type[grid->GetCellIndex(i + 1, j, k)] == CellType::Solid ||
						cell_type[grid->GetCellIndex(i, j - 1, k)] == CellType::Solid &&
						cell_type[grid->GetCellIndex(i, j + 1, k)] == CellType::Solid ||
						cell_type[grid->GetCellIndex(i, j, k - 1)] == CellType::Solid &&
						cell_type[grid->GetCellIndex(i, j, k + 1)] == CellType::Solid)
					{
						cell_type[iCell] = CellType::Solid;
						n_new_solid++;
					}
				}
			}
		}
		ZaranLog::info("n_new_solid: {}", n_new_solid);
	}

	ZaranLog::info("Tag cell type finished");
	ZaranLog::info("Start find fluid-solid interface");
	for (int iCell = 0;iCell < cell_type.size();++iCell)
	{
		grid->GetCellIndex(iCell, i, j, k);
		if (cell_type[iCell] == CellType::Fluid)
		{
			if (i > 0 && cell_type[grid->GetCellIndex(i - 1, j, k)] == CellType::Solid)
			{
				cell_type[iCell] = CellType::FluidSolid;
				continue;
			}
			if (i < ni - 2 && cell_type[grid->GetCellIndex(i + 1, j, k)] == CellType::Solid)
			{
				cell_type[iCell] = CellType::FluidSolid;
				continue;
			}
			if (j > 0 && cell_type[grid->GetCellIndex(i, j - 1, k)] == CellType::Solid)
			{
				cell_type[iCell] = CellType::FluidSolid;
				continue;
			}
			if (j < nj - 2 && cell_type[grid->GetCellIndex(i, j + 1, k)] == CellType::Solid)
			{
				cell_type[iCell] = CellType::FluidSolid;
				continue;
			}
			if (k > 0 && cell_type[grid->GetCellIndex(i, j, k - 1)] == CellType::Solid)
			{
				cell_type[iCell] = CellType::FluidSolid;
				continue;
			}
			if (k < nk - 2 && cell_type[grid->GetCellIndex(i, j, k + 1)] == CellType::Solid)
			{
				cell_type[iCell] = CellType::FluidSolid;
				continue;
			}
		}
	}
	int n_bad_cell = 0;
	for (int iCell = 0;iCell < cell_type.size();++iCell)
	{
		if (cell_type[iCell] != CellType::FluidSolid)
			continue;
		if (i > 0 && i < ni - 2)
		{
			if (cell_type[grid->GetCellIndex(i - 1, j, k)] == CellType::Solid && cell_type[grid->GetCellIndex(i + 1, j, k)] == CellType::Solid)
			{
				n_bad_cell++;
				continue;
			}
		}
		if (j > 0 && j < nj - 2)
		{
			if (cell_type[grid->GetCellIndex(i, j - 1, k)] == CellType::Solid && cell_type[grid->GetCellIndex(i, j + 1, k)] == CellType::Solid)
			{
				n_bad_cell++;
				continue;
			}
		}
		if (k > 0 && k < nk - 2)
		{
			if (cell_type[grid->GetCellIndex(i, j, k - 1)] == CellType::Solid && cell_type[grid->GetCellIndex(i, j, k + 1)] == CellType::Solid)
			{
				n_bad_cell++;
				continue;
			}
		}
	}
	ZaranLog::info("n_bad_cell: {}", n_bad_cell);
	ZaranLog::info("Find fluid-solid interface finished");
	ZaranLog::info("Output cell type");
	std::ofstream fout("cell_type.dat");
	fout << "variables=x,y,z,label\n";
	fout << "ZONE T= grid_" << grid->GetName() << std::endl;
	fout << "I=" << ni << ", J=" << nj << ", K=" << nk << ", DATAPACKING=BLOCK, VARLOCATION=([4]=CELLCENTERED)" << std::endl;
	int count = 0;
	for (k = 0; k < nk; ++k)
	{
		for (j = 0; j < nj; ++j)
		{
			for (i = 0; i < ni; ++i)
			{
				fout << x_min + i * dx << " ";
				if (++count % 10 == 0)
				{
					count = 0;
					fout << std::endl;
				}
			}
		}
	}
	for (k = 0; k < nk; ++k)
	{
		for (j = 0; j < nj; ++j)
		{
			for (i = 0; i < ni; ++i)
			{
				fout << y_min + j * dy << " ";
				if (++count % 10 == 0)
				{
					count = 0;
					fout << std::endl;
				}
			}
		}
	}
	for (k = 0; k < nk; ++k)
	{
		for (j = 0; j < nj; ++j)
		{
			for (i = 0; i < ni; ++i)
			{
				fout << z_min + k * dz << " ";
				if (++count % 10 == 0)
				{
					count = 0;
					fout << std::endl;
				}
			}
		}
	}
	for (k = 0; k < nk - 1; ++k)
	{
		for (j = 0; j < nj - 1; ++j)
		{
			for (i = 0; i < ni - 1; ++i)
			{
				fout << int(cell_type[grid->GetCellIndex(i, j, k)]) << " ";
				// fout << dist_cell_to_model[grid->GetCellIndex(i, j, k)] << " ";

				if (++count % 10 == 0)
				{
					count = 0;
					fout << std::endl;
				}
			}
		}
	}
	fout.close();
	ZaranLog::info("Output cell type finished");

}

void zaran::GridListFactoryZaran3D::CrateBoundPatch(Ptr<GridList>& gridList)
{
	Ptr<zaran::Grid_Zaran_3D> grid = std::static_pointer_cast<Grid_Zaran_3D>(gridList->GetGrid(0));
	auto cellTopo = std::static_pointer_cast<CellTopoInfoZaran>(grid->GetCellTopo());
	auto& cell_type = cellTopo->GetType();
	auto& cell_center = cellTopo->GetCenterCoord();
	int ni, nj, nk;
	grid->GetNodeNum(ni, nj, nk);
	int i, j, k;
	double x_min, x_max, y_min, y_max, z_min, z_max;
	grid->GetBox(x_min, x_max, y_min, y_max, z_min, z_max);
	double dx = (x_max - x_min) / (ni - 1);
	double dy = (y_max - y_min) / (nj - 1);
	double dz = (z_max - z_min) / (nk - 1);
	auto poly_data = m_polyDataModel->GetPolyData();
	vtkNew<vtkCellLocator> cellLocator;
	cellLocator->SetDataSet(poly_data);
	cellLocator->BuildLocator();
	ZaranBoundPatch& bound_patch = grid->GetBoundPatch();
	int n_bound_patch = 0;
	for (int iCell = 0;iCell < cell_type.size();iCell++)
	{
		if (cell_type[iCell] == CellType::FluidSolid)
			n_bound_patch++;
	}
	auto& mid_index = bound_patch.GetIndex();
	auto& bound_coord = bound_patch.GetCoordinate();
	auto& bound_norm = bound_patch.GetNormal();
	mid_index.resize(n_bound_patch);
	bound_coord.resize(n_bound_patch);
	bound_norm.resize(n_bound_patch);
	double closestPoint[3];   // the coordinates of the closest point will be
	// returned here
	double closestPointDist2; // the squared distance to the closest point will be
	// returned here
	vtkIdType cellId; // the cell id of the cell containing the closest point will
	// be returned here
	int subId;        // this is rarely used (in triangle strips only, I believe)
	int iBound = 0;
	double tol = 0.5 * sqrt(dx * dx + dy * dy + dz * dz);
	double dist;
	for (int iCell = 0;iCell < cell_type.size();iCell++)
	{
		if (cell_type[iCell] != CellType::FluidSolid)
			continue;
		grid->GetCellIndex(iCell, i, j, k);
		cellLocator->FindClosestPoint(cell_center[iCell].data(), closestPoint, cellId, subId,
			closestPointDist2);
		mid_index[iBound] = { i,j,k };
		bound_coord[iBound] = { closestPoint[0],closestPoint[1],closestPoint[2] };
		bound_norm[iBound] = { cell_center[iCell][0] - closestPoint[0],cell_center[iCell][1] - closestPoint[1],cell_center[iCell][2] - closestPoint[2] };
		dist = sqrt(bound_norm[iBound][0] * bound_norm[iBound][0] + bound_norm[iBound][1] * bound_norm[iBound][1] + bound_norm[iBound][2] * bound_norm[iBound][2]);
		bound_norm[iBound].normalize();
		iBound++;
	}
	std::ofstream fout("bound_patch.dat");
	fout << "variables=x,y,z\n";
	for (int iBound = 0;iBound < n_bound_patch;++iBound)
	{
		fout << bound_coord[iBound][0] << " " << bound_coord[iBound][1] << " " << bound_coord[iBound][2] << std::endl;
	}
}
void zaran::GridListFactoryZaran3D::ReadModel()
{
	string modelFileName = GlobalData::GetString("modelFileName");
	STLReader reader;
	reader.ReadSTLFile(modelFileName.c_str());
	m_polyDataModel = std::make_shared<PolyDataModel>();
	m_polyDataModel->SetPolyData(reader.GetMesh(), 1e-6);
	if (m_polyDataModel->IsClosed())
		ZaranLog::info("Import Model: {}, is closed!", modelFileName);
	else
		ZaranLog::info("Import Model: {}, is not closed!", modelFileName);
}
