#include "Visual.h"
#include "GlobalData.h"
#include "Log.h"
#include <TECIO.h>
#include <fstream>
#include <string>
#include <vector>
#include <cgnslib.h>
#include "NSFieldFN.h"
using namespace zaran;
void zaran::Visual::WriteTecplotBinary(shared_ptr<NSFieldFNFDM> field)
{
	auto data_manager = field->GetDataManager();
	auto grid = field->GetGrid();
	CellFN& cell = grid->GetCell();
	NodeFN& node = grid->GetNode();

	const double* density, * velocity_x, * velocity_y, * velocity_z, * pressure;
	density = data_manager->GetPrim(ID_DENSITY);
	velocity_x = data_manager->GetPrim(ID_VELOCITY_X);
	velocity_y = data_manager->GetPrim(ID_VELOCITY_Y);
	velocity_z = data_manager->GetPrim(ID_VELOCITY_Z);
	pressure = data_manager->GetPrim(ID_PRESSURE);

	INTEGER4 node_num = grid->GetTotalNodeNum();
	INTEGER4 cell_num = cell.GetCellNum();
	dynamic_array<double> x(node_num), y(node_num), z(node_num);
	for (int iNode = 0; iNode < node_num; ++iNode)
	{
		x[iNode] = node.GetCoord(iNode)[0];
		y[iNode] = node.GetCoord(iNode)[1];
		z[iNode] = node.GetCoord(iNode)[2];
	}
	INTEGER4 vIsDouble = 1;
	string zone_name = "grid_" + std::to_string(field->GetIdx());
	INTEGER4 zone_type = 5; // Brick
	INTEGER4 face_num = 6;
	INTEGER4 iCellMax = 0;
	INTEGER4 jCellMax = 0;
	INTEGER4 kCellMax = 0;
	double solution_time = GlobalData::GetDouble("currentTime");
	INTEGER4 strandID = 0;
	INTEGER4 parentZn = 0;
	INTEGER4 isBlock = 1;
	INTEGER4 nFConns = 0;
	INTEGER4 FNMode = 0;
	int valueLocation[] = { 1, 1, 1, 1, 1, 1, 1, 1 };
	int shrConn = 0;
	int i = TECZNE142((char*)zone_name.c_str(), &zone_type, &node_num, &cell_num, &face_num, &iCellMax, &jCellMax,
		&kCellMax, &solution_time, &strandID, &parentZn, &isBlock, &nFConns, &FNMode, 0, 0, 0, NULL,
		valueLocation, NULL, &shrConn);

	i = TECDAT142(&node_num, x.data(), &vIsDouble);
	i = TECDAT142(&node_num, y.data(), &vIsDouble);
	i = TECDAT142(&node_num, z.data(), &vIsDouble);
	i = TECDAT142(&node_num, density, &vIsDouble);
	i = TECDAT142(&node_num, velocity_x, &vIsDouble);
	i = TECDAT142(&node_num, velocity_y, &vIsDouble);
	i = TECDAT142(&node_num, velocity_z, &vIsDouble);
	i = TECDAT142(&node_num, pressure, &vIsDouble);
	INTEGER4 connectivityCount = cell_num * 8;
	dynamic_array<INTEGER4> cell_nodes(connectivityCount);
	for (int iCell = 0; iCell < cell_num; ++iCell)
	{
		int node_num = cell.GetNodeNum(iCell);
		auto cell2node = cell.GetNode(iCell);
		for (int iNode = 0; iNode < node_num; ++iNode)
		{
			cell_nodes[iCell * 8 + iNode] = cell2node[iNode] + 1;
		}
	}

	i = TECNODE142(&connectivityCount, cell_nodes.data());

	/// bound face
	auto face_topo = grid->GetFace();
	cell_num = face_topo.GetFaceNum();
	if (cell_num == 0)
	{
		return;
	}
	zone_name = "grid_" + grid->GetName() + "_bound";
	zone_type = 3; // Brick
	face_num = 6;
	iCellMax = 0;
	jCellMax = 0;
	kCellMax = 0;
	strandID = 0;
	parentZn = 0;
	isBlock = 1;
	nFConns = 0;
	FNMode = 0;
	i = TECZNE142((char*)zone_name.c_str(), &zone_type, &node_num, &cell_num, &face_num, &iCellMax, &jCellMax,
		&kCellMax, &solution_time, &strandID, &parentZn, &isBlock, &nFConns, &FNMode, 0, 0, 0, NULL,
		valueLocation, NULL, &shrConn);

	i = TECDAT142(&node_num, x.data(), &vIsDouble);
	i = TECDAT142(&node_num, y.data(), &vIsDouble);
	i = TECDAT142(&node_num, z.data(), &vIsDouble);
	i = TECDAT142(&node_num, density, &vIsDouble);
	i = TECDAT142(&node_num, velocity_x, &vIsDouble);
	i = TECDAT142(&node_num, velocity_y, &vIsDouble);
	i = TECDAT142(&node_num, velocity_z, &vIsDouble);
	i = TECDAT142(&node_num, pressure, &vIsDouble);
	int node_num_per_cell = 4;
	connectivityCount = cell_num * node_num_per_cell;
	dynamic_array<INTEGER4> face_nodes(connectivityCount);
	for (int iFace = 0; iFace < cell_num; ++iFace)
	{
		auto face2node = face_topo.GetFace2Node(iFace);
		int n_node = face_topo.GetFaceNodeNum(iFace);
		for (int iNode = 0; iNode < n_node; ++iNode)
		{
			face_nodes[iFace * node_num_per_cell + iNode] = face2node[iNode] + 1;
		}
		if (n_node < 4)
		{
			for (int i = n_node; i < 4; ++i)
			{
				face_nodes[iFace * node_num_per_cell + i] = face2node[0] + 1;
			}
		}
	}

	i = TECNODE142(&connectivityCount, face_nodes.data());

}

void Visual::WriteCGNS(shared_ptr<NSFieldFNFDM> field, cgsize_t index_file1, cgsize_t index_base1)
{
	std::string work_dir = GlobalData::GetString("work_dir");
	std::string file_name = "result\\" + std::to_string(GlobalData::GetInt("currentIter")) + "_bound" + ".cgns";
	file_name = work_dir + "\\" + file_name;
	int index_file;
	cg_open(file_name.c_str(), CG_MODE_WRITE, &index_file);
	int cell_dim = 3;
	int phys_dim = 3;
	int index_base;
	cg_base_write(index_file, "Base", cell_dim, phys_dim, &index_base);

	auto data_manager = field->GetDataManager();
	auto grid = field->GetGrid();
	auto cell = grid->GetCell();
	auto node = grid->GetNode();

	const double* density, * velocity_x, * velocity_y, * velocity_z, * pressure;
	density = data_manager->GetPrim(ID_DENSITY);
	velocity_x = data_manager->GetPrim(ID_VELOCITY_X);
	velocity_y = data_manager->GetPrim(ID_VELOCITY_Y);
	velocity_z = data_manager->GetPrim(ID_VELOCITY_Z);
	pressure = data_manager->GetPrim(ID_PRESSURE);

	int node_num = grid->GetTotalNodeNum();
	int cell_num = cell.GetCellNum();
	dynamic_array<double> x(node_num), y(node_num), z(node_num);
	for (int iNode = 0; iNode < node_num; ++iNode)
	{
		x[iNode] = node.GetCoord(iNode)[0];
		y[iNode] = node.GetCoord(iNode)[1];
		z[iNode] = node.GetCoord(iNode)[2];
	}
	int index_zone;
	cgsize_t isize[3];
	isize[0] = node_num;
	isize[1] = cell_num;
	isize[2] = 0;
	cg_zone_write(index_file, index_base, "FN_Zone", isize, CGNS_ENUMV(Unstructured), &index_zone);
	int index_coord;
	cg_coord_write(index_file, index_base, index_zone, CGNS_ENUMV(RealDouble), "CoordinateX", x.data(), &index_coord);
	cg_coord_write(index_file, index_base, index_zone, CGNS_ENUMV(RealDouble), "CoordinateY", y.data(), &index_coord);
	cg_coord_write(index_file, index_base, index_zone, CGNS_ENUMV(RealDouble), "CoordinateZ", z.data(), &index_coord);
	dynamic_array<cgsize_t>elements(cell_num * 8);
	for (int iCell = 0; iCell < cell_num; ++iCell)
	{
		int cell_node_num = cell.GetNodeNum(iCell);
		auto cell2node = cell.GetNode(iCell);
		for (int iNode = 0; iNode < cell_node_num; ++iNode)
		{
			elements[iCell * 8 + iNode] = cell2node[iNode] + 1;
		}
	}
	int index_section;
	cg_section_write(index_file, index_base, index_zone, "Section", CGNS_ENUMV(HEXA_8), 1, cell_num, 0, elements.data(), &index_section);
	int index_field;
	cg_sol_write(index_file, index_base, index_zone, "FlowSolution", CGNS_ENUMV(Vertex), &index_field);
	cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "Density", density, &index_field);
	cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "VelocityX", velocity_x, &index_field);
	cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "VelocityY", velocity_y, &index_field);
	cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "VelocityZ", velocity_z, &index_field);
	cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "Pressure", pressure, &index_field);
	//输出边界面
	int index_file2, index_base2;
	file_name = "result\\" + std::to_string(GlobalData::GetInt("currentIter")) + "_bound_face" + ".cgns";
	file_name = work_dir + "\\" + file_name;
	cg_open(file_name.c_str(), CG_MODE_WRITE, &index_file2);
	cg_base_write(index_file2, "Base", cell_dim, phys_dim, &index_base2);
	auto face_topo = grid->GetFace();
	int face_num = face_topo.GetFaceNum();
	dynamic_array<cgsize_t> faces(face_num * 4);
	for (int iFace = 0; iFace < face_num; ++iFace)
	{
		auto face2node = face_topo.GetFace2Node(iFace);
		for (int iNode = 0; iNode < 4; ++iNode)
		{
			faces[iFace * 4 + iNode] = face2node[iNode] + 1;
		}
	}
	int index_zone1;
	cgsize_t isize1[3];
	isize1[0] = node_num;
	isize1[1] = face_num;
	isize1[2] = 0;
	cg_zone_write(index_file2, index_base2, "bound_Zone", isize1, CGNS_ENUMV(Unstructured), &index_zone1);
	int index_coord1;
	cg_coord_write(index_file2, index_base2, index_zone1, CGNS_ENUMV(RealDouble), "CoordinateX", x.data(), &index_coord1);
	cg_coord_write(index_file2, index_base2, index_zone1, CGNS_ENUMV(RealDouble), "CoordinateY", y.data(), &index_coord1);
	cg_coord_write(index_file2, index_base2, index_zone1, CGNS_ENUMV(RealDouble), "CoordinateZ", z.data(), &index_coord1);
	int index_section1;
	cg_section_write(index_file2, index_base2, index_zone1, "bound_Section", CGNS_ENUMV(QUAD_4), 1, face_num, 0, faces.data(), &index_section1);
	int index_field1;
	cg_sol_write(index_file2, index_base2, index_zone1, "bound_FlowSolution", CGNS_ENUMV(Vertex), &index_field1);
	cg_field_write(index_file2, index_base2, index_zone1, index_field1, CGNS_ENUMV(RealDouble), "Density", density, &index_field1);
	cg_field_write(index_file2, index_base2, index_zone1, index_field1, CGNS_ENUMV(RealDouble), "VelocityX", velocity_x, &index_field1);
	cg_field_write(index_file2, index_base2, index_zone1, index_field1, CGNS_ENUMV(RealDouble), "VelocityY", velocity_y, &index_field1);
	cg_field_write(index_file2, index_base2, index_zone1, index_field1, CGNS_ENUMV(RealDouble), "VelocityZ", velocity_z, &index_field1);
	cg_field_write(index_file2, index_base2, index_zone1, index_field1, CGNS_ENUMV(RealDouble), "Pressure", pressure, &index_field1);
}

void Visual::WriteCGNS(shared_ptr<NSFieldZaran> field, cgsize_t index_file1, cgsize_t index_base1)
{
	std::string work_dir = GlobalData::GetString("work_dir");
	std::string file_name = "result\\" + std::to_string(GlobalData::GetInt("currentIter")) + "_back" + ".cgns";
	file_name = work_dir + "\\" + file_name;
	int index_file;
	cg_open(file_name.c_str(), CG_MODE_WRITE, &index_file);
	int cell_dim = 3;
	int phys_dim = 3;
	int index_base;
	cg_base_write(index_file, "Base", cell_dim, phys_dim, &index_base);

	auto grid = field->GetGrid();
	auto node = grid->GetNode();
	auto cell = grid->GetCell();
	auto data_manager = field->GetDataManager();
	auto idx_proxy = grid->GetIdxProxy();
	int is, ie, js, je, ks, ke;
	grid->GetRange(is, ie, js, je, ks, ke);
	int ni = ie - is + 1;
	int nj = je - js + 1;
	int nk = ke - ks + 1;
	int grid_ni = grid->GetNi();
	int grid_nj = grid->GetNj();
	int grid_nk = grid->GetNk();
	int node_num = ni * nj * nk;
	// 创建 zone
	cgsize_t isize[3][3];
	isize[0][0] = ni;
	isize[0][1] = nj;
	isize[0][2] = nk;
	isize[1][0] = isize[0][0] - 1;
	isize[1][1] = isize[0][1] - 1;
	isize[1][2] = isize[0][2] - 1;
	isize[2][0] = 0;
	isize[2][1] = 0;
	isize[2][2] = 0;
	int index_zone;
	cg_zone_write(index_file, index_base, "FN_Zone", *isize, CGNS_ENUMV(Structured), &index_zone);
	dynamic_array<double> x(node_num), y(node_num), z(node_num);
	dynamic_array<double> density(node_num), velocity_x(node_num), velocity_y(node_num), velocity_z(node_num), pressure(node_num);
	for (index_type k = 0; k < nk; ++k)
	{
		for (index_type j = 0; j < nj; ++j)
		{
			for (index_type i = 0; i < ni; ++i)
			{
				index_type idx = i + ni * j + ni * nj * k;
				auto coord = node->GetCoord(i + is, j + js, k + ks);
				x[idx] = coord[0];
				y[idx] = coord[1];
				z[idx] = coord[2];
				int idx0 = idx_proxy->GetIdx(i + is, j + js, k + ks);
				density[idx] = data_manager->GetPrim(ID_DENSITY, idx0);
				velocity_x[idx] = data_manager->GetPrim(ID_VELOCITY_X, idx0);
				velocity_y[idx] = data_manager->GetPrim(ID_VELOCITY_Y, idx0);
				velocity_z[idx] = data_manager->GetPrim(ID_VELOCITY_Z, idx0);
				pressure[idx] = data_manager->GetPrim(ID_PRESSURE, idx0);
			}
		}
	}
	int index_coord;
	cg_coord_write(index_file, index_base, index_zone, CGNS_ENUMV(RealDouble), "CoordinateX", x.data(), &index_coord);
	cg_coord_write(index_file, index_base, index_zone, CGNS_ENUMV(RealDouble), "CoordinateY", y.data(), &index_coord);
	cg_coord_write(index_file, index_base, index_zone, CGNS_ENUMV(RealDouble), "CoordinateZ", z.data(), &index_coord);

	int index_field;
	cg_sol_write(index_file, index_base, index_zone, "FlowSolution", CGNS_ENUMV(Vertex), &index_field);
	cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "Density", density.data(), &index_field);
	cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "VelocityX", velocity_x.data(), &index_field);
	cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "VelocityY", velocity_y.data(), &index_field);
	cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "VelocityZ", velocity_z.data(), &index_field);
	cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "Pressure", pressure.data(), &index_field);
	cg_close(index_file);



}

void Visual::WriteCGNS(shared_ptr<NSFieldStruct> field, cgsize_t index_file1, cgsize_t index_base1)
{
	std::string work_dir = GlobalData::GetString("work_dir");
	std::string file_name = "result\\" + std::to_string(GlobalData::GetInt("currentIter")) + "_back" + ".cgns";
	file_name = work_dir + "\\" + file_name;
	int index_file;
	cg_open(file_name.c_str(), CG_MODE_WRITE, &index_file);
	int cell_dim = 3;
	int phys_dim = 3;
	int index_base;
	cg_base_write(index_file, "Base", cell_dim, phys_dim, &index_base);

	auto grid = field->GetGrid();
	auto node = grid->GetNode();
	auto cell = grid->GetCell();
	auto data_manager = field->GetDataManager();
	auto idx_proxy = grid->GetIdxProxy();
	int is, ie, js, je, ks, ke;
	grid->GetRange(is, ie, js, je, ks, ke);
	int ni = ie - is + 1;
	int nj = je - js + 1;
	int nk = ke - ks + 1;
	int grid_ni = grid->GetNi();
	int grid_nj = grid->GetNj();
	int grid_nk = grid->GetNk();
	int node_num = ni * nj * nk;
	// 创建 zone
	cgsize_t isize[3][3];
	isize[0][0] = nk;
	isize[0][1] = nj;
	isize[0][2] = ni;
	isize[1][0] = isize[0][0] - 1;
	isize[1][1] = isize[0][1] - 1;
	isize[1][2] = isize[0][2] - 1;
	isize[2][0] = 0;
	isize[2][1] = 0;
	isize[2][2] = 0;
	int index_zone;
	cg_zone_write(index_file, index_base, "Structured_Zone", *isize, CGNS_ENUMV(Structured), &index_zone);
	dynamic_array<double> x(node_num), y(node_num), z(node_num);
	dynamic_array<double> density(node_num), velocity_x(node_num), velocity_y(node_num), velocity_z(node_num), pressure(node_num);
	for (index_type k = 0; k < nk; ++k)
	{
		for (index_type j = 0; j < nj; ++j)
		{
			for (index_type i = 0; i < ni; ++i)
			{
				index_type idx = i + ni * j + ni * nj * k;
				x[idx] = node->GetCoord(i + is, j + js, k + ks)[0];
				y[idx] = node->GetCoord(i + is, j + js, k + ks)[1];
				z[idx] = node->GetCoord(i + is, j + js, k + ks)[2];
				int idx0 = idx_proxy->GetIdx(i + is, j + js, k + ks);
				density[idx] = data_manager->GetPrim(ID_DENSITY, idx0);
				velocity_x[idx] = data_manager->GetPrim(ID_VELOCITY_X, idx0);
				velocity_y[idx] = data_manager->GetPrim(ID_VELOCITY_Y, idx0);
				velocity_z[idx] = data_manager->GetPrim(ID_VELOCITY_Z, idx0);
				pressure[idx] = data_manager->GetPrim(ID_PRESSURE, idx0);
			}
		}
	}
	int index_coord;
	cg_coord_write(index_file, index_base, index_zone, CGNS_ENUMV(RealDouble), "CoordinateX", x.data(), &index_coord);
	cg_coord_write(index_file, index_base, index_zone, CGNS_ENUMV(RealDouble), "CoordinateY", y.data(), &index_coord);
	cg_coord_write(index_file, index_base, index_zone, CGNS_ENUMV(RealDouble), "CoordinateZ", z.data(), &index_coord);

	int index_field;
	cg_sol_write(index_file, index_base, index_zone, "FlowSolution", CGNS_ENUMV(Vertex), &index_field);
	cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "Density", density.data(), &index_field);
	cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "VelocityX", velocity_x.data(), &index_field);
	cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "VelocityY", velocity_y.data(), &index_field);
	cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "VelocityZ", velocity_z.data(), &index_field);
	cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "Pressure", pressure.data(), &index_field);
	cg_close(index_file);

}

void Visual::WriteCGNS(shared_ptr<FieldManager> field_manager)
{
	//std::string work_dir = GlobalData::GetString("work_dir");
	//std::string file_name = "result\\" + std::to_string(GlobalData::GetInt("currentIter")) + ".cgns";
	//file_name = work_dir + "\\" + file_name;
	int index_file;
	//cg_open(file_name.c_str(), CG_MODE_WRITE, &index_file);
	int index_base;
	//int cell_dim = 3;
	//int phys_dim = 3;
	//cg_base_write(index_file, "Base", cell_dim, phys_dim, &index_base);
	for (size_t iter_field = 0; iter_field < field_manager->GetFieldNum(); iter_field++)
	{
		auto field = field_manager->GetField(iter_field);
		auto field_type = field->GetFieldType();
		if (field_type == FieldType::NS_Structured)
		{
			auto field_struct = std::dynamic_pointer_cast<NSFieldStruct>(field);
			WriteCGNS(field_struct, index_file, index_base);
		}
		else if (field_type == FieldType::NS_Zaran)
		{
			auto field_zaran = std::dynamic_pointer_cast<NSFieldZaran>(field);
			WriteCGNS(field_zaran, index_file, index_base);
		}
		else if (field_type == FieldType::NS_FlexibleNode)
		{
			auto field_fn = std::dynamic_pointer_cast<NSFieldFNFDM>(field);
			WriteCGNS(field_fn, index_file, index_base);
		}
		else
		{
			Log::warn("Field type is not supported!");
		}
	}
	//cg_close(index_file);



}

void Visual::WriteVtkASCII(shared_ptr<NSFieldZaran> field, std::ostream& os)
{
	auto grid = field->GetGrid();
	auto node = grid->GetNode();
	auto data_manager = field->GetDataManager();
	auto idx_proxy = grid->GetIdxProxy();

	int is, ie, js, je, ks, ke;
	grid->GetRange(is, ie, js, je, ks, ke);
	int ni = ie - is + 1;
	int nj = je - js + 1;
	int nk = ke - ks + 1;
	int node_num = ni * nj * nk;

	// 写入 VTK 文件头
	os << "# vtk DataFile Version 3.0\n";
	os << "Flow Field Data\n";
	os << "ASCII\n";
	os << "DATASET STRUCTURED_GRID\n";
	os << "DIMENSIONS " << ni << " " << nj << " " << nk << "\n";
	os << "POINTS " << node_num << " double\n";

	// 写入点坐标
	for (int k = ks; k <= ke; ++k)
	{
		for (int j = js; j <= je; ++j)
		{
			for (int i = is; i <= ie; ++i)
			{
				auto coord = node->GetCoord(i, j, k);
				os << coord[0] << " " << coord[1] << " " << coord[2] << "\n";
			}
		}
	}

	// 写入点数据
	os << "\nPOINT_DATA " << node_num << "\n";

	// 写入标量数据：密度
	os << "SCALARS Density double 1\n";
	os << "LOOKUP_TABLE default\n";
	for (int k = ks; k <= ke; ++k)
	{
		for (int j = js; j <= je; ++j)
		{
			for (int i = is; i <= ie; ++i)
			{
				int idx = idx_proxy->GetIdx(i, j, k);
				double density = data_manager->GetPrim(ID_DENSITY, idx);
				os << density << "\n";
			}
		}
	}

	// 写入向量数据：速度
	os << "\nVECTORS Velocity double\n";
	for (int k = ks; k <= ke; ++k)
	{
		for (int j = js; j <= je; ++j)
		{
			for (int i = is; i <= ie; ++i)
			{
				int idx = idx_proxy->GetIdx(i, j, k);
				double vx = data_manager->GetPrim(ID_VELOCITY_X, idx);
				double vy = data_manager->GetPrim(ID_VELOCITY_Y, idx);
				double vz = data_manager->GetPrim(ID_VELOCITY_Z, idx);
				os << vx << " " << vy << " " << vz << "\n";
			}
		}
	}

	// 写入标量数据：压力
	os << "\nSCALARS Pressure double 1\n";
	os << "LOOKUP_TABLE default\n";
	for (int k = ks; k <= ke; ++k)
	{
		for (int j = js; j <= je; ++j)
		{
			for (int i = is; i <= ie; ++i)
			{
				int idx = idx_proxy->GetIdx(i, j, k);
				double pressure = data_manager->GetPrim(ID_PRESSURE, idx);
				os << pressure << "\n";
			}
		}
	}
}

void Visual::WriteVtkASCII(shared_ptr<NSFieldFNFDM> field, std::ostream& os)
{
	auto grid = field->GetGrid();
	NodeFN& node = grid->GetNode();
	CellFN& cell = grid->GetCell();
	auto data_manager = field->GetDataManager();
	auto node_num = grid->GetTotalNodeNum();
	const double* density = data_manager->GetPrim(ID_DENSITY);
	const double* velocity_x = data_manager->GetPrim(ID_VELOCITY_X);
	const double* velocity_y = data_manager->GetPrim(ID_VELOCITY_Y);
	const double* velocity_z = data_manager->GetPrim(ID_VELOCITY_Z);
	const double* pressure = data_manager->GetPrim(ID_PRESSURE);
	os << "# vtk DataFile Version 3.0\n";
	os << "Flow Field Data\n";
	os << "ASCII\n";
	os << "DATASET UNSTRUCTURED_GRID\n";
	os << "POINTS " << node_num << " double\n";
	for (int iNode = 0; iNode < node_num; ++iNode)
	{
		auto coord = node.GetCoord(iNode);
		os << coord[0] << " " << coord[1] << " " << coord[2] << "\n";
	}
	int cell_num = cell.GetCellNum();
	int cell_node_num = 8;
	int cell_count = cell_num * cell_node_num;
	os << "CELLS " << cell_num << " " << cell_count << "\n";
	for (int iCell = 0; iCell < cell_num; ++iCell)
	{
		auto cell2node = cell.GetNode(iCell);
		os << cell_node_num << " ";
		for (int i = 0; i < cell_node_num; ++i)
		{
			os << cell2node[i] << " ";
		}
		os << "\n";
	}
	os << "CELL_TYPES " << 12 << "\n";
	for (int iCell = 0; iCell < cell_num; ++iCell)
	{
		os << 12 << "\n";
	}
	os << "POINT_DATA " << node_num << "\n";
	os << "SCALARS Density double 1\n";
	os << "LOOKUP_TABLE default\n";
	for (int iNode = 0; iNode < node_num; ++iNode)
	{
		os << density[iNode] << "\n";
	}
	os << "VECTORS Velocity double\n";
	for (int iNode = 0; iNode < node_num; ++iNode)
	{
		os << velocity_x[iNode] << " " << velocity_y[iNode] << " " << velocity_z[iNode] << "\n";
	}
	os << "SCALARS Pressure double 1\n";
	os << "LOOKUP_TABLE default\n";
	for (int iNode = 0; iNode < node_num; ++iNode)
	{
		os << pressure[iNode] << "\n";
	}

}

void Visual::WriteVtkASCII(shared_ptr<NSFieldStruct> field, std::ostream& os)
{
	auto grid = field->GetGrid();
	auto node = grid->GetNode();
	auto data_manager = field->GetDataManager();
	auto idx_proxy = grid->GetIdxProxy();

	int is, ie, js, je, ks, ke;
	grid->GetRange(is, ie, js, je, ks, ke);
	int ni = ie - is + 1;
	int nj = je - js + 1;
	int nk = ke - ks + 1;
	int node_num = ni * nj * nk;

	// 写入 VTK 文件头
	os << "# vtk DataFile Version 3.0\n";
	os << "Flow Field Data\n";
	os << "ASCII\n";
	os << "DATASET STRUCTURED_GRID\n";
	os << "DIMENSIONS " << ni << " " << nj << " " << nk << "\n";
	os << "POINTS " << node_num << " double\n";

	// 写入点坐标
	for (int k = ks; k <= ke; ++k)
	{
		for (int j = js; j <= je; ++j)
		{
			for (int i = is; i <= ie; ++i)
			{
				auto coord = node->GetCoord(i, j, k);
				os << coord[0] << " " << coord[1] << " " << coord[2] << "\n";
			}
		}
	}

	// 写入点数据
	os << "\nPOINT_DATA " << node_num << "\n";

	// 写入标量数据：密度
	os << "SCALARS Density double 1\n";
	os << "LOOKUP_TABLE default\n";
	for (int k = ks; k <= ke; ++k)
	{
		for (int j = js; j <= je; ++j)
		{
			for (int i = is; i <= ie; ++i)
			{
				int idx = idx_proxy->GetIdx(i, j, k);
				double density = data_manager->GetPrim(ID_DENSITY, idx);
				os << density << "\n";
			}
		}
	}

	// 写入向量数据：速度
	os << "\nVECTORS Velocity double\n";
	for (int k = ks; k <= ke; ++k)
	{
		for (int j = js; j <= je; ++j)
		{
			for (int i = is; i <= ie; ++i)
			{
				int idx = idx_proxy->GetIdx(i, j, k);
				double vx = data_manager->GetPrim(ID_VELOCITY_X, idx);
				double vy = data_manager->GetPrim(ID_VELOCITY_Y, idx);
				double vz = data_manager->GetPrim(ID_VELOCITY_Z, idx);
				os << vx << " " << vy << " " << vz << "\n";
			}
		}
	}

	// 写入标量数据：压力
	os << "\nSCALARS Pressure double 1\n";
	os << "LOOKUP_TABLE default\n";
	for (int k = ks; k <= ke; ++k)
	{
		for (int j = js; j <= je; ++j)
		{
			for (int i = is; i <= ie; ++i)
			{
				int idx = idx_proxy->GetIdx(i, j, k);
				double pressure = data_manager->GetPrim(ID_PRESSURE, idx);
				os << pressure << "\n";
			}
		}
	}
}

void Visual::WriteVtkASCII(shared_ptr<FieldManager> field_manager)
{
	std::string work_dir = GlobalData::GetString("work_dir");
	std::string file_name = "result\\" + std::to_string(GlobalData::GetInt("currentIter")) + ".vtu";
	file_name = work_dir + "\\" + file_name;
	std::ofstream out(file_name);
	for (size_t iter_field = 0; iter_field < field_manager->GetFieldNum(); iter_field++)
	{
		auto field = field_manager->GetField(iter_field);
		auto field_type = field->GetFieldType();
		if (field_type == FieldType::NS_Structured)
		{
			auto field_struct = std::dynamic_pointer_cast<NSFieldStruct>(field);
			WriteVtkASCII(field_struct, out);
		}
		else if (field_type == FieldType::NS_Zaran)
		{
			auto field_zaran = std::dynamic_pointer_cast<NSFieldZaran>(field);
			WriteVtkASCII(field_zaran, out);
		}
		else if (field_type == FieldType::NS_FlexibleNode)
		{
			auto field_fn = std::dynamic_pointer_cast<NSFieldFNFDM>(field);
			WriteVtkASCII(field_fn, out);
		}
		else
		{
			Log::warn("Field type is not supported!");
		}
	}
}

void Visual::WriteVtkBinary(shared_ptr<FieldManager> field_manager)
{
	std::string work_dir = GlobalData::GetString("work_dir");
	std::string file_name = "result\\" + std::to_string(GlobalData::GetInt("currentIter")) + ".vtk";
	file_name = work_dir + "\\" + file_name;
	std::ofstream out(file_name);
	for (size_t iter_field = 0; iter_field < field_manager->GetFieldNum(); iter_field++)
	{
		auto field = field_manager->GetField(iter_field);
		auto field_type = field->GetFieldType();
		if (field_type == FieldType::NS_Structured)
		{
			auto field_struct = std::dynamic_pointer_cast<NSFieldStruct>(field);
			WriteVtkBinary(field_struct, out);
		}
		else
		{
			Log::warn("Field type is not supported!");
		}
	}
}

void Visual::WriteVtkBinary(shared_ptr<NSFieldStruct> field, std::ostream& os)
{
	auto grid = field->GetGrid();
	auto node = grid->GetNode();
	auto data_manager = field->GetDataManager();
	auto idx_proxy = grid->GetIdxProxy();

	int is, ie, js, je, ks, ke;
	grid->GetRange(is, ie, js, je, ks, ke);
	int ni = ie - is + 1;
	int nj = je - js + 1;
	int nk = ke - ks + 1;
	int node_num = ni * nj * nk;

	auto SwapEndian = [](double& value)
		{
			char* p = reinterpret_cast<char*>(&value);
			std::reverse(p, p + sizeof(double));
		};


	// 写入 VTK 文件头部信息
	std::string header =
		"# vtk DataFile Version 3.0\n"
		"Flow Field Data\n"
		"BINARY\n"
		"DATASET STRUCTURED_GRID\n"
		"DIMENSIONS " + std::to_string(ni) + " " + std::to_string(nj) + " " + std::to_string(nk) + "\n"
		"POINTS " + std::to_string(node_num) + " double\n";
	os.write(header.c_str(), header.size());

	// 写入节点坐标数据
	for (int k = ks; k <= ke; ++k)
	{
		for (int j = js; j <= je; ++j)
		{
			for (int i = is; i <= ie; ++i)
			{
				double coord[3];
				auto node_coord = node->GetCoord(i, j, k);
				coord[0] = node_coord[0];
				coord[1] = node_coord[1];
				coord[2] = node_coord[2];

				// 转换为大端字节序
				SwapEndian(coord[0]);
				SwapEndian(coord[1]);
				SwapEndian(coord[2]);

				// 写入坐标数据
				os.write(reinterpret_cast<char*>(coord), sizeof(coord));
			}
		}
	}
}

void zaran::Visual::WriteTecplotASCII(shared_ptr<NSFieldStruct> field, std::ostream& os)
{
	auto grid = field->GetGrid();
	auto node = grid->GetNode();
	auto data_manager = field->GetDataManager();
	auto solver = field->GetSolver();
	auto metrics = solver->GetNodeMetrics();
	int grid_ni = grid->GetNi();
	int grid_nj = grid->GetNj();
	int grid_nk = grid->GetNk();
	IdProxyStruct* idx_proxy = new IdProxyStruct(grid_ni, grid_nj, grid_nk);
	int is, ie, js, je, ks, ke;
	grid->GetRange(is, ie, js, je, ks, ke);
	int ni = ie - is + 1;
	int nj = je - js + 1;
	int nk = ke - ks + 1;
	int node_num = ni * nj * nk;
	const double* density = data_manager->GetPrim(ID_DENSITY);
	const double* velocity_x = data_manager->GetPrim(ID_VELOCITY_X);
	const double* velocity_y = data_manager->GetPrim(ID_VELOCITY_Y);
	const double* velocity_z = data_manager->GetPrim(ID_VELOCITY_Z);
	const double* pressure = data_manager->GetPrim(ID_PRESSURE);
	//double zeta_z_Linf = -LARGE_NUMBER;
	//double zeta_z_L2 = 0.0;
	//double zeta_z_error;
	//double zeta_z_exact = 100;
	//for (int k = ks; k <= ke; ++k)
	//{
	//	for (int j = js; j <= je; ++j)
	//	{
	//		for (int i = is; i <= ie; ++i)
	//		{
	//			int idx = idx_proxy->GetIdx(i, j, k);
	//			zeta_z_error = metrics->GetZeta(idx)[2] * metrics->GetJacobian(idx) - zeta_z_exact;
	//			zeta_z_error /= zeta_z_exact;
	//			zeta_z_error = std::abs(zeta_z_error);
	//			zeta_z_Linf = std::max(zeta_z_Linf, zeta_z_error);
	//			zeta_z_L2 += zeta_z_error* zeta_z_error;
	//		}
	//	}
	//}
	//zeta_z_L2 /= node_num;
	//zeta_z_L2 = std::sqrt(zeta_z_L2);

	os << "TITLE=\"NSFieldStruct Field\"\n";
	os << "VARIABLES=\"X\",\"Y\",\"Z\",\"Density\",\"Velocity_x\",\"Velocity_y\",\"Velocity_z\",\"Pressure\"\n";
	//os << "VARIABLES=\"X\",\"Y\",\"Z\",\"Density\",\"Velocity_x\",\"Velocity_y\",\"Velocity_z\",\"Pressure,\"zeta_z\n";
	os << "ZONE I=" << ni << ", J=" << nj << ", K=" << nk << ", F=POINT\n";
	for (int k = ks; k <= ke; ++k)
	{
		for (int j = js; j <= je; ++j)
		{
			for (int i = is; i <= ie; ++i)
			{
				int idx = idx_proxy->GetIdx(i, j, k);
				//double x = node->GetCoord(i, j, k)[0];
				//double y = node->GetCoord(i, j, k)[1];
				//double rho = 1.0 + 0.01 * y;
				os << node->GetCoord(i, j, k)[0] << " " << node->GetCoord(i, j, k)[1] << " " << node->GetCoord(i, j, k)[2] << " "
					<< density[idx] << " " << velocity_x[idx] << " " << velocity_y[idx] << " " << velocity_z[idx] << " "
					//<< pressure[idx] << " " << (metrics->GetZeta(idx)[2] * metrics->GetJacobian(idx)-zeta_z_exac)/zeta_z_exact << "\n";
					<< pressure[idx] << "\n";
			}
		}
	}
	delete idx_proxy;
	//Log::info("zeta_z L_inf norm:{:2E}, L_2 norm:{:2E}", zeta_z_Linf, zeta_z_L2);
}

void Visual::WriteTecplotASCII(shared_ptr<NSFieldFNFDM> field, std::ostream& os)
{
	auto grid = field->GetGrid();
	NodeFN& node = grid->GetNode();
	CellFN& cell = grid->GetCell();
	FaceFN& face = grid->GetFace();
	auto data_manager = field->GetDataManager();
	auto solver = field->GetSolver();
	auto node_num = grid->GetTotalNodeNum();
	const double* density = data_manager->GetPrim(ID_DENSITY);
	const double* velocity_x = data_manager->GetPrim(ID_VELOCITY_X);
	const double* velocity_y = data_manager->GetPrim(ID_VELOCITY_Y);
	const double* velocity_z = data_manager->GetPrim(ID_VELOCITY_Z);
	const double* pressure = data_manager->GetPrim(ID_PRESSURE);
	os << "VARIABLES=\"X\",\"Y\",\"Z\",\"Density\",\"Velocity_x\",\"Velocity_y\",\"Velocity_z\",\"Pressure\",\"iBlank\"\n";
	//非结构网格
	os << "ZONE T=\"NSFieldFNFDM Field\", N=" << node_num << ", E=" << cell.GetCellNum() << ", F=FEPOINT, ET=BRICK\n";
	for (int iNode = 0; iNode < node_num; ++iNode)
	{
		auto coord = node.GetCoord(iNode);
		os << coord[0] << " " << coord[1] << " " << coord[2] << " " << density[iNode] << " " << velocity_x[iNode] << " " << velocity_y[iNode] << " " << velocity_z[iNode] << " " << pressure[iNode] << " " << 0 << "\n";
	}
	for (index_type iCell = 0; iCell < cell.GetCellNum(); ++iCell)
	{
		auto cell2node = cell.GetNode(iCell);
		for (int iNode = 0; iNode < cell.GetNodeNum(iCell); ++iNode)
		{
			os << cell2node[iNode] + 1 << " ";
		}
		if (node_num < 8)
		{
			for (int i = node_num; i < 8; ++i)
			{
				os << cell2node[node_num] + 1 << " ";
			}
		}
		os << "\n";
	}
	os << "ZONE T=\"Bound\", N=" << node_num << ", E=" << face.GetFaceNum() << ", F=FEPOINT, ET=QUADRILATERAL\n";
	for (int iNode = 0; iNode < node_num; ++iNode)
	{
		auto coord = node.GetCoord(iNode);
		os << coord[0] << " " << coord[1] << " " << coord[2] << " " << density[iNode] << " " << velocity_x[iNode] << " " << velocity_y[iNode] << " " << velocity_z[iNode] << " " << pressure[iNode] << "  " << 0 << "\n";
	}
	for (int iFace = 0; iFace < face.GetFaceNum(); ++iFace)
	{
		auto face2node = face.GetFace2Node(iFace);
		auto n_node = face.GetFaceNodeNum(iFace);
		for (int iNode = 0; iNode < n_node; ++iNode)
		{
			os << face2node[iNode] + 1 << " ";
		}
		if (n_node < 4)
		{
			for (int i = n_node; i < 4; ++i)
			{
				os << face2node[0] + 1 << " ";
			}
		}
		os << "\n";
	}




}

void zaran::Visual::WriteTecplotASCII(shared_ptr<NSFieldZaran> field, std::ostream& os)
{
	auto grid = field->GetGrid();
	auto node = grid->GetNode();
	auto data_manager = field->GetDataManager();
	auto solver = field->GetSolver();
	int grid_ni = grid->GetNi();
	int grid_nj = grid->GetNj();
	int grid_nk = grid->GetNk();
	auto idx_proxy = grid->GetIdxProxy();
	int is, ie, js, je, ks, ke;
	grid->GetRange(is, ie, js, je, ks, ke);
	int ni = ie - is + 1;
	int nj = je - js + 1;
	int nk = ke - ks + 1;
	int node_num = ni * nj * nk;
	const double* density = data_manager->GetPrim(ID_DENSITY);
	const double* velocity_x = data_manager->GetPrim(ID_VELOCITY_X);
	const double* velocity_y = data_manager->GetPrim(ID_VELOCITY_Y);
	const double* velocity_z = data_manager->GetPrim(ID_VELOCITY_Z);
	const double* pressure = data_manager->GetPrim(ID_PRESSURE);
	// 以ASCII格式写入，后期可以改为二进制格式
	os << "TITLE=\"Flow Field\"\n";
	os << "VARIABLES=\"X\",\"Y\",\"Z\",\"Density\",\"Velocity_x\",\"Velocity_y\",\"Velocity_z\",\"Pressure\",\"iBlank\"\n";
	os << "ZONE T=\"block grid\", I=" << ni << ", J=" << nj << ", K=" << nk << ", F=POINT\n";
	for (int k = ks; k <= ke; ++k)
	{
		for (int j = js; j <= je; ++j)
		{
			for (int i = is; i <= ie; ++i)
			{
				int idx = idx_proxy->GetIdx(i, j, k);
				os << node->GetCoord(i, j, k)[0] << " " << node->GetCoord(i, j, k)[1] << " " << node->GetCoord(i, j, k)[2] << " " << density[idx] << " "
					<< velocity_x[idx] << " " << velocity_y[idx] << " " << velocity_z[idx] << " " << pressure[idx] << "  " << (int)grid->GetIBlank(i, j, k) << "\n";
			}
		}
	}
}
void zaran::Visual::WriteTecplotBinary(shared_ptr<NSFieldZaran> field)
{
	auto data_manager = field->GetDataManager();
	auto grid = field->GetGrid();
	auto node = grid->GetNode();
	int is, ie, js, je, ks, ke;
	grid->GetRange(is, ie, js, je, ks, ke);
	int ni = ie - is + 1;
	int nj = je - js + 1;
	int nk = ke - ks + 1;
	int grid_ni = grid->GetNi();
	int grid_nj = grid->GetNj();
	int grid_nk = grid->GetNk();
	auto idx_proxy = new IdProxyStruct(grid_ni, grid_nj, grid_nk);
	INTEGER4 node_num = ni * nj * nk;
	INTEGER4 cell_num = (ni - 1) * (nj - 1) * (nk - 1);
	dynamic_array<double> x(node_num), y(node_num), z(node_num), density(node_num),
		velocity_x(node_num), velocity_y(node_num), velocity_z(node_num),
		pressure(node_num);
	for (index_type k = 0; k < nk; ++k)
	{
		for (index_type j = 0; j < nj; ++j)
		{
			for (index_type i = 0; i < ni; ++i)
			{
				index_type idx = i + ni * j + ni * nj * k;
				x[idx] = node->GetCoord(i + is, j + js, k + ks)[0];
				y[idx] = node->GetCoord(i + is, j + js, k + ks)[1];
				z[idx] = node->GetCoord(i + is, j + js, k + ks)[2];
				index_type idx0 = idx_proxy->GetIdx(i + is, j + js, k + ks);
				density[idx] = data_manager->GetPrim(ID_DENSITY, idx0);
				velocity_x[idx] = data_manager->GetPrim(ID_VELOCITY_X, idx0);
				velocity_y[idx] = data_manager->GetPrim(ID_VELOCITY_Y, idx0);
				velocity_z[idx] = data_manager->GetPrim(ID_VELOCITY_Z, idx0);
				pressure[idx] = data_manager->GetPrim(ID_PRESSURE, idx0);
			}
		}
	}
	INTEGER4 vIsDouble = 1;
	string zone_name = grid->GetName() + std::to_string(field->GetIdx());
	INTEGER4 zone_type = 0; // Brick
	INTEGER4 face_num = 6;
	INTEGER4 iCellMax = 0;
	INTEGER4 jCellMax = 0;
	INTEGER4 kCellMax = 0;
	double solution_time = GlobalData::GetDouble("currentTime");
	INTEGER4 strandID = 1;
	INTEGER4 parentZn = 0;
	INTEGER4 isBlock = 1;
	INTEGER4 TotalNumFaceNodes = 1;
	INTEGER4 TotalNumBndryFaces = 1;
	INTEGER4 TotalNumBndryConnections = 1;
	INTEGER4 nFConns = 0;
	INTEGER4 FNMode = 0;
	int shrConn = 0;
	int  i = TECZNE142((char*)zone_name.c_str(),
		&zone_type,
		&ni,
		&nj,
		&nk,
		&iCellMax,
		&jCellMax,
		&kCellMax,
		&solution_time,
		&strandID,
		&parentZn,
		&isBlock,
		&nFConns,
		&FNMode,
		&TotalNumFaceNodes,
		&TotalNumBndryFaces,
		&TotalNumBndryConnections,
		NULL,
		NULL,
		NULL,
		&shrConn);

	i = TECDAT142(&node_num, x.data(), &vIsDouble);
	i = TECDAT142(&node_num, y.data(), &vIsDouble);
	i = TECDAT142(&node_num, z.data(), &vIsDouble);
	i = TECDAT142(&node_num, density.data(), &vIsDouble);
	i = TECDAT142(&node_num, velocity_x.data(), &vIsDouble);
	i = TECDAT142(&node_num, velocity_y.data(), &vIsDouble);
	i = TECDAT142(&node_num, velocity_z.data(), &vIsDouble);
	i = TECDAT142(&node_num, pressure.data(), &vIsDouble);
	delete idx_proxy;
}
void zaran::Visual::WriteTecplotBinary(shared_ptr<NSFieldStruct> field)
{
	auto data_manager = field->GetDataManager();
	auto grid = field->GetGrid();
	auto node = grid->GetNode();
	int is, ie, js, je, ks, ke;
	grid->GetRange(is, ie, js, je, ks, ke);
	int ni = ie - is + 1;
	int nj = je - js + 1;
	int nk = ke - ks + 1;
	int grid_ni = grid->GetNi();
	int grid_nj = grid->GetNj();
	int grid_nk = grid->GetNk();
	auto idx_proxy = new IdProxyStruct(grid_ni, grid_nj, grid_nk);
	INTEGER4 node_num = ni * nj * nk;
	INTEGER4 cell_num = (ni - 1) * (nj - 1) * (nk - 1);
	dynamic_array<double> x(node_num), y(node_num), z(node_num), density(node_num),
		velocity_x(node_num), velocity_y(node_num), velocity_z(node_num),
		pressure(node_num), density_error(node_num);
	for (int k = 0; k < nk; ++k)
	{
		for (int j = 0; j < nj; ++j)
		{
			for (int i = 0; i < ni; ++i)
			{
				int idx = i + ni * j + ni * nj * k;
				x[idx] = node->GetCoord(i + is, j + js, k + ks)[0];
				y[idx] = node->GetCoord(i + is, j + js, k + ks)[1];
				z[idx] = node->GetCoord(i + is, j + js, k + ks)[2];
				int idx0 = idx_proxy->GetIdx(i + is, j + js, k + ks);
				density[idx] = data_manager->GetPrim(ID_DENSITY, idx0);
				velocity_x[idx] = data_manager->GetPrim(ID_VELOCITY_X, idx0);
				velocity_y[idx] = data_manager->GetPrim(ID_VELOCITY_Y, idx0);
				velocity_z[idx] = data_manager->GetPrim(ID_VELOCITY_Z, idx0);
				pressure[idx] = data_manager->GetPrim(ID_PRESSURE, idx0);
				double gamma = 1.4;
				double beta = 5.0;
				double r2 = x[idx] * x[idx] + y[idx] * y[idx];
				double density_exact = pow(1.0 - (gamma - 1.0) * beta * beta * exp(1.0 - r2) / (8.0 * gamma * PI * PI), 1.0 / (gamma - 1.0));
				density_error[idx] = (density[idx] - density_exact) / density_exact;
			}
		}
	}
	INTEGER4 vIsDouble = 1;
	string zone_name = grid->GetName() + std::to_string(field->GetIdx());
	INTEGER4 zone_type = 0; // Brick
	INTEGER4 face_num = 6;
	INTEGER4 iCellMax = 0;
	INTEGER4 jCellMax = 0;
	INTEGER4 kCellMax = 0;
	double solution_time = GlobalData::GetDouble("currentTime");
	INTEGER4 strandID = 0;
	INTEGER4 parentZn = 0;
	INTEGER4 isBlock = 1;
	INTEGER4 TotalNumFaceNodes = 1;
	INTEGER4 TotalNumBndryFaces = 1;
	INTEGER4 TotalNumBndryConnections = 1;
	INTEGER4 nFConns = 0;
	INTEGER4 FNMode = 0;
	int shrConn = 0;
	int  i = TECZNE142((char*)zone_name.c_str(),
		&zone_type,
		&ni,
		&nj,
		&nk,
		&iCellMax,
		&jCellMax,
		&kCellMax,
		&solution_time,
		&strandID,
		&parentZn,
		&isBlock,
		&nFConns,
		&FNMode,
		&TotalNumFaceNodes,
		&TotalNumBndryFaces,
		&TotalNumBndryConnections,
		NULL,
		NULL,
		NULL,
		&shrConn);

	i = TECDAT142(&node_num, x.data(), &vIsDouble);
	i = TECDAT142(&node_num, y.data(), &vIsDouble);
	i = TECDAT142(&node_num, z.data(), &vIsDouble);
	i = TECDAT142(&node_num, density.data(), &vIsDouble);
	i = TECDAT142(&node_num, velocity_x.data(), &vIsDouble);
	i = TECDAT142(&node_num, velocity_y.data(), &vIsDouble);
	i = TECDAT142(&node_num, velocity_z.data(), &vIsDouble);
	i = TECDAT142(&node_num, pressure.data(), &vIsDouble);
	delete idx_proxy;
}
void Visual::WriteTecASCII(shared_ptr<FieldManager> field_manager)
{
	std::string work_dir = GlobalData::GetString("work_dir");
	std::string file_name = "result\\" + std::to_string(GlobalData::GetInt("currentIter")) + ".dat";
	file_name = work_dir + "\\" + file_name;
	std::ofstream out(file_name);
	//out << "TITLE=\"Flow Field\"\n";
	//out << "VARIABLES=\"X\",\"Y\",\"Z\",\"Density\",\"Velocity_x\",\"Velocity_y\",\"Velocity_z\",\"Pressure\",\"iBlank\"\n";
	double solution_time = GlobalData::GetDouble("currentTime");
	out << "SOLUTIONTIME=" << solution_time << "\n";
	for (size_t iter_field = 0; iter_field < field_manager->GetFieldNum(); iter_field++)
	{
		auto field = field_manager->GetField(iter_field);
		auto field_type = field->GetFieldType();
		if (field_type == FieldType::NS_Structured)
		{
			auto field_struct = std::dynamic_pointer_cast<NSFieldStruct>(field);
			WriteTecplotASCII(field_struct, out);
		}
		else if (field_type == FieldType::NS_Zaran)
		{
			auto field_zaran = std::dynamic_pointer_cast<NSFieldZaran>(field);
			WriteTecplotASCII(field_zaran, out);
		}
		else if (field_type == FieldType::NS_FlexibleNode)
		{
			auto field_fn = std::dynamic_pointer_cast<NSFieldFNFDM>(field);
			WriteTecplotASCII(field_fn, out);
		}
		else
		{
			Log::warn("Field type is not supported!");
		}
	}

}

void zaran::Visual::WriteTecplotBinary(shared_ptr<FieldManager> field_manager)
{
	INTEGER4 file_format = 0;
	INTEGER4 debug = 0;
	INTEGER4 vIsDouble = 1;
	INTEGER4 vIsInt = 0;
	INTEGER4 fileType = 0;
	string grid_name = "grid";
	string var_name = "x, y, z, density, velocity_x, velocity_y, velocity_z, pressure";
	std::string work_dir = GlobalData::GetString("work_dir");
	std::string file_name = "result\\" + std::to_string(GlobalData::GetInt("currentIter")) + ".plt";
	file_name = work_dir + "\\" + file_name;
	Log::info("Write Tecplot file: {}", file_name);
	int i = TECINI142(grid_name.c_str(), var_name.c_str(), file_name.c_str(), (char*)".", &file_format, &fileType,
		&debug, &vIsDouble);
	for (size_t iter_field = 0; iter_field < field_manager->GetFieldNum(); iter_field++)
	{
		auto field = field_manager->GetField(iter_field);
		auto field_type = field->GetFieldType();
		if (field_type == FieldType::NS_Zaran)
		{
			auto field_zaran = std::dynamic_pointer_cast<NSFieldZaran>(field);
			WriteTecplotBinary(field_zaran);
		}
		else if (field_type == FieldType::NS_FlexibleNode)
		{
			auto field_fn = std::dynamic_pointer_cast<NSFieldFNFDM>(field);
			WriteTecplotBinary(field_fn);
		}
		else if (field_type == FieldType::NS_Structured)
		{
			auto field_struct = std::dynamic_pointer_cast<NSFieldStruct>(field);
			WriteTecplotBinary(field_struct);
		}
		else
		{
			Log::warn("Field type is not supported!");
		}
	}
	i = TECEND142();
}
