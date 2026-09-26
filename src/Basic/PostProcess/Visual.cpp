#include "Visual.h"
#include "GlobalData.h"
#include "Log.h"
#include "NSFieldFN.h"
#include "DEMField.h"
#include "FastNumberFormat.h"
#include "ZaranError.h"
#include <TECIO.h>
#include <cgnslib.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <string>
#include <utility>
#include <vector>
using namespace zaran;
void zaran::Visual::WriteTecplotBinary(shared_ptr<NSFieldFNFDM> field)
{
    auto data_manager = field->GetDataManager();
    auto grid = field->GetGrid();
    CellFN &cell = grid->GetCell();
    NodeFN &node = grid->GetNode();

    const double *density, *velocity_x, *velocity_y, *velocity_z, *pressure;
    density = data_manager->GetPrim(ID_DENSITY);
    velocity_x = data_manager->GetPrim(ID_VELOCITY_X);
    velocity_y = data_manager->GetPrim(ID_VELOCITY_Y);
    velocity_z = data_manager->GetPrim(ID_VELOCITY_Z);
    pressure = data_manager->GetPrim(ID_PRESSURE);
    INTEGER4 node_num = grid->GetTotalNodeNum();
    INTEGER4 cell_num = cell.GetCellNum();
    dynamic_array<double> x(node_num), y(node_num), z(node_num), iblank(node_num);
    for (int iNode = 0; iNode < node_num; ++iNode)
    {
        x[iNode] = node.GetCoord(iNode)[0];
        y[iNode] = node.GetCoord(iNode)[1];
        z[iNode] = node.GetCoord(iNode)[2];
        iblank[iNode] = -1;//不显示
    }
    INTEGER4 vIsDouble = 1;
    string zone_name = "grid_" + std::to_string(field->GetIdx());
    INTEGER4 zone_type = 5; // Brick
    INTEGER4 face_num = 6;
    INTEGER4 iCellMax = 0;
    INTEGER4 jCellMax = 0;
    INTEGER4 kCellMax = 0;
    double solution_time = GlobalData::GetDouble("iteration.current_time");
    INTEGER4 strandID = 2;
    INTEGER4 parentZn = 0;
    INTEGER4 isBlock = 1;
    INTEGER4 nFConns = 0;
    INTEGER4 FNMode = 0;
    int valueLocation[] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    int shrConn = 0;
    int i = TECZNE142((char *)zone_name.c_str(), &zone_type, &node_num, &cell_num, &face_num, &iCellMax, &jCellMax,
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
    i = TECDAT142(&node_num, iblank.data(), &vIsDouble);
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
    strandID = 3;
    parentZn = 0;
    isBlock = 1;
    nFConns = 0;
    FNMode = 0;
    i = TECZNE142((char *)zone_name.c_str(), &zone_type, &node_num, &cell_num, &face_num, &iCellMax, &jCellMax,
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
    i = TECDAT142(&node_num, iblank.data(), &vIsDouble);
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
    std::string file_name = "result\\" + std::to_string(GlobalData::GetInt("iteration.current_iter")) + "_bound" + ".cgns";
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

    const double *density, *velocity_x, *velocity_y, *velocity_z, *pressure;
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
    dynamic_array<cgsize_t> elements(cell_num * 8);
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
    cg_section_write(index_file, index_base, index_zone, "Section", CGNS_ENUMV(HEXA_8), 1, cell_num, 0, elements.data(),
                     &index_section);
    int index_field;
    cg_sol_write(index_file, index_base, index_zone, "FlowSolution", CGNS_ENUMV(Vertex), &index_field);
    cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "Density", density,
                   &index_field);
    cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "VelocityX", velocity_x,
                   &index_field);
    cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "VelocityY", velocity_y,
                   &index_field);
    cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "VelocityZ", velocity_z,
                   &index_field);
    cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "Pressure", pressure,
                   &index_field);
    // 输出边界面
    int index_file2, index_base2;
    file_name = "result\\" + std::to_string(GlobalData::GetInt("iteration.current_iter")) + "_bound_face" + ".cgns";
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
    cg_coord_write(index_file2, index_base2, index_zone1, CGNS_ENUMV(RealDouble), "CoordinateX", x.data(),
                   &index_coord1);
    cg_coord_write(index_file2, index_base2, index_zone1, CGNS_ENUMV(RealDouble), "CoordinateY", y.data(),
                   &index_coord1);
    cg_coord_write(index_file2, index_base2, index_zone1, CGNS_ENUMV(RealDouble), "CoordinateZ", z.data(),
                   &index_coord1);
    int index_section1;
    cg_section_write(index_file2, index_base2, index_zone1, "bound_Section", CGNS_ENUMV(QUAD_4), 1, face_num, 0,
                     faces.data(), &index_section1);
    int index_field1;
    cg_sol_write(index_file2, index_base2, index_zone1, "bound_FlowSolution", CGNS_ENUMV(Vertex), &index_field1);
    cg_field_write(index_file2, index_base2, index_zone1, index_field1, CGNS_ENUMV(RealDouble), "Density", density,
                   &index_field1);
    cg_field_write(index_file2, index_base2, index_zone1, index_field1, CGNS_ENUMV(RealDouble), "VelocityX", velocity_x,
                   &index_field1);
    cg_field_write(index_file2, index_base2, index_zone1, index_field1, CGNS_ENUMV(RealDouble), "VelocityY", velocity_y,
                   &index_field1);
    cg_field_write(index_file2, index_base2, index_zone1, index_field1, CGNS_ENUMV(RealDouble), "VelocityZ", velocity_z,
                   &index_field1);
    cg_field_write(index_file2, index_base2, index_zone1, index_field1, CGNS_ENUMV(RealDouble), "Pressure", pressure,
                   &index_field1);
}

void Visual::WriteCGNS(shared_ptr<NSFieldZaran> field, cgsize_t index_file1, cgsize_t index_base1)
{
    std::string work_dir = GlobalData::GetString("work_dir");
    std::string file_name = "result\\" + std::to_string(GlobalData::GetInt("iteration.current_iter")) + "_back" + ".cgns";
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
    IdProxyStruct &idx_proxy = grid->GetIdxProxy();
    index_type is, ie, js, je, ks, ke;
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
    dynamic_array<double> density(node_num), velocity_x(node_num), velocity_y(node_num), velocity_z(node_num),
        pressure(node_num);
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
                int idx0 = idx_proxy(i + is, j + js, k + ks);
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
    cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "Density", density.data(),
                   &index_field);
    cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "VelocityX",
                   velocity_x.data(), &index_field);
    cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "VelocityY",
                   velocity_y.data(), &index_field);
    cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "VelocityZ",
                   velocity_z.data(), &index_field);
    cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "Pressure", pressure.data(),
                   &index_field);
    cg_close(index_file);
}

void Visual::WriteCGNS(shared_ptr<NSFieldStruct> field, cgsize_t index_file1, cgsize_t index_base1)
{
    std::string work_dir = GlobalData::GetString("work_dir");
    std::string file_name = "result\\" + std::to_string(GlobalData::GetInt("iteration.current_iter")) + "_back" + ".cgns";
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
    IdProxyStruct &idx_proxy = grid->GetIdxProxy();
    index_type is, ie, js, je, ks, ke;
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
    dynamic_array<double> density(node_num), velocity_x(node_num), velocity_y(node_num), velocity_z(node_num),
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
                int idx0 = idx_proxy(i + is, j + js, k + ks);
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
    cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "Density", density.data(),
                   &index_field);
    cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "VelocityX",
                   velocity_x.data(), &index_field);
    cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "VelocityY",
                   velocity_y.data(), &index_field);
    cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "VelocityZ",
                   velocity_z.data(), &index_field);
    cg_field_write(index_file, index_base, index_zone, index_field, CGNS_ENUMV(RealDouble), "Pressure", pressure.data(),
                   &index_field);
    cg_close(index_file);
}

void Visual::WriteCGNS(shared_ptr<FieldManager> field_manager)
{
    // std::string work_dir = GlobalData::GetString("work_dir");
    // std::string file_name = "result\\" + std::to_string(GlobalData::GetInt("iteration.current_iter")) + ".cgns";
    // file_name = work_dir + "\\" + file_name;
    int index_file;
    // cg_open(file_name.c_str(), CG_MODE_WRITE, &index_file);
    int index_base;
    // int cell_dim = 3;
    // int phys_dim = 3;
    // cg_base_write(index_file, "Base", cell_dim, phys_dim, &index_base);
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
    // cg_close(index_file);
}

void Visual::WriteVtkASCII(shared_ptr<NSFieldZaran> field, std::ostream &os)
{
    auto grid = field->GetGrid();
    auto node = grid->GetNode();
    auto data_manager = field->GetDataManager();
    IdProxyStruct &idx_proxy = grid->GetIdxProxy();

    index_type is, ie, js, je, ks, ke;
    grid->GetRange(is, ie, js, je, ks, ke);
    index_type ni = ie - is + 1;
    index_type nj = je - js + 1;
    index_type nk = ke - ks + 1;
    index_type node_num = ni * nj * nk;

    // 写入 VTK 文件头
    os << "# vtk DataFile Version 3.0\n";
    os << "Flow Field Data\n";
    os << "ASCII\n";
    os << "DATASET STRUCTURED_GRID\n";
    os << "DIMENSIONS " << ni << " " << nj << " " << nk << "\n";
    os << "POINTS " << node_num << " double\n";

    // 写入点坐标
    for (index_type k = ks; k <= ke; ++k)
    {
        for (index_type j = js; j <= je; ++j)
        {
            for (index_type i = is; i <= ie; ++i)
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
    for (index_type k = ks; k <= ke; ++k)
    {
        for (index_type j = js; j <= je; ++j)
        {
            for (index_type i = is; i <= ie; ++i)
            {
                index_type idx = idx_proxy(i, j, k);
                double density = data_manager->GetPrim(ID_DENSITY, idx);
                os << density << "\n";
            }
        }
    }

    // 写入向量数据：速度
    os << "\nVECTORS Velocity double\n";
    for (index_type k = ks; k <= ke; ++k)
    {
        for (index_type j = js; j <= je; ++j)
        {
            for (index_type i = is; i <= ie; ++i)
            {
                int idx = idx_proxy(i, j, k);
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
    for (index_type k = ks; k <= ke; ++k)
    {
        for (index_type j = js; j <= je; ++j)
        {
            for (index_type i = is; i <= ie; ++i)
            {
                index_type idx = idx_proxy(i, j, k);
                double pressure = data_manager->GetPrim(ID_PRESSURE, idx);
                os << pressure << "\n";
            }
        }
    }
}

void Visual::WriteVtkASCII(shared_ptr<NSFieldFNFDM> field, std::ostream &os)
{
    auto grid = field->GetGrid();
    NodeFN &node = grid->GetNode();
    CellFN &cell = grid->GetCell();
    auto data_manager = field->GetDataManager();
    auto node_num = grid->GetTotalNodeNum();
    const double *density = data_manager->GetPrim(ID_DENSITY);
    const double *velocity_x = data_manager->GetPrim(ID_VELOCITY_X);
    const double *velocity_y = data_manager->GetPrim(ID_VELOCITY_Y);
    const double *velocity_z = data_manager->GetPrim(ID_VELOCITY_Z);
    const double *pressure = data_manager->GetPrim(ID_PRESSURE);
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

void Visual::WriteVtkASCII(shared_ptr<NSFieldStruct> field, std::ostream &os)
{
    auto grid = field->GetGrid();
    auto node = grid->GetNode();
    auto data_manager = field->GetDataManager();
    IdProxyStruct &idx_proxy = grid->GetIdxProxy();

    index_type is, ie, js, je, ks, ke;
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
                int idx = idx_proxy(i, j, k);
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
                int idx = idx_proxy(i, j, k);
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
                int idx = idx_proxy(i, j, k);
                double pressure = data_manager->GetPrim(ID_PRESSURE, idx);
                os << pressure << "\n";
            }
        }
    }
}

void Visual::WriteVtkASCII(shared_ptr<FieldManager> field_manager)
{
    std::string work_dir = GlobalData::GetString("work_dir");
    std::string file_name = "result\\" + std::to_string(GlobalData::GetInt("iteration.current_iter")) + ".vtu";
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
    std::string file_name = "result\\" + std::to_string(GlobalData::GetInt("iteration.current_iter")) + ".vtk";
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

void Visual::WriteVtkBinary(shared_ptr<NSFieldStruct> field, std::ostream &os)
{
    auto grid = field->GetGrid();
    auto node = grid->GetNode();
    auto data_manager = field->GetDataManager();
    IdProxyStruct &idx_proxy = grid->GetIdxProxy();

    index_type is, ie, js, je, ks, ke;
    grid->GetRange(is, ie, js, je, ks, ke);
    index_type ni = ie - is + 1;
    index_type nj = je - js + 1;
    index_type nk = ke - ks + 1;
    index_type node_num = ni * nj * nk;

    auto SwapEndian = [](double &value) {
        char *p = reinterpret_cast<char *>(&value);
        std::reverse(p, p + sizeof(double));
    };

    // 写入 VTK 文件头部信息
    std::string header = "# vtk DataFile Version 3.0\n"
                         "Flow Field Data\n"
                         "BINARY\n"
                         "DATASET STRUCTURED_GRID\n"
                         "DIMENSIONS " +
                         std::to_string(ni) + " " + std::to_string(nj) + " " + std::to_string(nk) +
                         "\n"
                         "POINTS " +
                         std::to_string(node_num) + " double\n";
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
                os.write(reinterpret_cast<char *>(coord), sizeof(coord));
            }
        }
    }
}

bool zaran::Visual::HasVolumeFraction(const shared_ptr<Field>& field)
{
    if (!field)
    {
        return false;
    }
    auto data = field->GetData();
    if (!data)
    {
        return false;
    }
    // 只有两相场（DataManagerNSTwoPhase）才注册了这个量
    return data->HasData("volume_fraction");
}

void zaran::Visual::WriteTecplotASCII(shared_ptr<NSFieldStruct> field, std::ostream &os)
{
    WriteTecplotASCII(field, os, HasVolumeFraction(field));
}

void zaran::Visual::WriteTecplotASCII(shared_ptr<NSFieldStruct> field, std::ostream &os,
    bool with_volume_fraction)
{
    auto grid = field->GetGrid();
    auto node = grid->GetNode();
    auto data_manager = field->GetDataManager();
    auto solver = field->GetSolver();
    auto metrics = solver->GetNodeMetrics();
    index_type grid_ni = grid->GetNi();
    index_type grid_nj = grid->GetNj();
    index_type grid_nk = grid->GetNk();
    IdProxyStruct idx_proxy(grid_ni, grid_nj, grid_nk);
    index_type is, ie, js, je, ks, ke;
    grid->GetRange(is, ie, js, je, ks, ke);
    index_type ni = ie - is + 1;
    index_type nj = je - js + 1;
    index_type nk = ke - ks + 1;
    index_type node_num = ni * nj * nk;
    const double *density = data_manager->GetPrim(ID_DENSITY);
    const double *velocity_x = data_manager->GetPrim(ID_VELOCITY_X);
    const double *velocity_y = data_manager->GetPrim(ID_VELOCITY_Y);
    const double *velocity_z = data_manager->GetPrim(ID_VELOCITY_Z);
    const double *pressure = data_manager->GetPrim(ID_PRESSURE);
    // 两相场：额外输出气相体积分数（该量由 DataManagerNSTwoPhase 注册在 FieldData 里）
    const double *volume_fraction = nullptr;
    if (with_volume_fraction)
    {
        double *volume_fraction_buffer = nullptr;
        field->GetData()->GetData("volume_fraction", volume_fraction_buffer);
        volume_fraction = volume_fraction_buffer;
    }
    // double zeta_z_Linf = -LARGE_NUMBER;
    // double zeta_z_L2 = 0.0;
    // double zeta_z_error;
    // double zeta_z_exact = 100;
    // for (int k = ks; k <= ke; ++k)
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
    // }
    // zeta_z_L2 /= node_num;
    // zeta_z_L2 = std::sqrt(zeta_z_L2);

    // double jacobi_Linf = -LARGE_NUMBER;
    // double jacobi_L2 = 0.0;
    // double jacobi_error;
    // double jacobi_exact = 100;
    // int N = ke - ks;
    // double L = 5;
    // double h = L / N;
    // for (int k = ks+N/5 ; k <= ke- N / 5; ++k)
    //{
    //	for (int j = js + N / 5; j <= je -N / 5; ++j)
    //	{
    //		for (int i = is + N / 5; i <= ie- N / 5; ++i)
    //		{
    //			int idx = idx_proxy(i, j, k);
    //			double Ci = cos(PI * h * i);
    //			double Cj = cos(PI * h * j);
    //			double Ck = cos(PI * h * k);
    //			double Si = sin(PI * h * i);
    //			double Sj = sin(PI * h * j);
    //			double Sk = sin(PI * h * k);
    //			double x_xi = 2 * h;
    //			double x_eta = 2 * PI * h * h * Sk * Cj;
    //			double x_zeta = 2 * PI * h * h * Sj * Ck;
    //			double y_xi = 2 * PI * h * h * Si * Ck;
    //			double y_eta = 2 * h;
    //			double y_zeta = 2 * PI * h * h * Sk * Ci;
    //			double z_xi = 0 * PI * h * h * Sj * Ci;
    //			double z_eta = 0 * PI * h * h * Si * Cj;
    //			double z_zeta = 2 * h;
    //			double jacobi_exact = x_xi * (y_eta * z_zeta - z_eta * y_zeta)
    //								- y_xi * (x_eta * z_zeta - z_eta * x_zeta)
    //								+ z_xi * (x_eta * y_zeta - y_eta * x_zeta);
    //			jacobi_error = 1.0 / metrics->GetJacobian(idx) - jacobi_exact;
    //			jacobi_error /= jacobi_exact;
    //			jacobi_error = std::abs(jacobi_error);
    //			jacobi_Linf = std::max(jacobi_Linf, jacobi_error);
    //			jacobi_L2 += jacobi_error * jacobi_error;
    //		}
    //	}
    // }
    // jacobi_L2 /= pow(N - 2*N / 5, 3);
    // jacobi_L2 = std::sqrt(jacobi_L2);

    // double pressure_inf = -LARGE_NUMBER;
    // double pressure_L2 = 0.0;
    // double pressure_error;
    // double pressure_exact = 1.0 / 1.4;
    // for (int k = ks; k <= ke; ++k)
    //{
    //	for (int j = js; j <= je; ++j)
    //	{
    //		for (int i = is; i <= ie; ++i)
    //		{
    //			int idx = idx_proxy(i, j, k);
    //			pressure_error = pressure[idx] - pressure_exact;
    //			pressure_error /= pressure_exact;
    //			pressure_error = std::abs(pressure_error);
    //			pressure_inf = std::max(pressure_inf, pressure_error);
    //			pressure_L2 += pressure_error * pressure_error;
    //		}
    //	}
    // }
    // pressure_L2 /= node_num;
    // pressure_L2 = std::sqrt(pressure_L2);
    os << "TITLE=\"NSFieldStruct Field\"\n";
    if (volume_fraction != nullptr)
    {
        os << "VARIABLES=\"X\",\"Y\",\"Z\",\"Density\",\"Velocity_x\",\"Velocity_y\",\"Velocity_z\","
              "\"Pressure\",\"Volume_fraction\"\n";
    }
    else
    {
        os << "VARIABLES=\"X\",\"Y\",\"Z\",\"Density\",\"Velocity_x\",\"Velocity_y\",\"Velocity_z\",\"Pressure\"\n";
    }
    double solution_time = GlobalData::GetDouble("iteration.current_time");
    // ZONE 参数全部写在同一行（含 SOLUTIONTIME）：
    // Tecplot 允许控制行续行，但把参数集中在一行是"无歧义"的写法，
    // 也方便脚本按 key=value 直接解析。
    os << "ZONE T=\"Flow\", STRANDID=1, SOLUTIONTIME=" << solution_time
       << ", I=" << ni << ", J=" << nj << ", K=" << nk
       << ", ZONETYPE=ORDERED, DATAPACKING=POINT\n";
    // 该输出通道主要用于脚本后处理与精度验证，因此用能往返 double 的精度
    os.precision(17);

    for (int k = ks; k <= ke; ++k)
    {
        for (int j = js; j <= je; ++j)
        {
            for (int i = is; i <= ie; ++i)
            {

                int idx = idx_proxy(i, j, k);

                // double Ci = cos(PI * h * i);
                // double Cj = cos(PI * h * j);
                // double Ck = cos(PI * h * k);
                // double Si = sin(PI * h * i);
                // double Sj = sin(PI * h * j);
                // double Sk = sin(PI * h * k);
                // double x_xi = 2 * h;
                // double x_eta = 2 * PI * h * h * Sk * Cj;
                // double x_zeta = 2 * PI * h * h * Sj * Ck;
                // double y_xi = 2 * PI * h * h * Si * Ck;
                // double y_eta = 2 * h;
                // double y_zeta = 2 * PI * h * h * Sk * Ci;
                // double z_xi = 0 * PI * h * h * Sj * Ci;
                // double z_eta = 0 * PI * h * h * Si * Cj;
                // double z_zeta = 2 * h;
                // double jacobi_exact = x_xi * (y_eta * z_zeta - z_eta * y_zeta)
                //	- y_xi * (x_eta * z_zeta - z_eta * x_zeta)
                //	+ z_xi * (x_eta * y_zeta - y_eta * x_zeta);

                // jacobi_error = 1.0 / metrics->GetJacobian(idx) - jacobi_exact;
                // jacobi_error /= jacobi_exact;

                // double x = node->GetCoord(i, j, k)[0];
                // double y = node->GetCoord(i, j, k)[1];
                // double rho = 1.0 - 0.005 * x + 0.01 * y;
                // pressure_error = pressure[idx] - pressure_exact;
                os << node->GetCoord(i, j, k)[0] << " " << node->GetCoord(i, j, k)[1] << " "
                   << node->GetCoord(i, j, k)[2] << " " << density[idx] << " " << velocity_x[idx] << " "
                   << velocity_y[idx] << " " << velocity_z[idx]
                   << " "
                   //<< pressure[idx] << " " << (metrics->GetZeta(idx)[2] *
                   // metrics->GetJacobian(idx)-zeta_z_exac)/zeta_z_exact << "\n";
                   //<< pressure[idx] << " " << metrics->GetX(idx)[2] << " " << metrics->GetY(idx)[2] << " " <<
                   // metrics->GetZ(idx)[0] << " " << jacobi_error << "\n";
                   << pressure[idx];
                if (volume_fraction != nullptr)
                {
                    os << " " << volume_fraction[idx];
                }
                os << "\n";
            }
        }
    }
    // Log::info("zeta_z L_inf norm:{:2E}, L_2 norm:{:2E}", zeta_z_Linf, zeta_z_L2);
    // Log::info("presure L_inf norm:{:2E}, L_2 norm:{:2E}", pressure_inf, pressure_L2);
    // Log::info("jacobi L_inf norm:{:2E}, L_2 norm:{:2E}", jacobi_Linf, jacobi_L2);
}

void Visual::WriteTecplotASCII(shared_ptr<NSFieldFNFDM> field, std::ostream &os)
{
    auto grid = field->GetGrid();
    NodeFN &node = grid->GetNode();
    CellFN &cell = grid->GetCell();
    FaceFN &face = grid->GetFace();
    auto data_manager = field->GetDataManager();
    auto solver = field->GetSolver();
    auto node_num = grid->GetTotalNodeNum();
    const double *density = data_manager->GetPrim(ID_DENSITY);
    const double *velocity_x = data_manager->GetPrim(ID_VELOCITY_X);
    const double *velocity_y = data_manager->GetPrim(ID_VELOCITY_Y);
    const double *velocity_z = data_manager->GetPrim(ID_VELOCITY_Z);
    const double *pressure = data_manager->GetPrim(ID_PRESSURE);
    os << "VARIABLES=\"X\",\"Y\",\"Z\",\"Density\",\"Velocity_x\",\"Velocity_y\",\"Velocity_z\",\"Pressure\","
          "\"iBlank\"\n";
    // 非结构网格
    if (grid->GetDim() == 2)
        os << "ZONE T=\"NSFieldFNFDM Field\", N=" << node_num << ", E=" << cell.GetCellNum()
           << ", F=FEPOINT, ET=QUADRILATERAL\n";
    else
        os << "ZONE T=\"NSFieldFNFDM Field\", N=" << node_num << ", E=" << cell.GetCellNum()
           << ", F=FEPOINT, ET=BRICK\n";
    double solution_time = GlobalData::GetDouble("iteration.current_time");
    os << "SOLUTIONTIME=" << solution_time << "\n";
    for (int iNode = 0; iNode < node_num; ++iNode)
    {
        auto coord = node.GetCoord(iNode);
        os << coord[0] << " " << coord[1] << " " << coord[2] << " " << density[iNode] << " " << velocity_x[iNode] << " "
           << velocity_y[iNode] << " " << velocity_z[iNode] << " " << pressure[iNode] << " " << 0 << "\n";
    }
    if (grid->GetDim() == 2)
    {
        for (index_type iCell = 0; iCell < cell.GetCellNum(); ++iCell)
        {
            auto cell2node = cell.GetNode(iCell);
            for (int iNode = 0; iNode < cell.GetNodeNum(iCell); ++iNode)
            {
                os << cell2node[iNode] + 1 << " ";
            }
            if (node_num < 4)
            {
                for (int i = node_num; i < 4; ++i)
                {
                    os << cell2node[node_num] + 1 << " ";
                }
            }
            os << "\n";
        }
    }
    else
    {
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
    }
    if (grid->GetDim() == 2)
        os << "ZONE T=\"Bound\", N=" << node_num << ", E=" << face.GetFaceNum() << ", F=FEPOINT, ET=LINESEG\n";
    else
        os << "ZONE T=\"Bound\", N=" << node_num << ", E=" << face.GetFaceNum() << ", F=FEPOINT, ET=QUADRILATERAL\n";
    for (int iNode = 0; iNode < node_num; ++iNode)
    {
        auto coord = node.GetCoord(iNode);
        os << coord[0] << " " << coord[1] << " " << coord[2] << " " << density[iNode] << " " << velocity_x[iNode] << " "
           << velocity_y[iNode] << " " << velocity_z[iNode] << " " << pressure[iNode] << "  " << 0 << "\n";
    }
    if (grid->GetDim() == 2)
    {
        for (int iFace = 0; iFace < face.GetFaceNum(); ++iFace)
        {
            auto face2node = face.GetFace2Node(iFace);
            auto n_node = face.GetFaceNodeNum(iFace);
            for (int iNode = 0; iNode < n_node; ++iNode)
            {
                os << face2node[iNode] + 1 << " ";
            }
            if (n_node < 2)
            {
                for (int i = n_node; i < 2; ++i)
                {
                    os << face2node[0] + 1 << " ";
                }
            }
            os << "\n";
        }
    }
    else
    {
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
}

void zaran::Visual::WriteTecplotASCII(const shared_ptr<NSFieldZaran> &field, std::ostream &os)
{
    const auto grid = field->GetGrid();
    const auto node = grid->GetNode();
    const auto data_manager = field->GetDataManager();
    auto solver = field->GetSolver();
    index_type grid_ni = grid->GetNi();
    index_type grid_nj = grid->GetNj();
    index_type grid_nk = grid->GetNk();
    IdProxyStruct &idx_proxy = grid->GetIdxProxy();
    index_type is, ie, js, je, ks, ke;
    grid->GetRange(is, ie, js, je, ks, ke);
    index_type ni = ie - is + 1;
    index_type nj = je - js + 1;
    index_type nk = ke - ks + 1;
    index_type node_num = ni * nj * nk;
    const double *density = data_manager->GetPrim(ID_DENSITY);
    const double *velocity_x = data_manager->GetPrim(ID_VELOCITY_X);
    const double *velocity_y = data_manager->GetPrim(ID_VELOCITY_Y);
    const double *velocity_z = data_manager->GetPrim(ID_VELOCITY_Z);
    const double *pressure = data_manager->GetPrim(ID_PRESSURE);
    // 以ASCII格式写入，后期可以改为二进制格式
    os << "TITLE=\"Flow Field\"\n";
    os << "VARIABLES=\"X\",\"Y\",\"Z\",\"Density\",\"Velocity_x\",\"Velocity_y\",\"Velocity_z\",\"Pressure\","
          "\"iBlank\"\n";
    os << "ZONE T=\"block grid\", I=" << ni << ", J=" << nj << ", K=" << nk << ", F=POINT\n";
    double solution_time = GlobalData::GetDouble("iteration.current_time");
    os << "SOLUTIONTIME=" << solution_time << "\n";
    for (index_type k = ks; k <= ke; ++k)
    {
        for (index_type j = js; j <= je; ++j)
        {
            for (index_type i = is; i <= ie; ++i)
            {
                index_type idx = idx_proxy(i, j, k);
                os << node->GetCoord(i, j, k)[0] << " " << node->GetCoord(i, j, k)[1] << " "
                   << node->GetCoord(i, j, k)[2] << " " << density[idx] << " " << velocity_x[idx] << " "
                   << velocity_y[idx] << " " << velocity_z[idx] << " " << pressure[idx] << "  "
                   << (int)grid->GetIBlank(i, j, k) << "\n";
            }
        }
    }
}
void zaran::Visual::WriteTecplotBinary(shared_ptr<NSFieldZaran> field)
{
    auto data_manager = field->GetDataManager();
    auto grid = field->GetGrid();
    auto node = grid->GetNode();
    index_type is, ie, js, je, ks, ke;
    grid->GetRange(is, ie, js, je, ks, ke);
    INTEGER4 ni = ie - is + 1;
    INTEGER4 nj = je - js + 1;
    INTEGER4 nk = ke - ks + 1;
    index_type grid_ni = grid->GetNi();
    index_type grid_nj = grid->GetNj();
    index_type grid_nk = grid->GetNk();
    IdProxyStruct idx_proxy(grid_ni, grid_nj, grid_nk);
    INTEGER4 node_num = ni * nj * nk;
    INTEGER4 cell_num = (ni - 1) * (nj - 1) * (nk - 1);
    dynamic_array<double> x(node_num), y(node_num), z(node_num), density(node_num), velocity_x(node_num),
        velocity_y(node_num), velocity_z(node_num), pressure(node_num), iblank(node_num);
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
                index_type idx0 = idx_proxy(i + is, j + js, k + ks);
                density[idx] = data_manager->GetPrim(ID_DENSITY, idx0);
                velocity_x[idx] = data_manager->GetPrim(ID_VELOCITY_X, idx0);
                velocity_y[idx] = data_manager->GetPrim(ID_VELOCITY_Y, idx0);
                velocity_z[idx] = data_manager->GetPrim(ID_VELOCITY_Z, idx0);
                pressure[idx] = data_manager->GetPrim(ID_PRESSURE, idx0);
                iblank[idx] = (int)grid->GetIBlank(i + is, j + js, k + ks);
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
    double solution_time = GlobalData::GetDouble("iteration.current_time");
    INTEGER4 strandID = 1;
    INTEGER4 parentZn = 0;
    INTEGER4 isBlock = 1;
    INTEGER4 TotalNumFaceNodes = 1;
    INTEGER4 TotalNumBndryFaces = 1;
    INTEGER4 TotalNumBndryConnections = 1;
    INTEGER4 nFConns = 0;
    INTEGER4 FNMode = 0;
    int shrConn = 0;
    int i = TECZNE142((char *)zone_name.c_str(), &zone_type, &ni, &nj, &nk, &iCellMax, &jCellMax, &kCellMax,
                      &solution_time, &strandID, &parentZn, &isBlock, &nFConns, &FNMode, &TotalNumFaceNodes,
                      &TotalNumBndryFaces, &TotalNumBndryConnections, NULL, NULL, NULL, &shrConn);

    i = TECDAT142(&node_num, x.data(), &vIsDouble);
    i = TECDAT142(&node_num, y.data(), &vIsDouble);
    i = TECDAT142(&node_num, z.data(), &vIsDouble);
    i = TECDAT142(&node_num, density.data(), &vIsDouble);
    i = TECDAT142(&node_num, velocity_x.data(), &vIsDouble);
    i = TECDAT142(&node_num, velocity_y.data(), &vIsDouble);
    i = TECDAT142(&node_num, velocity_z.data(), &vIsDouble);
    i = TECDAT142(&node_num, pressure.data(), &vIsDouble);
    i = TECDAT142(&node_num, iblank.data(), &vIsDouble);
}
void zaran::Visual::WriteTecplotBinary(shared_ptr<NSFieldStruct> field)
{
    WriteTecplotBinary(field, HasVolumeFraction(field));
}

void zaran::Visual::WriteTecplotBinary(shared_ptr<NSFieldStruct> field, bool with_volume_fraction)
{
    auto data_manager = field->GetDataManager();
    auto grid = field->GetGrid();
    auto node = grid->GetNode();
    index_type is, ie, js, je, ks, ke;
    grid->GetRange(is, ie, js, je, ks, ke);
    INTEGER4 ni = ie - is + 1;
    INTEGER4 nj = je - js + 1;
    INTEGER4 nk = ke - ks + 1;
    index_type grid_ni = grid->GetNi();
    index_type grid_nj = grid->GetNj();
    index_type grid_nk = grid->GetNk();
    IdProxyStruct idx_proxy(grid_ni, grid_nj, grid_nk);
    INTEGER4 node_num = ni * nj * nk;
    INTEGER4 cell_num = (ni - 1) * (nj - 1) * (nk - 1);
    dynamic_array<double> x(node_num), y(node_num), z(node_num), density(node_num), velocity_x(node_num),
        velocity_y(node_num), velocity_z(node_num), pressure(node_num), iblank(node_num);
    dynamic_array<double> volume_fraction;
    if (with_volume_fraction)
    {
        volume_fraction.resize(node_num);
    }
    const double *volume_fraction_field = nullptr;
    if (with_volume_fraction)
    {
        double *volume_fraction_buffer = nullptr;
        field->GetData()->GetData("volume_fraction", volume_fraction_buffer);
        volume_fraction_field = volume_fraction_buffer;
    }
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
                index_type idx0 = idx_proxy(i + is, j + js, k + ks);
                density[idx] = data_manager->GetPrim(ID_DENSITY, idx0);
                velocity_x[idx] = data_manager->GetPrim(ID_VELOCITY_X, idx0);
                velocity_y[idx] = data_manager->GetPrim(ID_VELOCITY_Y, idx0);
                velocity_z[idx] = data_manager->GetPrim(ID_VELOCITY_Z, idx0);
                pressure[idx] = data_manager->GetPrim(ID_PRESSURE, idx0);
                if (volume_fraction_field != nullptr)
                {
                    volume_fraction[idx] = volume_fraction_field[idx0];
                }
                iblank[idx] =-1;
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
    double solution_time = GlobalData::GetDouble("iteration.current_time");
    INTEGER4 strandID = 1;
    INTEGER4 parentZn = 0;
    INTEGER4 isBlock = 1;
    INTEGER4 TotalNumFaceNodes = 1;
    INTEGER4 TotalNumBndryFaces = 1;
    INTEGER4 TotalNumBndryConnections = 1;
    INTEGER4 nFConns = 0;
    INTEGER4 FNMode = 0;
    int shrConn = 0;
    int i = TECZNE142((char *)zone_name.c_str(), &zone_type, &ni, &nj, &nk, &iCellMax, &jCellMax, &kCellMax,
                      &solution_time, &strandID, &parentZn, &isBlock, &nFConns, &FNMode, &TotalNumFaceNodes,
                      &TotalNumBndryFaces, &TotalNumBndryConnections, NULL, NULL, NULL, &shrConn);

    i = TECDAT142(&node_num, x.data(), &vIsDouble);
    i = TECDAT142(&node_num, y.data(), &vIsDouble);
    i = TECDAT142(&node_num, z.data(), &vIsDouble);
    i = TECDAT142(&node_num, density.data(), &vIsDouble);
    i = TECDAT142(&node_num, velocity_x.data(), &vIsDouble);
    i = TECDAT142(&node_num, velocity_y.data(), &vIsDouble);
    i = TECDAT142(&node_num, velocity_z.data(), &vIsDouble);
    i = TECDAT142(&node_num, pressure.data(), &vIsDouble);
    i = TECDAT142(&node_num, iblank.data(), &vIsDouble);
    if (!volume_fraction.empty())
    {
        i = TECDAT142(&node_num, volume_fraction.data(), &vIsDouble);
    }
}
void Visual::WriteTecASCII(shared_ptr<FieldManager> field_manager)
{
    std::string work_dir = GlobalData::GetString("work_dir");
    std::string file_name = "result\\" + std::to_string(GlobalData::GetInt("iteration.current_iter")) + ".dat";
    file_name = work_dir + "\\" + file_name;
    std::ofstream out(file_name);
    // out << "TITLE=\"Flow Field\"\n";
    // out <<
    // "VARIABLES=\"X\",\"Y\",\"Z\",\"Density\",\"Velocity_x\",\"Velocity_y\",\"Velocity_z\",\"Pressure\",\"iBlank\"\n";
    for (size_t iter_field = 0; iter_field < field_manager->GetFieldNum(); iter_field++)
    {
        auto field = field_manager->GetField(iter_field);
        auto field_type = field->GetFieldType();
        if (field_type == FieldType::NS_Structured)
        {
            auto field_struct = std::dynamic_pointer_cast<NSFieldStruct>(field);
            // 两相场多输出一个气相体积分数（单相场不输出，保持既有格式不变）
            WriteTecplotASCII(field_struct, out, HasVolumeFraction(field));
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
    // 变量名必须与各场实际写出的数据条数一致：两相场会多写一个体积分数，
    // 因此先扫一遍场列表，只要有一个 NS_Structured 场带体积分数就切换表头。
    bool with_volume_fraction = false;
    for (size_t iter_field = 0; iter_field < field_manager->GetFieldNum(); iter_field++)
    {
        auto probe = field_manager->GetField(iter_field);
        if (probe->GetFieldType() == FieldType::NS_Structured && HasVolumeFraction(probe))
        {
            with_volume_fraction = true;
            break;
        }
    }
    string var_name = "x, y, z, density, velocity_x, velocity_y, velocity_z, pressure, iBlank";
    if (with_volume_fraction)
    {
        var_name += ", volume_fraction";
    }
    std::string work_dir = GlobalData::GetString("work_dir");
    std::string file_name = std::to_string(GlobalData::GetInt("iteration.current_iter")) + ".plt";
    std::string temp_file_name = std::to_string(GlobalData::GetInt("iteration.current_iter")) + ".tmp";
    std::string result_dir = work_dir + "/result";
    int ierr = TECINI142(grid_name.c_str(), var_name.c_str(), temp_file_name.c_str(), (char *)".", &file_format,
                         &fileType, &debug, &vIsDouble);
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
            WriteTecplotBinary(field_struct, with_volume_fraction && HasVolumeFraction(field));
        }
        else
        {
            Log::warn("Field type is not supported!");
        }
    }
    ierr = TECEND142();
    if (ierr != 0)
    {
        Log::error("Tecplot file write error!");
        std::filesystem::remove(temp_file_name);
    }
    else
    {
        // 检查result目录是否存在
        if (!std::filesystem::exists(result_dir))
        {
            std::filesystem::create_directory(result_dir);
        }
        // 将临时文件拷贝到result目录
        std::filesystem::path temp_file_path(temp_file_name);
        std::filesystem::path file_path(result_dir + "\\" + file_name);
        std::filesystem::rename(temp_file_path, file_path);
        Log::info("Tecplot file write success!");
        Log::info("File name: {}", file_path.string());
    }
}

namespace zaran
{

void Visual::WriteParticleVTP(shared_ptr<FieldManager> field_manager, int iter)
{
    // 找到第一个 DEMField
    for (size_t i = 0; i < field_manager->GetFieldNum(); ++i)
    {
        auto f = std::dynamic_pointer_cast<DEMField>(field_manager->GetField(i));
        if (!f) continue;

        std::string work_dir    = GlobalData::GetString("work_dir");
        std::string result_folder = GlobalData::IsExist("output.result_folder")
                                        ? GlobalData::GetString("output.result_folder")
                                        : "result";
        std::string dir = work_dir + "/" + result_folder;
        std::string filename = dir + "/particles_" + std::to_string(iter) + ".vtp";
        WriteParticleVTP(*f->GetDEMData(), filename);
        return;
    }
    Log::warn("Visual::WriteParticleVTP: no DEMField found in FieldManager");
}

void Visual::WriteParticleVTP(const DEMFieldData& dem_data, const std::string& filename)
{
    const auto& particles = dem_data.GetParticles();
    const index_type N = particles.size();
    const auto& bonds = dem_data.GetBonds();
    const index_type B = bonds.size();

    // 性能：改为"整块字符串缓冲 + std::to_chars + 单次落盘"。
    // 旧实现逐字段 `fout <<`，万级粒子/键合下每帧需数十万次流格式化（实测约 200 ms/帧）。
    // 数值文本与旧实现保持可比（general 格式、6 位有效数字）。
    std::string buf;
    buf.reserve(static_cast<std::size_t>(N) * 520 + static_cast<std::size_t>(B) * 240 + 4096);

    auto num = [&](auto v) { AppendNumber(buf, v); };
    auto sp  = [&] { buf.push_back(' '); };
    auto nl  = [&] { buf.push_back('\n'); };
    auto raw = [&](const char* s) { buf.append(s); };

    auto point_scalar = [&](const char* type, const char* name, auto getter)
    {
        raw("        <DataArray type=\""); raw(type); raw("\" Name=\""); raw(name);
        raw("\" format=\"ascii\">\n");
        for (const auto& p : particles) { raw("          "); num(getter(p)); nl(); }
        raw("        </DataArray>\n");
    };
    auto point_vec3 = [&](const char* name, auto getter)
    {
        raw("        <DataArray type=\"Float64\" Name=\""); raw(name);
        raw("\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for (const auto& p : particles)
        {
            const Eigen::Vector3d v = getter(p);
            raw("          "); num(v.x()); sp(); num(v.y()); sp(); num(v.z()); nl();
        }
        raw("        </DataArray>\n");
    };
    auto cell_head = [&](const char* type, const char* name, const char* nc)
    {
        raw("        <DataArray type=\""); raw(type); raw("\" Name=\""); raw(name);
        raw("\""); raw(nc); raw(" format=\"ascii\">\n");
    };
    auto cell_tail = [&] { raw("        </DataArray>\n"); };

    raw("<?xml version=\"1.0\"?>\n");
    raw("<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    raw("  <PolyData>\n");
    raw("    <Piece NumberOfPoints=\""); num(N);
    raw("\" NumberOfVerts=\""); num(N);
    raw("\" NumberOfLines=\""); num(B);
    raw("\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n");

    // --- 点坐标 ---
    raw("      <Points>\n");
    raw("        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (const auto& p : particles)
    {
        raw("          "); num(p.pos.x()); sp(); num(p.pos.y()); sp(); num(p.pos.z()); nl();
    }
    raw("        </DataArray>\n");
    raw("      </Points>\n");

    // --- Verts (每个点作为一个 Vertex cell) ---
    raw("      <Verts>\n");
    raw("        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    raw("          ");
    for (index_type i = 0; i < N; ++i) { num(i); sp(); }
    raw("\n        </DataArray>\n");
    raw("        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    raw("          ");
    for (index_type i = 1; i <= N; ++i) { num(i); sp(); }
    raw("\n        </DataArray>\n");
    raw("      </Verts>\n");

    // 永久弹性键合以线单元输出。
    raw("      <Lines>\n");
    raw("        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for (const auto& bond : bonds)
    {
        raw("          "); num(bond.idx_a); sp(); num(bond.idx_b); nl();
    }
    raw("        </DataArray>\n");
    raw("        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    for (index_type i = 1; i <= B; ++i) { raw("          "); num(2 * i); nl(); }
    raw("        </DataArray>\n");
    raw("      </Lines>\n");

    // --- 弹簧/键合连接信息开关（output.bond_details，默认开启）---
    // 目的：让 ParaView 中能直接分辨"哪些线是弹簧连接、哪些已断裂、当前受力与伸长"。
    // 置 false 时输出与历史版本完全一致，便于超大规模键数时压缩文件体积。
    const bool bond_details =
        !GlobalData::IsExist("output.bond_details") || GlobalData::GetBool("output.bond_details");

    // 粒子级聚合量：与该粒子相连的键数、其中已断裂的条数及比例。
    // 用 bonds_broken_ratio 着色即可让裂纹带/断裂面在云图上直接凸显。
    std::vector<int> bond_incident(N, 0);
    std::vector<int> bond_broken_count(N, 0);
    if (bond_details)
    {
        for (const auto& bond : bonds)
        {
            if (bond.idx_a >= N || bond.idx_b >= N) continue;
            ++bond_incident[bond.idx_a];
            ++bond_incident[bond.idx_b];
            if (!bond.active)
            {
                ++bond_broken_count[bond.idx_a];
                ++bond_broken_count[bond.idx_b];
            }
        }
    }

    // --- 点数据 ---
    raw("      <PointData>\n");
    point_scalar("Float64", "radius", [](const DEMParticle& p) { return p.radius; });
    point_vec3("velocity", [](const DEMParticle& p) { return p.vel; });
    point_vec3("rotation", [](const DEMParticle& p) { return p.rotation; });
    point_vec3("omega", [](const DEMParticle& p) { return p.omega; });
    point_vec3("force", [](const DEMParticle& p) { return p.force; });
    point_scalar("Int64", "id", [](const DEMParticle& p) { return p.id; });
    point_scalar("Int32", "group", [](const DEMParticle& p) { return p.group; });
    point_scalar("Float64", "temperature", [](const DEMParticle& p) { return p.temperature; });
    point_scalar("Float64", "reaction_progress", [](const DEMParticle& p) { return p.reaction_progress; });
    point_scalar("Float64", "reaction_rate", [](const DEMParticle& p) { return p.reaction_rate; });
    point_scalar("Int32", "phase", [](const DEMParticle& p) { return p.phase; });
    point_scalar("Float64", "gas_temperature", [](const DEMParticle& p) { return p.gas_temperature; });
    point_scalar("Float64", "gas_pressure", [](const DEMParticle& p) { return p.gas_pressure; });
    point_scalar("Float64", "volume_ratio", [](const DEMParticle& p) { return p.volume_ratio; });
    point_scalar("Float64", "solid_volume", [](const DEMParticle& p) { return p.solid_volume; });
    point_scalar("Float64", "gas_volume", [](const DEMParticle& p) { return p.gas_volume; });
    point_scalar("Float64", "solid_core_radius", [](const DEMParticle& p) { return p.solid_core_radius; });
    point_scalar("Float64", "gas_internal_energy", [](const DEMParticle& p) { return p.gas_internal_energy; });
    point_scalar("Float64", "internal_heat_transfer", [](const DEMParticle& p) { return p.internal_heat_transfer; });
    point_scalar("Float64", "body_reaction_increment", [](const DEMParticle& p) { return p.body_reaction_increment; });
    point_scalar("Float64", "core_burn_increment", [](const DEMParticle& p) { return p.core_burn_increment; });
    point_scalar("Float64", "neighbor_burn_increment", [](const DEMParticle& p) { return p.neighbor_burn_increment; });
    if (bond_details)
    {
        raw("        <DataArray type=\"Int32\" Name=\"bonds_incident\" format=\"ascii\">\n");
        for (index_type i = 0; i < N; ++i) { raw("          "); num(bond_incident[i]); nl(); }
        raw("        </DataArray>\n");
        raw("        <DataArray type=\"Int32\" Name=\"bonds_broken\" format=\"ascii\">\n");
        for (index_type i = 0; i < N; ++i) { raw("          "); num(bond_broken_count[i]); nl(); }
        raw("        </DataArray>\n");
        raw("        <DataArray type=\"Float64\" Name=\"bonds_broken_ratio\" format=\"ascii\">\n");
        for (index_type i = 0; i < N; ++i)
        {
            raw("          ");
            num(bond_incident[i] > 0
                    ? static_cast<double>(bond_broken_count[i]) / static_cast<double>(bond_incident[i])
                    : 0.0);
            nl();
        }
        raw("        </DataArray>\n");
    }
    raw("      </PointData>\n");

    // --- 单元数据（前 N 个对应 Verts：占位；其后 B 个对应 Lines）---
    raw("      <CellData>\n");

    cell_head("Int64", "bond_id", "");
    for (index_type i = 0; i < N; ++i) raw("          -1\n");
    for (const auto& bond : bonds) { raw("          "); num(bond.id); nl(); }
    cell_tail();

    cell_head("Float64", "bond_extension", "");
    for (index_type i = 0; i < N; ++i) raw("          0\n");
    for (const auto& bond : bonds) { raw("          "); num(bond.extension); nl(); }
    cell_tail();

    cell_head("Float64", "bond_force_a", " NumberOfComponents=\"3\"");
    for (index_type i = 0; i < N; ++i) raw("          0 0 0\n");
    for (const auto& bond : bonds)
    {
        raw("          "); num(bond.force_a.x()); sp(); num(bond.force_a.y()); sp(); num(bond.force_a.z()); nl();
    }
    cell_tail();

    cell_head("Float64", "bond_heat_flow_a", "");
    for (index_type i = 0; i < N; ++i) raw("          0\n");
    for (const auto& bond : bonds) { raw("          "); num(bond.heat_flow_a); nl(); }
    cell_tail();

    cell_head("Float64", "bond_damage", "");
    for (index_type i = 0; i < N; ++i) raw("          0\n");
    for (const auto& bond : bonds) { raw("          "); num(bond.damage); nl(); }
    cell_tail();

    cell_head("Float64", "bond_elastic_energy", "");
    for (index_type i = 0; i < N; ++i) raw("          0\n");
    for (const auto& bond : bonds) { raw("          "); num(bond.elastic_energy); nl(); }
    cell_tail();

    cell_head("Float64", "bond_fracture_energy", "");
    for (index_type i = 0; i < N; ++i) raw("          0\n");
    for (const auto& bond : bonds) { raw("          "); num(bond.fracture_energy); nl(); }
    cell_tail();

    if (bond_details)
    {
        // 是否仍在承载：1 = 完好，0 = 已断裂。
        // ParaView 中对该数组做 Threshold(0,0) 即只显示断裂的连接。
        cell_head("Int32", "bond_active", "");
        for (index_type i = 0; i < N; ++i) raw("          0\n");
        for (const auto& bond : bonds) { raw("          "); num(bond.active ? 1 : 0); nl(); }
        cell_tail();

        // 连接来源：0 = 输入 bonds.csv 给定，1 = t=0 按几何自动生成的弹簧连接网络。
        cell_head("Int32", "bond_source", "");
        for (index_type i = 0; i < N; ++i) raw("          0\n");
        for (const auto& bond : bonds) { raw("          "); num(bond.source); nl(); }
        cell_tail();

        // 法向刚度 k_n (N/m)
        cell_head("Float64", "bond_stiffness", "");
        for (index_type i = 0; i < N; ++i) raw("          0\n");
        for (const auto& bond : bonds) { raw("          "); num(bond.normal_stiffness); nl(); }
        cell_tail();

        // 无应力长度 L0 (m)
        cell_head("Float64", "bond_rest_length", "");
        for (index_type i = 0; i < N; ++i) raw("          0\n");
        for (const auto& bond : bonds) { raw("          "); num(bond.rest_length); nl(); }
        cell_tail();

        // 当前轴向应变 extension/L0（正 = 拉伸，负 = 压缩）
        cell_head("Float64", "bond_strain", "");
        for (index_type i = 0; i < N; ++i) raw("          0\n");
        for (const auto& bond : bonds)
        {
            raw("          ");
            num(bond.rest_length > 0.0 ? bond.extension / bond.rest_length : 0.0);
            nl();
        }
        cell_tail();

        // 历史最大拉伸应变：与 dem.bond_break_strain / dem.spring_network_fracture_strain
        // 直接可比，用来判断"离断裂还有多远"。
        cell_head("Float64", "bond_max_strain", "");
        for (index_type i = 0; i < N; ++i) raw("          0\n");
        for (const auto& bond : bonds) { raw("          "); num(bond.maximum_tensile_strain); nl(); }
        cell_tail();

        // 键合力大小 |F_a| (N)：断裂后恒为 0，可直接作为着色变量
        cell_head("Float64", "bond_force_magnitude", "");
        for (index_type i = 0; i < N; ++i) raw("          0\n");
        for (const auto& bond : bonds) { raw("          "); num(bond.force_a.norm()); nl(); }
        cell_tail();

        // 断裂时耗散的能量 (J)：> 0 即真实发生过断裂
        cell_head("Float64", "bond_energy_loss", "");
        for (index_type i = 0; i < N; ++i) raw("          0\n");
        for (const auto& bond : bonds) { raw("          "); num(bond.dissipated_fracture_energy); nl(); }
        cell_tail();

        // 历史最大压缩应变：与"压缩断裂阈值 = ratio × 抗拉阈值 × strength_scale"直接可比，
        // 用来判断破坏是否由压缩（粉碎）而非拉伸（劈裂）主导。
        cell_head("Float64", "bond_max_compression", "");
        for (index_type i = 0; i < N; ++i) raw("          0\n");
        for (const auto& bond : bonds) { raw("          "); num(bond.maximum_compressive_strain); nl(); }
        cell_tail();

        // 逐键强度折减系数（Weibull 异质性；1 = 均质）。有效阈值为名义值乘以该系数。
        cell_head("Float64", "bond_strength_scale", "");
        for (index_type i = 0; i < N; ++i) raw("          1\n");
        for (const auto& bond : bonds) { raw("          "); num(bond.strength_scale); nl(); }
        cell_tail();

        // 切向弹簧是否已独立断裂（LSM breakmod=4）：1 = 切向已断、只剩法向承载。
        // 与 bond_active 不同 —— 切向断裂**不**断开整条键，只是退出切向通道并把
        // 剩余法向断裂能按 UtIII 折减。用该数组着色可看到"切向先坏"的损伤前驱区。
        cell_head("Int32", "bond_shear_broken", "");
        for (index_type i = 0; i < N; ++i) raw("          0\n");
        for (const auto& bond : bonds) { raw("          "); num(bond.shear_broken ? 1 : 0); nl(); }
        cell_tail();

        // 失效模式（LSM break_mod）：0 = 未失效，1 = 拉伸混合，2 = 压剪，3 = 压缩压溃，
        // 4 = 切向独立断裂。用它可以直接回答"这片区域是被拉开的还是被剪/压碎的"。
        cell_head("Int32", "bond_break_mode", "");
        for (index_type i = 0; i < N; ++i) raw("          0\n");
        for (const auto& bond : bonds) { raw("          "); num(bond.break_mode); nl(); }
        cell_tail();
    }

    raw("      </CellData>\n");
    raw("    </Piece>\n");
    raw("  </PolyData>\n");
    raw("</VTKFile>\n");

    std::ofstream fout(filename, std::ios::binary);
    if (!fout.is_open())
    {
        Log::warn("Visual::WriteParticleVTP: cannot open '{}'", filename);
        return;
    }
    fout.write(buf.data(), static_cast<std::streamsize>(buf.size()));

    Log::info("Visual::WriteParticleVTP: {} particles written to '{}'", N, filename);
}

// ==================================================================
// DEM 输出：Tecplot ASCII
//
// 为什么拆成两个文件：
//   Tecplot ASCII 的文件头 `VARIABLES` 是**整个文件共享**的，不是每个 zone 各自一套。
//   粒子是点云（惰性算例 25 个变量、含能算例 45 个），键合是线元（3 个节点量 + 20 个
//   单元中心量，共 23 个），两者的变量表不同；硬塞进一个文件只能用大量空列去凑，
//   既浪费体积又难读。所以粒子、键合各写一个 .dat，流场（网格）本来就有自己的 .dat。
// ==================================================================
namespace
{
	/// @brief 解析粒子点云的 ZONE 行写法（控制文件 output.particle_zone_keywords）
	/// @details Tecplot ASCII 有两代关键字风格：
	///            - 传统（Tecplot 10 / Focus）：有限元 zone 用 `F=FEPOINT` + `N=` + `E=`
	///              元素类型 `ET=TRIANGLE|QUADRILATERAL|TETRAHEDRON|BRICK`；
	///            - 新版（Tecplot 360 之后）：`ZONETYPE=FEPOINT` + `DATAPACKING=POINT`。
	///          两代混写时有些 Tecplot 版本会直接拒绝该 zone（甚至整个文件），
	///          因此这里给三档：
	///            "none"    （默认）ZONE 行只留 T= 与 SOLUTIONTIME=，最小化关键字；
	///            "modern"  ZONETYPE=FEPOINT + DATAPACKING=POINT（360 系）；
	///            "classic" F=FEPOINT + N= + E=0（传统 FE 写法）。
	///          默认取 "none" —— 用户实测其 Tecplot 只认这一种。
	enum class ParticleZoneStyle
	{
		None,
		Modern,
		Classic,
	};

	ParticleZoneStyle ResolveParticleZoneStyle()
	{
		std::string mode = "none";
		if (GlobalData::IsExist("output.particle_zone_keywords"))
		{
			mode = GlobalData::GetString("output.particle_zone_keywords");
		}
		if (mode == "none") { return ParticleZoneStyle::None; }
		if (mode == "modern") { return ParticleZoneStyle::Modern; }
		if (mode == "classic") { return ParticleZoneStyle::Classic; }
		Log::warn("output.particle_zone_keywords='{}' 无法识别（应为 none|modern|classic），"
			"按 none 处理", mode);
		return ParticleZoneStyle::None;
	}

	/// @brief 解析输出精度（控制文件 output.tecplot_precision，默认 12 位有效数字）
	int ResolveTecplotPrecision()
	{
		if (GlobalData::IsExist("output.tecplot_precision"))
		{
			const int p = GlobalData::GetInt("output.tecplot_precision");
			if (p >= 6 && p <= 17)
			{
				return p;
			}
			Log::warn("output.tecplot_precision={} 超出允许范围 [6,17]，改用默认值 12", p);
		}
		return 12;
	}

	/// @brief 解析粒子变量集（控制文件 output.particle_fields = auto | core | full）
	/// @details auto（默认）= 只有出现含能粒子时才附加热化学列。惰性算例因此少 20 列，
	///          既省体积，也不会出现一整列常数的无用字段。
	bool ResolveParticleFieldsFull(const DEMFieldData& dem_data)
	{
		std::string mode = "auto";
		if (GlobalData::IsExist("output.particle_fields"))
		{
			mode = GlobalData::GetString("output.particle_fields");
		}
		if (mode == "full") { return true; }
		if (mode == "core") { return false; }
		if (mode != "auto")
		{
			Log::warn("output.particle_fields='{}' 无法识别（应为 auto|core|full），按 auto 处理", mode);
		}
		const auto& particles = dem_data.GetParticles();
		return std::any_of(particles.begin(), particles.end(),
			[](const DEMParticle& p) { return p.energetic; });
	}

	/// @brief Tecplot ASCII 缓冲区：统一处理分隔符、换行与非有限值
	class TecplotBuffer
	{
	public:
		explicit TecplotBuffer(int precision) : m_precision(precision) {}

		void Reserve(std::size_t n) { m_buf.reserve(n); }
		void Raw(const char* s) { m_buf.append(s); }
		void Raw(const std::string& s) { m_buf.append(s); }
		void Sp() { m_buf.push_back(' '); }
		void Nl() { m_buf.push_back('\n'); }

		/// @brief 写一个浮点数（general 格式、m_precision 位有效数字）
		/// @details 非有限值（NaN/Inf）必须替换掉：Tecplot 解析不了 "nan"/"inf"，
		///          一旦写出会让**整个文件**打不开，而不只是那一个值不对。
		void Real(double v)
		{
			if (!std::isfinite(v))
			{
				++m_non_finite;
				v = 0.0;
			}
			char tmp[64];
			const auto r = std::to_chars(tmp, tmp + sizeof(tmp), v,
				std::chars_format::general, m_precision);
			m_buf.append(tmp, static_cast<std::size_t>(r.ptr - tmp));
		}

		void Int(long long v)
		{
			char tmp[32];
			const auto r = std::to_chars(tmp, tmp + sizeof(tmp), v);
			m_buf.append(tmp, static_cast<std::size_t>(r.ptr - tmp));
		}

		int NonFinite() const { return m_non_finite; }
		const std::string& Str() const { return m_buf; }

	private:
		std::string m_buf;
		int m_precision = 12;
		int m_non_finite = 0;
	};

	/// @brief 写粒子点云的 ZONE 行（三种风格共用一份代码，避免两处漂移）
	/// @details `none`（默认）只写 T= 与 SOLUTIONTIME= —— 用户实测其 Tecplot 只认这一种；
	///          需要 FE 关键字时用 `output.particle_zone_keywords = modern|classic`。
	void WriteParticleZoneLine(TecplotBuffer& out, const std::string& title,
		double solution_time, long long num_points, ParticleZoneStyle style)
	{
		out.Raw("ZONE T=\""); out.Raw(title); out.Raw("\"");
		if (style == ParticleZoneStyle::Modern)
		{
			out.Raw(", STRANDID=1");
		}
		out.Raw(", SOLUTIONTIME="); out.Real(solution_time);
		if (style == ParticleZoneStyle::Classic)
		{
			out.Raw(", N="); out.Int(num_points);
			out.Raw(", E=0, F=FEPOINT");
		}
		else if (style == ParticleZoneStyle::Modern)
		{
			out.Raw(", N="); out.Int(num_points);
			out.Raw(", E=0, ZONETYPE=FEPOINT, DATAPACKING=POINT");
		}
		out.Raw("\n");
	}

	/// @brief 写一个 BLOCK 数据块（每 per_line 个值换行；长常量段压缩成 `Rep*Num`）
	/// @details Tecplot ASCII 支持 `Rep*Num` 记法（如 `20000*1`）。键合的 active /
	///          source / 常刚度这类整块常量用它能把文件缩小一个量级。
	///          粒子走的是 POINT 打包（一行一粒子），不经这里，以保证脚本可直接解析。
	template <typename Getter>
	void WriteTecplotBlock(TecplotBuffer& out, std::size_t count, Getter&& value,
		std::size_t per_line = 8, std::size_t min_run = 8)
	{
		std::size_t col = 0;
		auto sep = [&]()
		{
			if (col == per_line) { out.Nl(); col = 0; }
			else if (col > 0) { out.Sp(); }
		};
		std::size_t i = 0;
		while (i < count)
		{
			const double v = value(i);
			std::size_t j = i + 1;
			while (j < count && value(j) == v) { ++j; }
			if (j - i >= min_run)
			{
				sep();
				out.Int(static_cast<long long>(j - i));
				out.Raw("*");
				out.Real(v);
				++col;
			}
			else
			{
				for (std::size_t k = i; k < j; ++k)
				{
					sep();
					out.Real(value(k));
					++col;
				}
			}
			i = j;
		}
		if (col > 0) { out.Nl(); }
	}

	/// @brief 由变量名列表生成 Tecplot 的 VARIABLES 行
	std::string TecplotVariableLine(const std::vector<std::string>& names)
	{
		std::string s = "VARIABLES=";
		for (std::size_t i = 0; i < names.size(); ++i)
		{
			if (i > 0) { s += ","; }
			s += "\"";
			s += names[i];
			s += "\"";
		}
		s += "\n";
		return s;
	}

	/// @brief 落盘（文本模式，与流场 ASCII 输出保持一致）
	bool FlushTecplotFile(const std::string& filename, const std::string& buf)
	{
		std::ofstream out(filename);
		if (!out.is_open())
		{
			Log::warn("Visual: cannot open '{}' for writing", filename);
			return false;
		}
		out.write(buf.data(), static_cast<std::streamsize>(buf.size()));
		return out.good();
	}
} // namespace

void Visual::WriteParticleTecplotASCII(shared_ptr<FieldManager> field_manager, int iter,
	double solution_time, bool with_bonds)
{
	// 找到第一个 DEMField
	for (size_t i = 0; i < field_manager->GetFieldNum(); ++i)
	{
		auto f = std::dynamic_pointer_cast<DEMField>(field_manager->GetField(i));
		if (!f) continue;

		std::string work_dir = GlobalData::GetString("work_dir");
		std::string result_folder = GlobalData::IsExist("output.result_folder")
			? GlobalData::GetString("output.result_folder")
			: "result";
		std::string dir = work_dir + "/" + result_folder;
		WriteParticleTecplotASCII(*f->GetDEMData(),
			dir + "/particles_" + std::to_string(iter) + ".dat",
			dir + "/bonds_" + std::to_string(iter) + ".dat",
			solution_time, with_bonds);
		return;
	}
	Log::warn("Visual::WriteParticleTecplotASCII: no DEMField found in FieldManager");
}

void Visual::WriteParticleTecplotASCII(const DEMFieldData& dem_data,
	const std::string& particle_file, const std::string& bond_file,
	double solution_time, bool with_bonds)
{
	const auto& particles = dem_data.GetParticles();
	const auto& bonds = dem_data.GetBonds();
	const std::size_t n_particle = particles.size();
	const std::size_t n_bond = bonds.size();
	const int precision = ResolveTecplotPrecision();
	const bool full_fields = ResolveParticleFieldsFull(dem_data);

	if (n_particle == 0)
	{
		Log::warn("Visual::WriteParticleTecplotASCII: 粒子数为 0，跳过输出（N=0 的 zone 在 Tecplot 里非法）");
		return;
	}

	// ---- 粒子级聚合量：相邻键数 / 已断键数（与 VTP 输出口径一致）----
	std::vector<int> bond_incident(n_particle, 0);
	std::vector<int> bond_broken(n_particle, 0);
	for (const auto& bond : bonds)
	{
		const index_type ends[2] = { bond.idx_a, bond.idx_b };
		for (const index_type e : ends)
		{
			if (e >= static_cast<index_type>(n_particle)) { continue; }
			++bond_incident[e];
			if (!bond.active) { ++bond_broken[e]; }
		}
	}

	// ---- 粒子变量表：表头与数据出自同一份列表，杜绝列错位 ----
	using ParticleGetter = std::function<double(const DEMParticle&)>;
	std::vector<std::pair<std::string, ParticleGetter>> fields;
	fields.reserve(full_fields ? 48 : 32);
	auto add = [&fields](const char* name, ParticleGetter g)
	{
		fields.emplace_back(name, std::move(g));
	};
	add("X", [](const DEMParticle& p) { return p.pos.x(); });
	add("Y", [](const DEMParticle& p) { return p.pos.y(); });
	add("Z", [](const DEMParticle& p) { return p.pos.z(); });
	add("id", [](const DEMParticle& p) { return static_cast<double>(p.id); });
	add("group", [](const DEMParticle& p) { return static_cast<double>(p.group); });
	add("radius", [](const DEMParticle& p) { return p.radius; });
	add("mass", [](const DEMParticle& p) { return p.mass; });
	add("active", [](const DEMParticle& p) { return p.active ? 1.0 : 0.0; });
	add("kinematic", [](const DEMParticle& p) { return p.kinematic ? 1.0 : 0.0; });
	add("velocity_x", [](const DEMParticle& p) { return p.vel.x(); });
	add("velocity_y", [](const DEMParticle& p) { return p.vel.y(); });
	add("velocity_z", [](const DEMParticle& p) { return p.vel.z(); });
	add("omega_x", [](const DEMParticle& p) { return p.omega.x(); });
	add("omega_y", [](const DEMParticle& p) { return p.omega.y(); });
	add("omega_z", [](const DEMParticle& p) { return p.omega.z(); });
	add("force_x", [](const DEMParticle& p) { return p.force.x(); });
	add("force_y", [](const DEMParticle& p) { return p.force.y(); });
	add("force_z", [](const DEMParticle& p) { return p.force.z(); });
	add("torque_x", [](const DEMParticle& p) { return p.torque.x(); });
	add("torque_y", [](const DEMParticle& p) { return p.torque.y(); });
	add("torque_z", [](const DEMParticle& p) { return p.torque.z(); });
	add("temperature", [](const DEMParticle& p) { return p.temperature; });

	if (full_fields)
	{
		add("energetic", [](const DEMParticle& p) { return p.energetic ? 1.0 : 0.0; });
		add("phase", [](const DEMParticle& p) { return static_cast<double>(p.phase); });
		add("rotation_x", [](const DEMParticle& p) { return p.rotation.x(); });
		add("rotation_y", [](const DEMParticle& p) { return p.rotation.y(); });
		add("rotation_z", [](const DEMParticle& p) { return p.rotation.z(); });
		add("reaction_progress", [](const DEMParticle& p) { return p.reaction_progress; });
		add("reaction_rate", [](const DEMParticle& p) { return p.reaction_rate; });
		add("gas_temperature", [](const DEMParticle& p) { return p.gas_temperature; });
		add("gas_pressure", [](const DEMParticle& p) { return p.gas_pressure; });
		add("volume_ratio", [](const DEMParticle& p) { return p.volume_ratio; });
		add("solid_volume", [](const DEMParticle& p) { return p.solid_volume; });
		add("gas_volume", [](const DEMParticle& p) { return p.gas_volume; });
		add("solid_core_radius", [](const DEMParticle& p) { return p.solid_core_radius; });
		add("gas_internal_energy", [](const DEMParticle& p) { return p.gas_internal_energy; });
		add("internal_heat_transfer", [](const DEMParticle& p) { return p.internal_heat_transfer; });
		add("body_reaction_increment", [](const DEMParticle& p) { return p.body_reaction_increment; });
		add("core_burn_increment", [](const DEMParticle& p) { return p.core_burn_increment; });
		add("neighbor_burn_increment", [](const DEMParticle& p) { return p.neighbor_burn_increment; });
	}

	// 末三列由下标决定，单独补上（保证它们总在末尾）
	std::vector<std::string> extra_names{ "bonds_incident", "bonds_broken", "bonds_broken_ratio" };

	// ---- 写粒子文件 ----
	{
		TecplotBuffer out(precision);
		out.Reserve(n_particle * ((fields.size() + extra_names.size()) * 16 + 32) + 1024);

		out.Raw("TITLE=\"Zaran3 DEM particles\"\n");
		std::vector<std::string> names;
		names.reserve(fields.size() + extra_names.size());
		for (const auto& f : fields) { names.push_back(f.first); }
		for (const auto& n : extra_names) { names.push_back(n); }
		out.Raw(TecplotVariableLine(names));

		WriteParticleZoneLine(out, "Particles", solution_time,
			static_cast<long long>(n_particle), ResolveParticleZoneStyle());

		// 一行一个粒子：全部变量按 VARIABLES 顺序排列
		for (index_type i = 0; i < static_cast<index_type>(n_particle); ++i)
		{
			const DEMParticle& p = particles[i];
			for (const auto& f : fields)
			{
				out.Real(f.second(p));
				out.Sp();
			}
			out.Real(static_cast<double>(bond_incident[i]));
			out.Sp();
			out.Real(static_cast<double>(bond_broken[i]));
			out.Sp();
			out.Real(bond_incident[i] > 0
				? static_cast<double>(bond_broken[i]) / static_cast<double>(bond_incident[i])
				: 0.0);
			out.Nl();
		}

		if (FlushTecplotFile(particle_file, out.Str()))
		{
			Log::info("Visual::WriteParticleTecplotASCII: {} particles, {} variables -> '{}'",
				n_particle, fields.size() + extra_names.size(), particle_file);
		}
		if (out.NonFinite() > 0)
		{
			Log::warn("Visual: 粒子输出中有 {} 个非有限值被写成 0（Tecplot 无法解析 nan/inf）",
				out.NonFinite());
		}
	}

	// ---- 写键合文件 ----
	if (!with_bonds || n_bond == 0)
	{
		return;
	}

	// 只保留两端都有效的键；给它们涉及的粒子重新编号，节点数远小于 2E
	std::vector<index_type> node_of(n_particle, static_cast<index_type>(-1));
	std::vector<index_type> node_particle;
	std::vector<std::size_t> valid;
	valid.reserve(n_bond);
	node_particle.reserve(std::min<std::size_t>(n_particle, 2 * n_bond));
	for (std::size_t b = 0; b < n_bond; ++b)
	{
		const DEMBond& bond = bonds[b];
		if (bond.idx_a >= static_cast<index_type>(n_particle)) { continue; }
		if (bond.idx_b >= static_cast<index_type>(n_particle)) { continue; }
		valid.push_back(b);
		const index_type ends[2] = { bond.idx_a, bond.idx_b };
		for (const index_type e : ends)
		{
			if (node_of[e] == static_cast<index_type>(-1))
			{
				node_of[e] = static_cast<index_type>(node_particle.size());
				node_particle.push_back(e);
			}
		}
	}
	const std::size_t n_elem = valid.size();
	const std::size_t n_node = node_particle.size();
	if (n_elem == 0)
	{
		Log::warn("Visual::WriteParticleTecplotASCII: 没有两端都有效的键，跳过键合文件");
		return;
	}

	// ---- 单元中心变量表（顺序 = VARIABLES 里第 4 个起的顺序）----
	using BondGetter = std::function<double(const DEMBond&)>;
	std::vector<std::pair<std::string, BondGetter>> cell;
	cell.reserve(24);
	auto addc = [&cell](const char* name, BondGetter g)
	{
		cell.emplace_back(name, std::move(g));
	};
	addc("bond_id", [](const DEMBond& b) { return static_cast<double>(b.id); });
	addc("active", [](const DEMBond& b) { return b.active ? 1.0 : 0.0; });
	addc("source", [](const DEMBond& b) { return static_cast<double>(b.source); });
	addc("normal_stiffness", [](const DEMBond& b) { return b.normal_stiffness; });
	addc("tangential_stiffness", [](const DEMBond& b) { return b.tangential_stiffness; });
	addc("rest_length", [](const DEMBond& b) { return b.rest_length; });
	addc("extension", [](const DEMBond& b) { return b.extension; });
	addc("strain", [](const DEMBond& b)
		{ return b.rest_length > 0.0 ? b.extension / b.rest_length : 0.0; });
	addc("max_strain", [](const DEMBond& b) { return b.maximum_tensile_strain; });
	addc("max_compression", [](const DEMBond& b) { return b.maximum_compressive_strain; });
	addc("strength_scale", [](const DEMBond& b) { return b.strength_scale; });
	addc("force_x", [](const DEMBond& b) { return b.force_a.x(); });
	addc("force_y", [](const DEMBond& b) { return b.force_a.y(); });
	addc("force_z", [](const DEMBond& b) { return b.force_a.z(); });
	addc("force_magnitude", [](const DEMBond& b) { return b.force_a.norm(); });
	addc("damage", [](const DEMBond& b) { return b.damage; });
	addc("elastic_energy", [](const DEMBond& b) { return b.elastic_energy; });
	addc("fracture_energy", [](const DEMBond& b) { return b.fracture_energy; });
	addc("energy_loss", [](const DEMBond& b) { return b.dissipated_fracture_energy; });
	addc("heat_flow_a", [](const DEMBond& b) { return b.heat_flow_a; });
	// γ 控制内聚律的两个特征分离度（未启用时为 0）：用于从输出直接复核
	// "曲线下面积 == γ·A" 与 "δ_f = 2·Efract/F_p"，不必再去猜内部参数。
	addc("cohesive_delta_p", [](const DEMBond& b) { return b.cohesive_peak_separation; });
	addc("cohesive_delta_f", [](const DEMBond& b) { return b.cohesive_failure_separation; });
	// 多通道失效诊断（默认关闭时恒为 0）：切向独立断裂标志与失效模式。
	// 用 break_mode 着色可以直接分辨"被拉开的"（1）/"被剪坏的"（2）/"被压碎的"（3）。
	addc("shear_broken", [](const DEMBond& b) { return b.shear_broken ? 1.0 : 0.0; });
	addc("shear_energy_loss", [](const DEMBond& b) { return b.dissipated_shear_energy; });
	addc("break_mode", [](const DEMBond& b) { return static_cast<double>(b.break_mode); });
	// 累积切向位移 |δ_t| (m)：切向独立断裂判据 Et = ½·k_t·|δ_t|² ≥ UtIII 的直接观测量
	// （切向断裂后恒为 0）。
	addc("tangential_displacement", [](const DEMBond& b) { return b.delta_t.norm(); });
	// 法向断裂能预算的折减系数（LSM 的 GBratio，初值 1）。切向独立断裂时乘 (1−α)。
	// 有效法向断裂能 = fracture_energy × cohesive_budget_scale。
	addc("cohesive_budget_scale", [](const DEMBond& b) { return b.cohesive_budget_scale; });

	TecplotBuffer out(precision);
	out.Reserve(n_elem * (cell.size() * 16 + 24) + n_node * 64 + 1024);

	out.Raw("TITLE=\"Zaran3 DEM bonds\"\n");
	std::vector<std::string> names{ "X", "Y", "Z" };
	names.reserve(3 + cell.size());
	for (const auto& c : cell) { names.push_back(c.first); }
	out.Raw(TecplotVariableLine(names));

	out.Raw("ZONE T=\"Bonds\", STRANDID=1, SOLUTIONTIME=");
	out.Real(solution_time);
	out.Raw(", N="); out.Int(static_cast<long long>(n_node));
	out.Raw(", E="); out.Int(static_cast<long long>(n_elem));
	out.Raw(", ZONETYPE=FELINESEG, DATAPACKING=BLOCK");
	// Tecplot 规定：带单元中心量的 zone 必须用 BLOCK，且数据块顺序为
	// 「全部节点变量（按 VARIABLES 顺序）→ 全部单元中心变量 → 连接表」
	out.Raw(", VARLOCATION=([4-");
	out.Int(static_cast<long long>(3 + cell.size()));
	out.Raw("]=CELLCENTERED)");
	out.Nl();

	// 节点变量块：X、Y、Z
	WriteTecplotBlock(out, n_node,
		[&](std::size_t k) { return particles[node_particle[k]].pos.x(); });
	WriteTecplotBlock(out, n_node,
		[&](std::size_t k) { return particles[node_particle[k]].pos.y(); });
	WriteTecplotBlock(out, n_node,
		[&](std::size_t k) { return particles[node_particle[k]].pos.z(); });

	// 单元中心变量块：顺序必须与 VARIABLES 严格一致
	for (const auto& c : cell)
	{
		WriteTecplotBlock(out, n_elem,
			[&](std::size_t k) { return c.second(bonds[valid[k]]); });
	}

	// 连接表（1 基）必须跟在所有数据之后
	for (const std::size_t b : valid)
	{
		out.Int(static_cast<long long>(node_of[bonds[b].idx_a]) + 1);
		out.Sp();
		out.Int(static_cast<long long>(node_of[bonds[b].idx_b]) + 1);
		out.Nl();
	}

	if (FlushTecplotFile(bond_file, out.Str()))
	{
		Log::info("Visual::WriteParticleTecplotASCII: {} bonds ({} nodes, {} cell variables) -> '{}'",
			n_elem, n_node, cell.size(), bond_file);
	}
	if (out.NonFinite() > 0)
	{
		Log::warn("Visual: 键合输出中有 {} 个非有限值被写成 0", out.NonFinite());
	}
}

void Visual::WriteUniformGridScalarTecplotASCII(const std::string& filename,
	const std::string& title, const std::string& zone_name,
	int ni, int nj, int nk,
	double x0, double y0, double z0,
	double dx, double dy, double dz,
	const std::vector<UniformScalarField>& fields,
	double solution_time)
{
	const std::size_t n_node = static_cast<std::size_t>(ni) * nj * nk;
	const int precision = ResolveTecplotPrecision();

	// 变量表：X,Y,Z + 各标量场。表头与数据出自同一份列表 ⇒ 不会列错位。
	std::vector<std::string> names{ "X", "Y", "Z" };
	for (const auto& f : fields)
	{
		if (!f.values || f.values->size() != n_node)
		{
			throw ZaranError("Visual::WriteUniformGridScalarTecplotASCII: field '" + f.name
				+ "' has " + std::to_string(f.values ? f.values->size() : 0)
				+ " values but the grid has " + std::to_string(n_node) + " nodes");
		}
		names.push_back(f.name);
	}

	TecplotBuffer out(precision);
	out.Reserve(n_node * (names.size() * 16 + 24) + 1024);
	out.Raw("TITLE=\""); out.Raw(title); out.Raw("\"\n");
	out.Raw(TecplotVariableLine(names));
	out.Raw("ZONE T=\""); out.Raw(zone_name);
	out.Raw("\", STRANDID=1, SOLUTIONTIME=");
	out.Real(solution_time);
	out.Raw(", I="); out.Int(ni);
	out.Raw(", J="); out.Int(nj);
	out.Raw(", K="); out.Int(nk);
	out.Raw(", ZONETYPE=ORDERED, DATAPACKING=POINT\n");

	// ORDERED + POINT 的数据顺序：k 最慢、i 最快
	for (int k = 0; k < nk; ++k)
	{
		for (int j = 0; j < nj; ++j)
		{
			for (int i = 0; i < ni; ++i)
			{
				const std::size_t idx = static_cast<std::size_t>(i)
					+ static_cast<std::size_t>(ni) * (static_cast<std::size_t>(j)
						+ static_cast<std::size_t>(nj) * static_cast<std::size_t>(k));
				out.Real(x0 + i * dx); out.Sp();
				out.Real(y0 + j * dy); out.Sp();
				out.Real(z0 + k * dz);
				for (const auto& f : fields)
				{
					out.Sp();
					out.Real((*f.values)[idx]);
				}
				out.Nl();
			}
		}
	}

	if (FlushTecplotFile(filename, out.Str()))
	{
		Log::info("Visual::WriteUniformGridScalarTecplotASCII: {} nodes, {} variables -> '{}'",
			n_node, names.size(), filename);
	}
	if (out.NonFinite() > 0)
	{
		Log::warn("Visual: '{}' 中有 {} 个非有限值被写成 0", filename, out.NonFinite());
	}
}

void Visual::WritePointsTecplotASCII(const std::string& filename,
	const std::string& title, const std::string& zone_name,
	int num_points, const std::vector<PointScalarField>& fields,
	double solution_time)
{
	if (num_points <= 0)
	{
		Log::warn("Visual::WritePointsTecplotASCII: 点数为 0，跳过输出（N=0 的 zone 在 Tecplot 里非法）");
		return;
	}
	if (fields.empty())
	{
		Log::warn("Visual::WritePointsTecplotASCII: 没有任何变量，跳过输出");
		return;
	}
	const int precision = ResolveTecplotPrecision();
	std::vector<std::string> names;
	names.reserve(fields.size());
	for (const auto& f : fields)
	{
		if (!f.values || f.values->size() != static_cast<std::size_t>(num_points))
		{
			throw ZaranError("Visual::WritePointsTecplotASCII: field '" + f.name + "' has "
				+ std::to_string(f.values ? f.values->size() : 0) + " values but there are "
				+ std::to_string(num_points) + " points");
		}
		names.push_back(f.name);
	}

	TecplotBuffer out(precision);
	out.Reserve(static_cast<std::size_t>(num_points) * (names.size() * 18 + 8) + 512);
	out.Raw("TITLE=\""); out.Raw(title); out.Raw("\"\n");
	out.Raw(TecplotVariableLine(names));
	WriteParticleZoneLine(out, zone_name, solution_time, num_points, ResolveParticleZoneStyle());

	// POINT 打包：一行一个点；每 512 个值插一个换行，保证单行远小于 32000 字节
	const std::size_t kMaxValuesPerLine = 512;
	std::size_t column = 0;
	for (int p = 0; p < num_points; ++p)
	{
		for (std::size_t v = 0; v < fields.size(); ++v)
		{
			if (column == 0)
			{
				// 行首不写分隔空格
			}
			else if (column >= kMaxValuesPerLine)
			{
				out.Nl();
				column = 0;
			}
			else
			{
				out.Sp();
			}
			out.Real((*fields[v].values)[static_cast<std::size_t>(p)]);
			++column;
		}
		out.Nl();
		column = 0;
	}

	if (FlushTecplotFile(filename, out.Str()))
	{
		Log::info("Visual::WritePointsTecplotASCII: {} points, {} variables -> '{}'",
			num_points, names.size(), filename);
	}
	if (out.NonFinite() > 0)
	{
		Log::warn("Visual: '{}' 中有 {} 个非有限值被写成 0", filename, out.NonFinite());
	}
}


} // namespace zaran
