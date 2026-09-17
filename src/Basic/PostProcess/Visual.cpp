#include "Visual.h"
#include "GlobalData.h"
#include "Log.h"
#include "NSFieldFN.h"
#include "DEMField.h"
#include "FastNumberFormat.h"
#include <TECIO.h>
#include <cgnslib.h>
#include <filesystem>
#include <fstream>
#include <string>
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

void zaran::Visual::WriteTecplotASCII(shared_ptr<NSFieldStruct> field, std::ostream &os)
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
    // os <<
    // "VARIABLES=\"X\",\"Y\",\"Z\",\"Density\",\"Velocity_x\",\"Velocity_y\",\"Velocity_z\",\"Pressure\",\"x_zeta\",\"y_zeta\",\"z_xi\",\"ERR_J\"\n";
    os << "VARIABLES=\"X\",\"Y\",\"Z\",\"Density\",\"Velocity_x\",\"Velocity_y\",\"Velocity_z\",\"Pressure\"\n"; // os
                                                                                                                 // <<
                                                                                                                 // "VARIABLES=\"X\",\"Y\",\"Z\",\"Density\",\"Velocity_x\",\"Velocity_y\",\"Velocity_z\",\"Pressure,\"zeta_z\n";
    os << "ZONE I=" << ni << ", J=" << nj << ", K=" << nk << ", F=POINT\n";
    double solution_time = GlobalData::GetDouble("iteration.current_time");
    os << "SOLUTIONTIME=" << solution_time << "\n";

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
                   << pressure[idx] << "\n";
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
    string var_name = "x, y, z, density, velocity_x, velocity_y, velocity_z, pressure, iBlank";
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
            WriteTecplotBinary(field_struct);
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


} // namespace zaran
