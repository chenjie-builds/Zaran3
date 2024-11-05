#include "Visual.h"
#include "GlobalData.h"
#include "Log.h"
#include <TECIO.h>
#include <fstream>
#include <string>
#include <vector>
#include "NSFieldFN.h"
using namespace zaran;
void zaran::Visual::WriteTecplotBinary(NSFieldFNFDM* field)
{
    auto data_manager = field->GetDataManager();
    auto grid = field->GetGrid();
    auto cell = grid->GetCell();
    auto node = grid->GetNode();

    const double* density, * velocity_x, * velocity_y, * velocity_z, * pressure;
    density = data_manager->GetDensity();
    velocity_x = data_manager->GetVelocity(0);
    velocity_y = data_manager->GetVelocity(1);
    velocity_z = data_manager->GetVelocity(2);
    pressure = data_manager->GetPressure();

    INTEGER4 node_num = grid->GetTotalNodeNum();
    INTEGER4 cell_num = cell->GetCellNum();
    DArray x(node_num), y(node_num), z(node_num);
    for (int iNode = 0; iNode < node_num; ++iNode)
    {
        x[iNode] = node->GetCoord(iNode)[0];
        y[iNode] = node->GetCoord(iNode)[1];
        z[iNode] = node->GetCoord(iNode)[2];
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
    std::vector<INTEGER4> cell_nodes(connectivityCount);
    for (int iCell = 0; iCell < cell_num; ++iCell)
    {
        int node_num = cell->GetNodeNum(iCell);
        auto cell2node = cell->GetNode(iCell);
        for (int iNode = 0; iNode < node_num; ++iNode)
        {
            cell_nodes[iCell * 8 + iNode] = cell2node[iNode] + 1;
        }
    }

    i = TECNODE142(&connectivityCount, cell_nodes.data());

    /// bound face
    auto face_topo = grid->GetFace();
    cell_num = face_topo->GetFaceNum();
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
    Array<INTEGER4> face_nodes(connectivityCount);
    for (int iFace = 0; iFace < cell_num; ++iFace)
    {
        int* face2node = face_topo->GetFace2Node(iFace);
        int n_node = face_topo->GetFaceNodeNum(iFace);
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

    i = TECEND142();
}

void zaran::Visual::WriteTecplotASCII(NSFieldStruct* field, std::ostream& os)
{
    auto grid = field->GetGrid();
    auto node = grid->GetNode();
    auto data_manager = field->GetDataManager();
    auto solver = field->GetSolver();
    auto metrics = solver->GetNodeMetrics();
    int grid_ni = grid->GetNi();
    int grid_nj = grid->GetNj();
    int grid_nk = grid->GetNk();
    StructIdxProxy* idx_proxy = new StructIdxProxy(grid_ni, grid_nj, grid_nk);
    int is, ie, js, je, ks, ke;
    grid->GetRange(is, ie, js, je, ks, ke);
    int ni = ie - is + 1;
    int nj = je - js + 1;
    int nk = ke - ks + 1;
    int node_num = ni * nj * nk;
    const double* density = data_manager->GetDensity();
    const double* velocity_x = data_manager->GetVelocity(0);
    const double* velocity_y = data_manager->GetVelocity(1);
    const double* velocity_z = data_manager->GetVelocity(2);
    const double* pressure = data_manager->GetPressure();

    os << "TITLE=\"Flow Field\"\n";
    os << "VARIABLES=\"X\",\"Y\",\"Z\",\"Density\",\"Velocity_x\",\"Velocity_y\",\"Velocity_z\",\"Pressure\"\n";
    os << "ZONE I=" << ni << ", J=" << nj << ", K=" << nk << ", F=POINT\n";
    for (int k = ks; k <= ke; ++k)
    {
        for (int j = js; j <= je; ++j)
        {
            for (int i = is; i <= ie; ++i)
            {
                int idx = idx_proxy->GetIdx(i, j, k);
                os << node->GetCoord(i, j, k)[0] << " " << node->GetCoord(i, j, k)[1] << " " << node->GetCoord(i, j, k)[2] << " " << density[idx] << " " << velocity_x[idx] << " " << velocity_y[idx] << " " << velocity_z[idx] << " " << pressure[idx] << "\n";
            }
        }
    }
    delete idx_proxy;
}

void zaran::Visual::WriteTecplotASCII(NSFieldZaran* field, std::ostream& os)
{
    auto grid = field->GetGrid();
    auto node = grid->GetNode();
    auto data_manager = field->GetDataManager();
    auto solver = field->GetSolver();
    int grid_ni = grid->GetNi();
    int grid_nj = grid->GetNj();
    int grid_nk = grid->GetNk();
    StructIdxProxy* idx_proxy = new StructIdxProxy(grid_ni, grid_nj, grid_nk);
    int is, ie, js, je, ks, ke;
    grid->GetRange(is, ie, js, je, ks, ke);
    int ni = ie - is + 1;
    int nj = je - js + 1;
    int nk = ke - ks + 1;
    int node_num = ni * nj * nk;
    const double* density = data_manager->GetDensity();
    const double* velocity_x = data_manager->GetVelocity(0);
    const double* velocity_y = data_manager->GetVelocity(1);
    const double* velocity_z = data_manager->GetVelocity(2);
    const double* pressure = data_manager->GetPressure();
    // 以ASCII格式写入，后期可以改为二进制格式
    os << "TITLE=\"Flow Field\"\n";
    os << "VARIABLES=\"X\",\"Y\",\"Z\",\"Density\",\"Velocity_x\",\"Velocity_y\",\"Velocity_z\",\"Pressure\",\"iBlank\"\n";
    os << "ZONE I=" << ni << ", J=" << nj << ", K=" << nk << ", F=POINT\n";
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
    delete idx_proxy;
}
void zaran::Visual::WriteTecplotBinary(NSFieldZaran* field)
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
    auto idx_proxy = new StructIdxProxy(grid_ni, grid_nj, grid_nk);
    INTEGER4 node_num = ni * nj * nk;
    INTEGER4 cell_num = (ni - 1) * (nj - 1) * (nk - 1);
    DArray x(node_num), y(node_num), z(node_num), density(node_num), velocity_x(node_num), velocity_y(node_num), velocity_z(node_num), pressure(node_num);

    std::vector<int> iBlank(node_num);
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
                density[idx] = data_manager->GetDensity(idx0);
                velocity_x[idx] = data_manager->GetVelocity(0, idx0);
                velocity_y[idx] = data_manager->GetVelocity(1, idx0);
                velocity_z[idx] = data_manager->GetVelocity(2, idx0);
                pressure[idx] = data_manager->GetPressure(idx0);
                iBlank[idx] = (int)grid->GetIBlank(i + is, j + js, k + ks);
            }
        }
    }
    INTEGER4 vIsDouble = 1;
    INTEGER4 vIsInt = 0;
    string zone_name = grid->GetName() + std::to_string(field->GetIdx());
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
    int valueLocation[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    int shrConn = 0;
    int i = TECZNE142((char*)zone_name.c_str(), &zone_type, &node_num, &cell_num, &face_num, &iCellMax, &jCellMax,
        &kCellMax, &solution_time, &strandID, &parentZn, &isBlock, &nFConns, &FNMode, 0, 0, 0, NULL,
        valueLocation, NULL, &shrConn);

    i = TECDAT142(&node_num, x.data(), &vIsDouble);
    i = TECDAT142(&node_num, y.data(), &vIsDouble);
    i = TECDAT142(&node_num, z.data(), &vIsDouble);
    i = TECDAT142(&node_num, density.data(), &vIsDouble);
    i = TECDAT142(&node_num, velocity_x.data(), &vIsDouble);
    i = TECDAT142(&node_num, velocity_y.data(), &vIsDouble);
    i = TECDAT142(&node_num, velocity_z.data(), &vIsDouble);
    i = TECDAT142(&node_num, pressure.data(), &vIsDouble);
    i = TECDAT142(&node_num, iBlank.data(), &vIsInt);
    INTEGER4 connectivityCount = cell_num * 8;
    std::vector<INTEGER4> cell_nodes(connectivityCount);
    for (int k = 0; k < nk - 1; k++)
    {
        for (int j = 0; j < nj - 1; j++)
        {
            for (int i = 0; i < ni - 1; i++)
            {
                int idx = i + (ni - 1) * j + (ni - 1) * (nj - 1) * k;
                int idx0 = i + ni * j + ni * nj * k;
                cell_nodes[idx * 8] = idx0 + 1;
                cell_nodes[idx * 8 + 1] = idx0 + 2;
                cell_nodes[idx * 8 + 2] = idx0 + ni + 2;
                cell_nodes[idx * 8 + 3] = idx0 + ni + 1;
                cell_nodes[idx * 8 + 4] = idx0 + ni * nj + 1;
                cell_nodes[idx * 8 + 5] = idx0 + ni * nj + 2;
                cell_nodes[idx * 8 + 6] = idx0 + ni * nj + ni + 2;
                cell_nodes[idx * 8 + 7] = idx0 + ni * nj + ni + 1;
            }
        }
    }
    i = TECNODE142(&connectivityCount, cell_nodes.data());
    delete idx_proxy;
    i = TECEND142();
}
void zaran::Visual::WriteTecplotBinary(NSFieldStruct* field)
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
    auto idx_proxy = new StructIdxProxy(grid_ni, grid_nj, grid_nk);
    INTEGER4 node_num = ni * nj * nk;
    INTEGER4 cell_num = (ni - 1) * (nj - 1) * (nk - 1);
    DArray x(node_num), y(node_num), z(node_num), density(node_num), velocity_x(node_num), velocity_y(node_num), velocity_z(node_num), pressure(node_num);
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
                density[idx] = data_manager->GetDensity(idx0);
                velocity_x[idx] = data_manager->GetVelocity(0, idx0);
                velocity_y[idx] = data_manager->GetVelocity(1, idx0);
                velocity_z[idx] = data_manager->GetVelocity(2, idx0);
                pressure[idx] = data_manager->GetPressure(idx0);
            }
        }
    }
    INTEGER4 file_format = 0;
    INTEGER4 debug = 0;
    INTEGER4 vIsDouble = 1;
    INTEGER4 vIsInt = 0;
    INTEGER4 fileType = 0;
    string grid_name = "grid_" + grid->GetName();
    string var_name = "x, y, z, density, velocity_x, velocity_y, velocity_z, pressure";
    std::string file_name = "result/" + std::to_string(GlobalData::GetInt("currentIter")) + "-block.plt";
    int i = TECINI142(grid_name.c_str(), var_name.c_str(), file_name.c_str(), (char*)".", &file_format, &fileType,
        &debug, &vIsDouble);
    string zone_name = grid->GetName() + std::to_string(field->GetIdx());
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
    int valueLocation[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    int shrConn = 0;
    i = TECZNE142((char*)zone_name.c_str(), &zone_type, &node_num, &cell_num, &face_num, &iCellMax, &jCellMax,
        &kCellMax, &solution_time, &strandID, &parentZn, &isBlock, &nFConns, &FNMode, 0, 0, 0, NULL,
        valueLocation, NULL, &shrConn);

    i = TECDAT142(&node_num, x.data(), &vIsDouble);
    i = TECDAT142(&node_num, y.data(), &vIsDouble);
    i = TECDAT142(&node_num, z.data(), &vIsDouble);
    i = TECDAT142(&node_num, density.data(), &vIsDouble);
    i = TECDAT142(&node_num, velocity_x.data(), &vIsDouble);
    i = TECDAT142(&node_num, velocity_y.data(), &vIsDouble);
    i = TECDAT142(&node_num, velocity_z.data(), &vIsDouble);
    i = TECDAT142(&node_num, pressure.data(), &vIsDouble);
    INTEGER4 connectivityCount = cell_num * 8;
    std::vector<INTEGER4> cell_nodes(connectivityCount);
    for (int k = 0; k < nk - 1; k++)
    {
        for (int j = 0; j < nj - 1; j++)
        {
            for (int i = 0; i < ni - 1; i++)
            {
                int idx = i + (ni - 1) * j + (ni - 1) * (nj - 1) * k;
                int idx0 = i + ni * j + ni * nj * k;
                cell_nodes[idx * 8] = idx0 + 1;
                cell_nodes[idx * 8 + 1] = idx0 + 2;
                cell_nodes[idx * 8 + 2] = idx0 + ni + 2;
                cell_nodes[idx * 8 + 3] = idx0 + ni + 1;
                cell_nodes[idx * 8 + 4] = idx0 + ni * nj + 1;
                cell_nodes[idx * 8 + 5] = idx0 + ni * nj + 2;
                cell_nodes[idx * 8 + 6] = idx0 + ni * nj + ni + 2;
                cell_nodes[idx * 8 + 7] = idx0 + ni * nj + ni + 1;
            }
        }
    }
    i = TECNODE142(&connectivityCount, cell_nodes.data());
    delete idx_proxy;
}
void Visual::WriteTecASCII(FieldManager* field_manager)
{
    std::ofstream out("result/" + std::to_string(GlobalData::GetInt("currentIter")) + ".dat");
    out << "VARIABLES=\"X\",\"Y\",\"Z\",\"Density\",\"Velocity_x\",\"Velocity_y\",\"Velocity_z\",\"Pressure\"\n";
    out << "TITLE=\"Flow Field\"\n";
    double solution_time = GlobalData::GetDouble("currentTime");
    out << "SOLUTIONTIME=" << solution_time << "\n";
    for (size_t iter_field = 0; iter_field < field_manager->GetFieldNum(); iter_field++)
    {
        auto field = field_manager->GetField(iter_field);
        auto field_type = field->GetFieldType();
        if (field_type == FieldType::NS_Structured)
        {
            auto field_struct = dynamic_cast<NSFieldStruct*>(field);
            WriteTecplotASCII(field_struct, out);
        }
        else if (field_type == FieldType::NS_Zaran)
        {
            auto field_zaran = dynamic_cast<NSFieldZaran*>(field);
            WriteTecplotASCII(field_zaran, out);
        }
        else
        {
            Log::warn("Field type is not supported!");
        }

    }

}

void zaran::Visual::WriteTecplotBinary(FieldManager* field_manager)
{
    INTEGER4 file_format = 0;
    INTEGER4 debug = 0;
    INTEGER4 vIsDouble = 1;
    INTEGER4 vIsInt = 0;
    INTEGER4 fileType = 0;
    string grid_name = "grid";
    string var_name = "x, y, z, density, velocity_x, velocity_y, velocity_z, pressure";
    std::string file_name = "result/" + std::to_string(GlobalData::GetInt("currentIter")) + ".plt";
    int i = TECINI142(grid_name.c_str(), var_name.c_str(), file_name.c_str(), (char*)".", &file_format, &fileType,
        &debug, &vIsDouble);
    for (size_t iter_field = 0; iter_field < field_manager->GetFieldNum(); iter_field++)
    {
        auto field = field_manager->GetField(iter_field);
        auto field_type = field->GetFieldType();
        if (field_type == FieldType::NS_Zaran)
        {
            auto field_zaran = dynamic_cast<NSFieldZaran*>(field);
            WriteTecplotBinary(field_zaran);
        }
        else if (field_type == FieldType::NS_FlexibleNode)
        {
            auto field_fn = dynamic_cast<NSFieldFNFDM*>(field);
            WriteTecplotBinary(field_fn);
        }
        else if (field_type == FieldType::NS_Structured)
        {
            auto field_struct = dynamic_cast<NSFieldStruct*>(field);
            WriteTecplotBinary(field_struct);
        }
        else
        {
            Log::warn("Field type is not supported!");
        }
    }
    i = TECEND142();
}
