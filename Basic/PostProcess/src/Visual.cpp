#include "Visual.h"
#include "GlobalData.h"
#include "Log.h"
#include <TECIO.h>
#include <fstream>
#include <string>
#include <vector>
#include "FieldNS.h"
using namespace zaran;
void zaran::Visual::WriteTecplotBinary(Field *field)
{
    FieldNS_FNFDM *fieldNS = dynamic_cast<FieldNS_FNFDM *>(field);
    auto data_manager = fieldNS->GetDataManager();
    auto grid = fieldNS->GetGrid();
    auto cell = grid->GetCellTopo();
    auto node = grid->GetNodeTopo();

    const double *density, *velocity_x, *velocity_y, *velocity_z, *pressure;
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
    INTEGER4 file_format = 0;
    INTEGER4 debug = 0;
    INTEGER4 vIsDouble = 1;
    INTEGER4 fileType = 0;
    string grid_name = "grid_" + grid->GetName();
    string var_name = "x y z rho u v w p";
    std::string file_name = "result/" + std::to_string(GlobalData::GetInt("currentIter")) + ".plt";
    int i = TECINI142(grid_name.c_str(), var_name.c_str(), file_name.c_str(), (char *)".", &file_format, &fileType,
                      &debug, &vIsDouble);
    string zone_name = "grid_" + grid->GetName() + "_zone";
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
    int valueLocation[] = {1, 1, 1, 1, 1, 1, 1, 1};
    int shrConn = 0;
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
    INTEGER4 connectivityCount = cell_num * 8;
    std::vector<INTEGER4> cell_nodes(connectivityCount);
    for (int iCell = 0; iCell < cell_num; ++iCell)
    {
		int node_num=cell->GetNodeNum(iCell);
		auto cell2node = cell->GetNode(iCell);
        for (int iNode = 0; iNode < node_num; ++iNode)
        {
            cell_nodes[iCell * 8 + iNode] = cell2node[iNode] + 1;
        }
    }

    i = TECNODE142(&connectivityCount, cell_nodes.data());

    /// bound face
    auto face_topo = grid->GetFaceTopo();
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
    int node_num_per_cell = 4;
    connectivityCount = cell_num * node_num_per_cell;
    Array<INTEGER4> face_nodes(connectivityCount);
    for (int iFace = 0; iFace < cell_num; ++iFace)
    {
        int *face2node = face_topo->GetFace2Node(iFace);
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
void Visual::WriteVTK(Field *field)
{
    // TO DO
}
