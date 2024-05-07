#include "controller.h"
#include "log.h"
#include <fstream>
#include"FlowSolver.h"
#include"MathBasic.h"
#include "SpecialField.h"
#include"File.h"
using namespace zaran;



zaran::Controller::Controller(Field** field, int field_size)
{
    m_field = field;
    m_field_size = field_size;
}

Controller::~Controller()
{
}

void Controller::Initialize()
{
    if (!GlobalData::IsExist("currentIter"))
    {
        GlobalData::Update("currentIter", 0);
    }
    if (!GlobalData::IsExist("currentTime"))
    {
        double startTime = GlobalData::GetDouble("startTime");
        GlobalData::Update("currentTime", startTime);
    }
    for (size_t iField = 0; iField < m_field_size; iField++)
    {
        m_field[iField]->GetSolver()->Init();
    }
}
void Controller::SaveWallNode()
{
    std::ofstream fout("boundNode.dat");
    fout << "variables=x,y,p" << std::endl;
    for (size_t iField = 0; iField < m_field_size; iField++)
    {
        Grid* currentGrid = m_field[iField]->GetGrid();
        NodeTopo* nodeTopo = currentGrid->GetNodeTopo();
        BoundaryMap* boundMap = currentGrid->GetBoundaryMap();
    }
    fout.close();
}
void Controller::SaveDataTecplot()
{
    for (size_t iField = 0; iField < m_field_size; iField++)
    {
        FieldSolver* solver = m_field[iField]->GetSolver();
        m_visual->WriteTecplot(solver);
        // m_visual->WriteTecplot2D(std::dynamic_pointer_cast<FieldSolver> (currentSolver));
        //  m_visual->WriteTecplotZaran3D(std::dynamic_pointer_cast<FieldSolver> (currentSolver));
        // m_visual->WriteTecplotZaran3DBinary(std::dynamic_pointer_cast<FieldSolver> (currentSolver));
    }
}
void Controller::SaveDataVTK(std::ostream& os)
{
    /*os << "# vtk DataFile Version 4.2\n";
    for (size_t iGrid = 0; iGrid < grid_.size(); iGrid++)
    {
        auto& currentGrid = grid_[iGrid];
        os << "grid" << iGrid << "\n"
            << "ASCII\n"
            << "DATASET UNSTRUCTURED_GRID\n"
            << "POINTS " << currentGrid->GetNodeNum() << " double\n";
        for (size_t iNode = 0; iNode < currentGrid->GetNodeNum(); iNode++)
        {
            auto& currentNode = currentGrid->GetNode(iNode);
            auto& currentNodeCoord = currentNode.GetCoord();
            os << currentNodeCoord.x() << "  " << currentNodeCoord.y() << " " << currentNodeCoord.z() << "\n";
        }
        size_t elementNodeNum = 0;
        for (size_t iElem = 0; iElem < currentGrid->GetElementNum(); ++iElem)
        {
            elementNodeNum += currentGrid->GetElement(iElem).GetNode().size();
        }
        os << "CELLS " << currentGrid->GetElementNum() << " " << currentGrid->GetElementNum() + elementNodeNum << "\n";
        for (size_t iElem = 0; iElem < currentGrid->GetElementNum(); ++iElem)
        {
            auto& currentElement = currentGrid->GetElement(iElem);
            os << currentElement.GetNode().size() << " ";
            for (auto& iNode : currentElement.GetNode())
            {
                os << iNode->GetIndex() << "  ";
            }
            os << "\n";
        }
        os << "CELL_TYPES " << currentGrid->GetElementNum() << "\n";
        for (size_t iElem = 0; iElem < currentGrid->GetElementNum(); ++iElem)
        {
            os << 7 << "\n";
        }
        os << "\n";
    }*/
}

void Controller::SolveField()
{
    Log::info("Start to solve field!");
    Initialize();
    Log::info("Initialize finished!");
    Log::info("Save init field data...");
    SaveFieldData();
    Log::info("Save init field data finished!");
    SaveResidual();
    while (true)
    {
        PreSolve();
        SolveFieldOneStep();
        PostSolve();
        if (IsStopSolve())
            break;
    }
    PostSolve();
}

double Controller::CalcMaxAveResidual()
{
    double maxResidual = 0.0;
    for (size_t iField = 0; iField < m_field_size; iField++)
    {
        FlowSolver* solver = dynamic_cast<FlowSolver*>(m_field[iField]->GetSolver());
        maxResidual = Max(maxResidual, solver->ComputeMaxResidual());
    }
    return maxResidual;
}

void zaran::Controller::CalcResidual()
{
    for (size_t iField = 0; iField < m_field_size; iField++)
    {
        FlowSolver* solver = dynamic_cast<FlowSolver*>(m_field[iField]->GetSolver());
        Grid* grid = m_field[iField]->GetGrid();
        FieldData* fieldData = m_field[iField]->GetFieldData();
        int n_data;
        fieldData->GetDataSize("density", n_data);
        NodeTopo* nodeTopo = grid->GetNodeTopo();

        auto& nodeCoord = nodeTopo->GetCoordinate();
        auto& nodeType = nodeTopo->GetType();
        double mar_res = 0.0;
        double ave_res = 0.0;
#pragma omp parallel for reduction(max:mar_res) reduction(+:ave_res)
        for (int iNode = 0;iNode < n_data;++iNode)
        {
            mar_res = Max(mar_res, abs(solver->GetResidual(iNode, 0)));
            ave_res += pow(solver->GetResidual(iNode, 0), 2);
        }
        ave_res /= n_data;
        ave_res = sqrt(ave_res);
        m_max_res = mar_res;
        m_ave_res = ave_res;
    }
    m_res_flag = true;
}

void Controller::SaveFieldData()
{
    SaveDataTecplot();
    std::string back_dir = GlobalData::GetString("backupFieldFolder");
    int currentIter = GlobalData::GetInt("currentIter");
    back_dir += "/iter=" + std::to_string(currentIter);
    CreateFolder(back_dir);
    BackupFieldData(back_dir);
    BackupResidual(back_dir);
    BackupLog(back_dir);
    BackupGlobalData(back_dir);
}

void Controller::BackupFieldData(std::string& back_folder)
{
    for (size_t iField = 0; iField < m_field_size; iField++)
    {
        FieldSolver* solver = m_field[iField]->GetSolver();
        solver->BackupField(back_folder);
    }
}

void Controller::BackupResidual(std::string& back_folder)
{
    std::string residual_file = GlobalData::GetString("residualFileName");
    std::string residual_file_back = back_folder + "/" + residual_file;
    if (IsFileExist(residual_file) == false)
    {
        Log::warn("residual file:{} is not exist!", residual_file);
        return;
    }
    if (IsFileExist(residual_file_back) == true)
    {
        DeleteFile(residual_file_back);
    }
    CopyFile(residual_file, residual_file_back);
}

void Controller::BackupLog(std::string& back_folder)
{
    std::string log_file = "log.txt";
    std::string log_file_back = back_folder + "/" + log_file;
    if (IsFileExist(log_file) == false)
    {
        Log::warn("logFile:{} is not exist!", log_file);
        return;
    }
    if (IsFileExist(log_file_back) == true)
    {
        DeleteFile(log_file_back);
    }
    CopyFile(log_file, log_file_back);
}

void zaran::Controller::BackupGlobalData(std::string& back_folder)
{
    std::string global_file = "zaran.ini";
    std::string global_file_back = back_folder + "/" + global_file;
    if (IsFileExist(global_file) == false)
    {
        Log::warn("gloabal data file:{} is not exist!", global_file);
        return;
    }
    if (IsFileExist(global_file_back) == true)
    {
        DeleteFile(global_file_back);
    }
    GlobalData::Backup(back_folder);
}

bool Controller::IsStopSolve()
{
    double endTime = GlobalData::GetDouble("endTime");
    int currentIter = GlobalData::GetInt("currentIter");
    int calResidualIter = GlobalData::GetInt("calResidualIter");
    int maxIter = GlobalData::GetInt("maxIter");
    double minResidual = GlobalData::GetDouble("minResidual");
    double currentTime = GlobalData::GetDouble("currentTime");
    if (m_res_flag == true)
    {
        if (m_max_res < minResidual)
        {
            Log::info("Max Residual is small than {}, stop compute!", minResidual);
            return true;
        }
    }
    if (currentIter > maxIter)
    {
        Log::info("Max Iter={}, stop compute!", maxIter);
        return true;
    }
    if (currentTime > endTime || abs(currentTime - endTime) < SMALL_NUMBER)
    {
        Log::info("Max time={}, stop compute!", endTime);
        return true;
    }
    return false;
}
void Controller::SaveResidual()
{
    int currentIter = GlobalData::GetInt("currentIter");
    string residual_file = GlobalData::GetString("residualFileName");


    if (currentIter == 0)
    {
        std::ofstream fout(residual_file);
        fout << "variables=step, time, MaxRes, AveRes\n";
        fout.close();
    }
    else
    {
        std::ofstream fout(residual_file, std::ios::app);
        fout << currentIter << "\t\t" << GlobalData::GetDouble("currentTime") << "\t\t" << m_max_res << "\t\t" << m_ave_res << std::endl;
        fout.close();
    }
}

void Controller::SolveFieldOneStep()
{
    for (size_t iField = 0; iField < m_field_size; iField++)
    {
        FieldSolver* solver = m_field[iField]->GetSolver();
        solver->Solve();
    }
}

void Controller::PreSolve()
{
    int currentIter = GlobalData::GetInt("currentIter");
    GlobalData::Update("currentIter", ++currentIter);
    for (size_t iField = 0; iField < m_field_size; iField++)
    {
        FieldSolver* solver = m_field[iField]->GetSolver();
        solver->Preprocess();
    }
}

void Controller::PostSolve()
{
    for (size_t iField = 0; iField < m_field_size; iField++)
    {
        FieldSolver* solver = m_field[iField]->GetSolver();
        solver->Postprocess();
    }
    CommInterNodeData();
    int currentIter = GlobalData::GetInt("currentIter");
    int calResidualIter = GlobalData::GetInt("calResidualIter");
    int writeFieldIter = GlobalData::GetInt("writeFieldIter");
    if (currentIter % calResidualIter == 0 || IsStopSolve())
    {
        CalcResidual();
        Log::info("currentIter={}, dt={:e}, max_res={:e}, ave_res={:e}", GlobalData::GetInt("currentIter"), GlobalData::GetDouble("dt"), m_max_res, m_ave_res);
        SaveResidual();
    }
    if (currentIter % writeFieldIter == 0 || IsStopSolve())
    {
        SaveFieldData();
        // SaveWallNode();
    }
}

void Controller::CommInterNodeData()
{
    //for (size_t iSolver = 0; iSolver < solverVec_->GetSolverNumber(); iSolver++)
    //{
    //	solverVec_->GetSolver(iSolver)->CommInterNodeData();
    //}
}
