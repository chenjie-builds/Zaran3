#include "controller.h"
#include "log.h"
#include <fstream>
#include"FlowSolver.h"
#include"MathBasic.h"
#include "SpecialField.h"
using namespace zaran;



zaran::Controller::Controller(Array<Ptr<Field>>& field)
{
    m_field = field;
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
    for (size_t iField = 0; iField < m_field.size(); iField++)
    {
        m_field[iField]->GetSolver()->Init();
    }
}
void Controller::SaveWallNode()
{
    std::ofstream fout("boundNode.dat");
    fout << "variables=x,y,p" << std::endl;
    for (size_t iField = 0; iField < m_field.size(); iField++)
    {
        auto& currentGrid = m_field[iField]->GetGrid();
        auto& nodeTopo = currentGrid->GetNodeTopo();
        auto& boundNode = currentGrid->GetBoundaryMap();
    }
    fout.close();
}
void Controller::SaveDataTecplot()
{
    for (size_t iField = 0; iField < m_field.size(); iField++)
    {
        auto& currentSolver = m_field[iField]->GetSolver();
        m_visual->WriteTecplot(std::dynamic_pointer_cast<FieldSolver> (currentSolver));
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
    while (!IsStopSolve())
    {
        PreSolve();
        SolveFieldOneStep();
        PostSolve();
    }
    PostSolve();
}

double Controller::CalcMaxAveResidual()
{
    double maxResidual = 0.0;
    for (size_t iField = 0; iField < m_field.size(); iField++)
    {
        auto& currentSolver = std::dynamic_pointer_cast<FlowSolver>(m_field[iField]->GetSolver());
        maxResidual = Max(maxResidual, currentSolver->ComputeMaxResidual());
    }
    return maxResidual;
}

void zaran::Controller::CalcResidual()
{
    for (size_t iField = 0; iField < m_field.size(); iField++)
    {
        auto& currentSolver = std::dynamic_pointer_cast<FlowSolver>(m_field[iField]->GetSolver());
        auto& currentGrid = m_field[iField]->GetGrid();
        auto& fieldData = m_field[iField]->GetFieldData();
        double* rho, * res0;
        fieldData->GetData("rho", rho);
        fieldData->GetData("res0", res0);
        int n_data;
        fieldData->GetDataSize("rho", n_data);
        auto& nodeTopo = currentGrid->GetNodeTopo();
        auto& nodeCoord = nodeTopo->GetCoordinate();
        auto& nodeType = nodeTopo->GetType();
        double maxResidual = 0.0;
        double aveResidual = 0.0;
#pragma omp parallel for reduction(max:maxResidual) reduction(+:aveResidual)
        for (int iNode = 0;iNode < n_data;++iNode)
        {
            maxResidual = Max(maxResidual, abs(res0[iNode]));
            aveResidual += pow(res0[iNode], 2);
        }
        aveResidual /= n_data;
        aveResidual = sqrt(aveResidual);
        maxResidual_ = maxResidual;
        aveResidual_ = aveResidual;
    }
}

void Controller::SaveFieldData()
{
    SaveDataTecplot();
}

bool Controller::IsStopSolve()
{
    double endTime = GlobalData::GetDouble("endTime");
    int currentIter = GlobalData::GetInt("currentIter");
    int calResidualIter = GlobalData::GetInt("calResidualIter");
    int maxIter = GlobalData::GetInt("maxIter");
    double minResidual = GlobalData::GetDouble("minResidual");
    double currentTime = GlobalData::GetDouble("currentTime");
    //达到最大迭代次数
    if (currentIter > maxIter || maxResidual_ < minResidual && currentIter > calResidualIter)
    {
        Log::info("Max Residual is small than {}, stop compute!", minResidual);
        return true;
    }
    //达到最大计算时间
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
    if (currentIter == 0)
    {
        std::ofstream fout("res.dat");
        fout << "variables=step,time,MaxRes,AveRes\n";
        fout.close();
    }
    else
    {
        std::ofstream fout("res.dat", std::ios::app);
        fout << currentIter << "\t\t" << GlobalData::GetDouble("currentTime") << "\t\t" << maxResidual_ << "\t\t" << aveResidual_ << std::endl;
        fout.close();
    }
}

void Controller::SolveFieldOneStep()
{
    for (size_t iField = 0; iField < m_field.size(); iField++)
    {
        auto& currentSolver = m_field[iField]->GetSolver();
        currentSolver->Solve();
    }
}

void Controller::PreSolve()
{
    int currentIter = GlobalData::GetInt("currentIter");
    GlobalData::Update("currentIter", ++currentIter);
}

void Controller::PostSolve()
{

    CommInterNodeData();
    int currentIter = GlobalData::GetInt("currentIter");
    int calResidualIter = GlobalData::GetInt("calResidualIter");
    int writeFieldIter = GlobalData::GetInt("writeFieldIter");
    if (currentIter % calResidualIter == 0 || IsStopSolve())
    {
        CalcResidual();
        Log::info("currentIter={}, dt={:e}, max_res={:e}, ave_res={:e}", GlobalData::GetInt("currentIter"), GlobalData::GetDouble("dt"), maxResidual_, aveResidual_);
        if (GlobalData::IsExist("min_dt_index"))
        {
            Log::info("min_dt_index={}", GlobalData::GetInt("min_dt_index"));
        }
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
