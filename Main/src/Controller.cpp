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
    if (!GlobalData::IsExist("step"))
    {
        GlobalData::Update("step", 0);
    }
    if (!GlobalData::IsExist("globalTime"))
    {
        double startTime = GlobalData::GetDouble("startTime");
        GlobalData::Update("globalTime", startTime);
    }
    for (size_t iField = 0; iField < m_field.size(); iField++)
    {
        m_field[iField]->GetSolver()->Init();
    }
    SaveFieldData();
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
        // visual_->WriteTecplot(std::dynamic_pointer_cast<FieldSolver> (currentSolver));
        m_visual->WriteTecplot2D(std::dynamic_pointer_cast<FieldSolver> (currentSolver));
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
    ZaranLog::info("Start to solve field!");
    Initialize();
    ZaranLog::info("Initialize finished!");
    SaveFieldData();
    ZaranLog::info("Save data");
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
        auto& rho = fieldData->GetData("rho");
        auto& nodeTopo = currentGrid->GetNodeTopo();
        auto& nodeCoord = nodeTopo->GetCoordinate();
        auto& nodeType = nodeTopo->GetType();
        double x, y;
        DVector prim;
        double theory_rho;
        maxResidual_ = aveResidual_ = 0.0;
        int n_data = rho.size();
        for (int iNode = 0;iNode < n_data;++iNode)
        {
            x = nodeCoord[iNode].x();
            y = nodeCoord[iNode].y();
            CalcIsentropicVortex(x, y, 5.0, prim);
            theory_rho = prim[0];
            maxResidual_ = Max(maxResidual_, abs(rho[iNode] - theory_rho));
            aveResidual_ += pow(rho[iNode] - theory_rho, 2);
        }
        aveResidual_ /= n_data;
        aveResidual_ = sqrt(aveResidual_);
    }
}

void Controller::SaveFieldData()
{
    SaveDataTecplot();
}

bool Controller::IsStopSolve()
{
    double endTime = GlobalData::GetDouble("endTime");
    int iterStep = GlobalData::GetInt("step");
    int calResidualStep = GlobalData::GetInt("calResidualStep");
    double minResidual = GlobalData::GetDouble("minResidual");
    double currentTime = GlobalData::GetDouble("globalTime");
    //达到要求的最小残差
    if (iterStep > calResidualStep && maxResidual_ < minResidual)
    {
        ZaranLog::info("Max Residual is small than {}, stop compute!", minResidual);
        return true;
    }
    //达到最大计算时间
    if (currentTime > endTime || abs(currentTime - endTime) < SMALL_NUMBER)
    {
        ZaranLog::info("Max time={}, stop compute!", endTime);
        return true;
    }
    return false;
}
void Controller::SaveResidual()
{
    int step = GlobalData::GetInt("step");
    if (step == 0)
    {
        std::ofstream fout("res.dat");
        fout << "variables=step,time,MaxRes,AveRes\n";
        fout.close();
    }
    else
    {
        std::ofstream fout("res.dat", std::ios::app);
        fout << step << "\t\t" << GlobalData::GetDouble("globalTime") << "\t\t" << maxResidual_ << "\t\t" << aveResidual_ << std::endl;
        fout.close();
    }
}

void Controller::SolveFieldOneStep()
{
    for (size_t iField = 0; iField < m_field.size(); iField++)
    {
        auto& currentSolver =m_field[iField]->GetSolver();
        currentSolver->Solve();
    }
}

void Controller::PreSolve()
{
    int iterStep = GlobalData::GetInt("step");
    GlobalData::Update("step", ++iterStep);
}

void Controller::PostSolve()
{

    CommInterNodeData();
    int iterStep = GlobalData::GetInt("step");
    int calResidualStep = GlobalData::GetInt("calResidualStep");
    int writeFieldStep = GlobalData::GetInt("writeFieldStep");
    if (iterStep % calResidualStep == 0 || IsStopSolve())
    {
        CalcResidual();
        ZaranLog::info("step={}, dt={:e}, maxRes={:e}, aveRes={:e}", GlobalData::GetInt("step"), GlobalData::GetDouble("dt"), maxResidual_, aveResidual_);
        SaveResidual();
    }
    if (iterStep % writeFieldStep == 0 || IsStopSolve())
    {
        SaveFieldData();
        SaveWallNode();
    }
}

void Controller::CommInterNodeData()
{
    //for (size_t iSolver = 0; iSolver < solverVec_->GetSolverNumber(); iSolver++)
    //{
    //	solverVec_->GetSolver(iSolver)->CommInterNodeData();
    //}
}
