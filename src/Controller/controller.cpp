#include "controller.h"
#include "log.h"
#include <fstream>
#include"FlowSolver.h"
#include"MathBasic.h"
using namespace zaran;
Controller::Controller(Ptr<GridList>& gridList, Ptr<SolverVec>& solverVec)
{
	gridList_ = gridList;
	solverVec_ = solverVec;
	visual_ = std::make_shared<Visual>();
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
	for (size_t iSolver = 0; iSolver < solverVec_->GetSolverNumber(); iSolver++)
	{
		solverVec_->GetSolverPtr(iSolver)->Init();
	}
	SaveFieldData();
}
void Controller::SaveWallNode()
{
	std::ofstream fout("boundNode.dat");
	fout << "variables=x,y,p" << std::endl;
	for (size_t iGrid = 0; iGrid < gridList_->GetGridNumber(); iGrid++)
	{
		auto& currentGrid = gridList_->GetGrid(iGrid);
		auto& nodeTopo = currentGrid->GetNodeTopoInfo();
		auto& boundNode = currentGrid->GetBoundaryMap();
	}
	fout.close();
}
void Controller::SaveDataTecplot()
{
	for (int iSolver = 0; iSolver < solverVec_->GetSolverNumber(); iSolver++)
	{
		auto& currentSolver = solverVec_->GetSolverPtr(iSolver);
		visual_->WriteTecplot(std::dynamic_pointer_cast<FieldSolver> (currentSolver));
		//visual_->WriteTecplotPoint(std::dynamic_pointer_cast<FieldSolver> (currentSolver));
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
}

double Controller::CalcMaxAveResidual()
{
	double maxResidual = 0.0;
	for (size_t iSolver = 0; iSolver < solverVec_->GetSolverNumber(); iSolver++)
	{
		auto& currentSolver = std::dynamic_pointer_cast<FlowSolver>(solverVec_->GetSolverPtr(iSolver));
		maxResidual = Max(maxResidual, currentSolver->ComputeMaxResidual());
	}
	return maxResidual;
}

void Controller::SaveFieldData()
{
	SaveDataTecplot();
}

bool Controller::IsStopSolve()
{
	//double endTime = GlobalData::GetDouble("endTime");
	//double currentTime = physicTime_->GetCurrentPhysicsTime();
	//int iterStep = GlobalData::GetInt("step");
	//int calResidualStep = GlobalData::GetInt("calResidualStep");
	//double minResidual = GlobalData::GetDouble("minResidual");
	////达到要求的最小残差
	//if (iterStep > calResidualStep && maxResidual_ < minResidual)
	//{
	//	ZaranLog::info("Max Residual is small than {}, stop compute!", minResidual);
	//	return true;
	//}
	////达到最大计算时间
	//if (currentTime > endTime || abs(currentTime - endTime) < SMALL_NUMBER)
	//{
	//	ZaranLog::info("Max time={}, stop compute!", endTime);
	//	return true;
	//}
	return false;
}
void Controller::SaveResidual()
{
	/*int step = GlobalData::GetInt("step");
	if (step == 0)
	{
		std::ofstream fout("res.dat");
		fout << "variables=step,time,MaxRes,AveRes\n";
		fout.close();
	}
	else
	{
		std::ofstream fout("res.dat", std::ios::app);
		fout << step << "\t\t" << physicTime_->GetCurrentPhysicsTime() << "\t\t" << maxResidual_ << "\t\t" << aveResidual_ << endl;
		fout.close();
	}*/
}

void Controller::SolveFieldOneStep()
{
	for (size_t iSolver = 0; iSolver < solverVec_->GetSolverNumber(); iSolver++)
	{
		solverVec_->GetSolverPtr(iSolver)->Solve();
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
	if (iterStep % calResidualStep == 0)
	{
		double maxResidual = CalcMaxAveResidual();
		ZaranLog::info("step={}, dt={:e}, maxRes={:e}", GlobalData::GetInt("step"), GlobalData::GetDouble("dt"), maxResidual);
		SaveResidual();
	}
	if (iterStep % writeFieldStep == 0)
	{
		SaveFieldData();
		SaveWallNode();
	}
}

void Controller::CommInterNodeData()
{
	//for (size_t iSolver = 0; iSolver < solverVec_->GetSolverNumber(); iSolver++)
	//{
	//	solverVec_->GetSolverPtr(iSolver)->CommInterNodeData();
	//}
}
