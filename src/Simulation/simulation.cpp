#include<fstream>
#include"simulation.h"
#include "GlobalData.h"
#include "log.h"
#include "GridListFactory.h"
#include"FNFDM3D.h"
#include"FNFDM2D.h"
#include <filesystem>
using namespace zaran;
void Simulation::Start()
{
	ReadGlobalData();
	InitSimulationTask();
	ShowInfo();
	if (task_ == SimulationTask::SOLVE_FIELD)
	{
		Ptr<GridList> gridList;
		//gridListFatory_ = std::make_shared<GridListFactoryFNFDM3D>();
		gridListFatory_ = std::make_shared<GridListFactoryFNFDM2D>();
		//gridListFatory_ = std::make_shared<GridListFactory>();
		gridListFatory_->Create(gridList);
		Ptr<SolverVec> solverVec = std::make_shared<SolverVec>();
		solverFactory_ = std::make_shared<SolverFactory>();
		solverFactory_->Create(gridList, solverVec);
		Ptr<Controller> controller = std::make_shared<Controller>(gridList, solverVec);
		controller->SolveField();
	}
	else if (task_ == SimulationTask::CONVERT_GRID)
	{

	}

}

void Simulation::ReadGlobalData()
{
	if (!std::filesystem::exists(globalDataFileName_))
	{
		ZaranLog::warn("Control File:{}, is NOT exist! Please Check!", globalDataFileName_);
		system("pause");
	}
	std::ifstream fin(globalDataFileName_);
	std::string line;
	std::string dataType;
	std::string dataName;
	std::string dataValue;

	while (std::getline(fin, line))
	{
		std::string separator = " =\r\n\t#$,;\"";
		if (line.empty())
			continue;
		//注释行
		if (line[0] == '!' || line[0] == '！' || line[0] == '#' || line[0] == '/')
			continue;
		//删除前方空格
		size_t id = line.find_first_not_of(' ');
		line.erase(line.begin(), line.begin() + id);
		id = line.find_first_of(separator);
		dataType = line.substr(0, id);
		line.erase(0, id);
		line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
		if (line.empty())
			continue;
		id = line.find_first_of("=");
		if (id == std::string::npos)
			continue;
		dataName = line.substr(0, id);
		dataValue = line.substr(id + 1);
		if (dataType == "string")
			GlobalData::Update(dataName, dataValue);
		else if (dataType == "double")
			GlobalData::Update(dataName, stod(dataValue));
		else if (dataType == "int")
			GlobalData::Update(dataName, stoi(dataValue));
		else
		{
			ZaranLog::warn("Unspported Data Type:{}, Name:{}, Value:{}", dataType, dataName, dataValue);
		}
	}
}


void Simulation::InitSimulationTask()
{
	std::string simuTask = GlobalData::GetString("simulationTask");
	if (simuTask == "SOLVE_FIELD")
	{
		task_ = SimulationTask::SOLVE_FIELD;
	}
	else if (simuTask == "CONVERT_GRID")
	{
		task_ = SimulationTask::CONVERT_GRID;
	}
	else
	{
		ZaranLog::warn("Unsupported Simulation Task:{}! Please Check!", simuTask);
		system("pause");
	}
}

void Simulation::ShowInfo()
{
	ZaranLog::info("********Zaran: A Totally Automatic CFD Solver!");
	ZaranLog::info(">>>>>>>>Software Version: {}", programVersion_);
	ZaranLog::info(">>>>>>>>Software Last Modify Date: {}", lastModifyDate_);
	ZaranLog::info(">>>>>>>>Control File Version: {}", GlobalData::GetString("version"));
	if (minSupportContrlFileVersion_ > GlobalData::GetString("version"))
	{
		ZaranLog::warn(">>>>>>>>Control File is too old, please use new Control File!");
		ZaranLog::warn(">>>>>>>>The mininus version Control File is:{}", minSupportContrlFileVersion_);
		system("pause");
	}
	ZaranLog::info(">>>>>>>>Simulation Task: {}", GlobalData::GetString("simulationTask"));
}

