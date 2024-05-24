#include<fstream>
#include"simulation.h"
#include "GlobalData.h"
#include "Log.h"
#include "FieldBuilder.h"
#include"ReadSTL.h"
#include"PolyData.h"
#include"File.h"
using namespace zaran;
void Simulation::Run()
{
	ReadGlobalData();
	InitSimulationTask();
	ShowInfo();
	if (task_ == SimulationTask::SOLVE_FIELD)
	{
		SolveField();
	}
	else if (task_ == SimulationTask::CONVERT_GRID)
	{
		ConvertGrid();
	}
	else if (task_ == SimulationTask::READ_MODEL)
	{
		ReadModel();
	}
	else
	{
		Log::warn("Unsupported Simulation Task! Please Check!");
		system("pause");
	}
}

void Simulation::ReadGlobalData()
{
	if (!IsFileExist(globalDataFileName_))
	{
		Log::warn("Control File:{}, is NOT exist! Please Check!", globalDataFileName_);
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
		//comments
		if (line[0] == '!' || line[0] == '#' || line[0] == '/')
			continue;
		//delete space
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
			Log::warn("Unsupported Data Type:{}, Name:{}, Value:{}", dataType, dataName, dataValue);
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
	else if (simuTask == "READ_MODEL")
	{
		task_ = SimulationTask::READ_MODEL;
	}
	else
	{
		Log::warn("Unsupported Simulation Task:{}! Please Check!", simuTask);
		system("pause");
	}
}

void Simulation::ShowInfo()
{
	Log::info("********Zaran: A Totally Automatic CFD Software!********");
	Log::info(">>>>>>>>Software Version: {}<<<<<<<<", programVersion_);
	Log::info(">>>>>>>>Software Last Modify Date: {}<<<<<<<<", lastModifyDate_);
	Log::info(">>>>>>>>Control File Version: {}<<<<<<<<", GlobalData::GetString("version"));
	if (minSupportCtrlFileVersion_ > GlobalData::GetString("version"))
	{
		Log::warn(">>>>>>>>Control File is too old, please use new Control File!<<<<<<<<");
		Log::warn(">>>>>>>>The minus version Control File is:{}<<<<<<<<", minSupportCtrlFileVersion_);
		system("pause");
	}
	Log::info(">>>>>>>>Simulation Task: {}<<<<<<<<", GlobalData::GetString("simulationTask"));
}

void zaran::Simulation::SolveField()
{
	std::string reslut_folder = GlobalData::GetString("resultFolder");
	std::string backup_folder = GlobalData::GetString("backupFieldFolder");
	CreateFolder(reslut_folder);
	CreateFolder(backup_folder);
	string grid_type_name = GlobalData::GetString("gridType");
	string solver_type_name = GlobalData::GetString("solverType");
	GridType grid_type;
	FieldSolverType solver_type;
	Dimension dim;
	if (GlobalData::GetInt("Dimension") == 2)
		dim = Dimension::two;
	else if (GlobalData::GetInt("Dimension") == 3)
		dim = Dimension::three;
	else
	{
		Log::warn("Unsupported Dimension! Please Check!");
		system("pause");
	}
	if (solver_type_name == "NS")
	{
		if (grid_type_name == "Structured")
		{
			grid_type = GridType::Structured;
			solver_type = FieldSolverType::NS_Struct;
		}
		else if (grid_type_name == "Flexible")
		{
			grid_type = GridType::Flexible;

				solver_type = FieldSolverType::NS_FNFDM;
		}
		else if (grid_type_name == "Zaran")
		{
			// grid_type = GridType::Zaran;
			// solver_type = FieldSolverType::NS_ZaRan_3D;
		}
		else
		{
			Log::warn("Unsupported Grid Type! Please Check!");
			system("pause");
		}
	}
	else
	{
		Log::warn("Unsupported Solver Type! Please Check!");
		system("pause");
	}
	FieldBuilder* fieldFactory = new FieldBuilder(grid_type, solver_type, dim);
	fieldFactory->Create();
	Field** field = fieldFactory->GetField();
	Controller* controller = new Controller(field, 1);
	controller->SolveField();
}

void zaran::Simulation::ConvertGrid()
{
	Log::info("Convert Grid!");
}

void zaran::Simulation::ReadModel()
{
	string modelFileName = GlobalData::GetString("modelFileName");
	STLReader reader;
	reader.ReadSTLFile(modelFileName.c_str());
	PolyDataModel model;
	model.SetPolyData(reader.GetMesh(), 1e-6);
	model.ShowModel();
	if (model.IsClosed())
		Log::info("Import Model: {}, is closed!", modelFileName);
	else
		Log::info("Import Model: {}, is not closed!", modelFileName);
}

