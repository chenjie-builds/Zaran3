#include<fstream>
#include"Application.h"
#include "GlobalData.h"
#include"GlobalField.h"
#include "Log.h"
#include "FieldGenerator.h"
#include"ReadSTL.h"
#include"PolyData.h"
#include"File.h"
#include "GridConvert.h"
#include "ZaranError.h"
#include "DEMFieldGenerator.h"
#include "DEMFieldSimulation.h"
#include <sstream>

namespace
{
	// 解析 "V主.次.修订" 或 "主.次.修订" 为数值三元组
	bool ParseVersion(const std::string& version, int& major, int& minor, int& patch)
	{
		if (version.empty())
			return false;
		size_t start = (version[0] == 'V' || version[0] == 'v') ? 1 : 0;
		char dot1 = 0, dot2 = 0;
		std::stringstream ss(version.substr(start));
		ss >> major >> dot1 >> minor >> dot2 >> patch;
		return !ss.fail() && dot1 == '.' && dot2 == '.';
	}
}

namespace zaran
{
	void Application::Run()
	{
		ReadGlobalData();
		InitTask();
		ShowInfo();
		if (m_task == TaskType::SOLVE_FIELD)
		{
			SolveField();
		}
		else if (m_task == TaskType::CONVERT_GRID)
		{
			ConvertGrid();
		}
		else if (m_task == TaskType::READ_MODEL)
		{
			ReadModel();
		}
		else
		{
			Log::warn("Unsupported Simulation Task! Please Check!");
			throw ZaranError("Unsupported Simulation Task");
		}
	}

	void Application::ReadGlobalData()
	{
		GlobalData::Load(m_control_file);
		GlobalData::Update("work_dir", m_work_dir); // 运行时键，加载后重新写入
	}


	void Application::InitTask()
	{
		std::string simuTask = GlobalData::GetString("task.simulation");
		if (simuTask == "SOLVE_FIELD")
		{
			m_task = TaskType::SOLVE_FIELD;
		}
		else if (simuTask == "CONVERT_GRID")
		{
			m_task = TaskType::CONVERT_GRID;
		}
		else if (simuTask == "READ_MODEL")
		{
			m_task = TaskType::READ_MODEL;
		}
		else
		{
			Log::warn("Unsupported Simulation Task:{}! Please Check!", simuTask);
			throw ZaranError("Unsupported Simulation Task: " + simuTask);
		}
	}

	void Application::ShowInfo()
	{
		Log::info("********Zaran: A Totally Automatic CFD Software!********");
		Log::info(">>>>>>>>Control File Version: {}<<<<<<<<", GlobalData::GetString("version"));
		if (IsVersionTooOld(GlobalData::GetString("version")))
		{
			Log::warn(">>>>>>>>Control File is too old, please use new Control File!<<<<<<<<");
			Log::warn(">>>>>>>>The minus version Control File is:{}<<<<<<<<", m_min_ctrol_file_version);
			throw ZaranError("Control file version is too old");
		}
		Log::info(">>>>>>>>Simulation Task: {}<<<<<<<<", GlobalData::GetString("task.simulation"));
	}

	bool Application::IsVersionTooOld(const std::string& version) const
	{
		int v_major = 0, v_minor = 0, v_patch = 0;
		int m_major = 0, m_minor = 0, m_patch = 0;
		if (!ParseVersion(version, v_major, v_minor, v_patch) ||
			!ParseVersion(m_min_ctrol_file_version, m_major, m_minor, m_patch))
		{
			// 解析失败时回退到字典序比较
			return m_min_ctrol_file_version > version;
		}
		if (v_major != m_major) return v_major < m_major;
		if (v_minor != m_minor) return v_minor < m_minor;
		return v_patch < m_patch;
	}

	void Application::SolveField() const
	{
		string solver_type_name = GlobalData::GetString("task.solver");

		// DEM 分支：不需要网格类型和空间维数
		if (solver_type_name == "DEM")
		{
			if (GlobalData::IsExist("output.result_folder"))
				CreateFolder(m_work_dir + "/" + GlobalData::GetString("output.result_folder"));
			if (GlobalData::IsExist("init.backup_folder"))
				CreateFolder(m_work_dir + "/" + GlobalData::GetString("init.backup_folder"));
			shared_ptr<DEMFieldGenerator> dem_factory = make_shared<DEMFieldGenerator>();
			shared_ptr<FieldManager> global_field = dem_factory->Create();
			shared_ptr<DEMFieldSimulation> controller = make_shared<DEMFieldSimulation>(global_field);
			controller->SolveField();
			return;
		}

		std::string reslut_folder = GlobalData::GetString("output.result_folder");
		std::string backup_folder = GlobalData::GetString("init.backup_folder");
		CreateFolder(m_work_dir + "/" + reslut_folder);
		CreateFolder(m_work_dir + "/" + backup_folder);

		string grid_type_name = GlobalData::IsExist("task.grid_type") ? GlobalData::GetString("task.grid_type") : "";
		GridType grid_type = GridType::Unkown;
		FieldSolverType solver_type = FieldSolverType::NS_Struct;
		Dimension dim = Dimension::three;
		if (GlobalData::GetInt("task.dimension") == 2)
			dim = Dimension::two;
		else if (GlobalData::GetInt("task.dimension") == 3)
			dim = Dimension::three;
		else
		{
			Log::warn("Unsupported Dimension! Please Check!");
			throw ZaranError("Unsupported Dimension");
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
				grid_type = GridType::Zaran;
				solver_type = FieldSolverType::NS_ZaRan;
			}
			else
			{
				Log::warn("Unsupported Grid Type! Please Check!");
				throw ZaranError("Unsupported Grid Type");
			}
		}
		else
		{
			Log::warn("Unsupported Solver Type! Please Check!");
			throw ZaranError("Unsupported Solver Type");
		}
		shared_ptr<FieldGenerator> fieldFactory = make_shared<FieldGenerator>(grid_type, solver_type, dim);
		//FieldGeneratorBuildingExplosion* fieldFactory = new FieldGeneratorBuildingExplosion(grid_type, solver_type, dim);
		shared_ptr<FieldManager> global_field = fieldFactory->Create();
		shared_ptr<NSFieldSimulation> controller = make_shared <NSFieldSimulation>(global_field);
		controller->SolveField();
	}

	void Application::ConvertGrid()
	{
		Log::info("Convert Grid Start!");
		string grid_file_name = "ZoneTestInput.dat";
		shared_ptr<GridConvertPiflow2Flexible> gridConvert = make_shared<GridConvertPiflow2Flexible>();
		gridConvert->SetFileName(grid_file_name);
		gridConvert->ReadFile();
		gridConvert->WriteFile();
	}

	void Application::ReadModel()
	{
		string modelFileName = GlobalData::GetString("zaran.model_file");
		STLReader reader;
		reader.ReadSTLFile(modelFileName.c_str());
		PolyDataModel model;
		model.SetPolyData(reader.GetPolyData(), 1e-6);
		if (model.IsClosed())
			Log::info("Import Model: {}, is closed!", modelFileName);
		else
			Log::info("Import Model: {}, is not closed!", modelFileName);
	}

}
