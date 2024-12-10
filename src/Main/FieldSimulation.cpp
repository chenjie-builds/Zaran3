#include "FieldSimulation.h"
#include "Log.h"
#include <fstream>
#include "FlowSolver.h"
#include "MathBasic.h"
#include "File.h"
#include "NSFieldFN.h"
using namespace zaran;

FieldSimulation::FieldSimulation(shared_ptr<FieldManager> global_Field)
{
	m_field_manager = global_Field;
	m_res_flag = false;
}

FieldSimulation::~FieldSimulation()
{

}

void FieldSimulation::Initialize()
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
	for (size_t iter_field = 0; iter_field < m_field_manager->GetFieldNum(); iter_field++)
	{
		m_field_manager->GetField(iter_field)->GetSolver()->Init();
	}
}

void FieldSimulation::SaveDataTecplot()
{
	//m_visual->WriteTecASCII(m_field_manager);
	m_visual->WriteTecplotBinary(m_field_manager);


}
void FieldSimulation::SaveDataVTK(std::ostream& os)
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

void FieldSimulation::SolveField()
{
	Log::info("Start to solve field!");
	Initialize();
	Log::info("Initialize finished!");
	Log::info("Save init field data...");
	SaveFieldData();
	Log::info("Save init field data finished!");
	SaveResidual();
	while (ContinueSolve())
	{
		PreSolve();
		SolveOneStep();
		PostSolve();
	}
	PostSolve();
	SaveFieldData();
}

void FieldSimulation::CalcResidual()
{
	m_res_max = -LARGE_NUMBER;
	m_res_ave = 0;
	for (size_t iter_field = 0; iter_field < m_field_manager->GetFieldNum(); iter_field++)
	{
		auto field = std::static_pointer_cast<FieldNS>(m_field_manager->GetField(iter_field));
		field->CalcResidual();
		auto res_info = field->GetResInfo();
		m_res_max = Max(m_res_max, res_info->GetInfNorm(0));
		m_res_ave = m_res_ave + res_info->GetL2Norm(0);
	}
	m_res_flag = true;
}

void FieldSimulation::SaveFieldData()
{
	SaveDataTecplot();
	std::string back_dir = GlobalData::GetString("backupFieldFolder");
	int currentIter = GlobalData::GetInt("currentIter");
	back_dir += "/iter=" + std::to_string(currentIter);
	std::string work_dir = GlobalData::GetString("work_dir");
	back_dir = work_dir + "/" + back_dir;
	CreateFolder(back_dir);
	BackupFieldData(back_dir);
	BackupResidual(back_dir);
	BackupLog(back_dir);
	BackupGlobalData(back_dir);
}

void FieldSimulation::BackupFieldData(std::string& back_folder)
{
	for (size_t iter_field = 0; iter_field < m_field_manager->GetFieldNum(); iter_field++)
	{
		auto solver = m_field_manager->GetField(iter_field)->GetSolver();
		solver->BackupField(back_folder);
	}
}

void FieldSimulation::BackupResidual(std::string& back_folder)
{
	std::string residual_file = GlobalData::GetString("residualFileName");
	std::string residual_file_back = back_folder + "/" + residual_file;
	if (IsFileExist(residual_file_back) == false)
	{
		Log::warn("residual file:{} is not exist!", residual_file_back);
		return;
	}
	if (IsFileExist(residual_file_back) == true)
	{
		DeleteFile(residual_file_back);
	}
	CopyFile(residual_file, residual_file_back);
}

void FieldSimulation::BackupLog(std::string& back_folder)
{
	std::string log_file = "log.txt";
	std::string log_file_back = back_folder + "/" + log_file;
	if (IsFileExist(log_file_back) == false)
	{
		Log::warn("logFile:{} is not exist!", log_file_back);
		return;
	}
	if (IsFileExist(log_file_back) == true)
	{
		DeleteFile(log_file_back);
	}
	CopyFile(log_file, log_file_back);
}

void FieldSimulation::BackupGlobalData(std::string& back_folder)
{
	std::string global_file = "zaran.ini";
	std::string global_file_back = back_folder + "/" + global_file;
	if (IsFileExist(global_file_back) == false)
	{
		Log::warn("gloabal data file:{} is not exist!", global_file_back);
		return;
	}
	if (IsFileExist(global_file_back) == true)
	{
		DeleteFile(global_file_back);
	}
	GlobalData::Backup(back_folder);
}

bool FieldSimulation::ContinueSolve()
{
	double end_time = GlobalData::GetDouble("endTime");
	int current_iter = GlobalData::GetInt("currentIter");
	int cal_res_iter = GlobalData::GetInt("calResidualIter");
	int max_iter = GlobalData::GetInt("maxIter");
	int min_res = GlobalData::GetInt("minResidual");
	double currentTime = GlobalData::GetDouble("currentTime");
	if (m_res_flag == true)
	{
		if (log10(m_res_max) < -min_res && current_iter > cal_res_iter)
		{
			Log::info("Max Residual is small than {}, stop compute!", min_res);
			return false;
		}
	}
	if (current_iter > max_iter)
	{
		Log::info("Max Iter={}, stop compute!", max_iter);
		return false;
	}
	if (currentTime > end_time || abs(currentTime - end_time) < SMALL_NUMBER)
	{
		Log::info("Max time={}, stop compute!", end_time);
		return false;
	}
	return true;
}
void FieldSimulation::SaveResidual()
{
	int current_iter = GlobalData::GetInt("currentIter");
	string residual_file = GlobalData::GetString("residualFileName");
	string work_dir = GlobalData::GetString("work_dir");
	residual_file = work_dir + "/" + residual_file;

	if (current_iter == 0)
	{
		std::ofstream fout(residual_file);
		fout << "variables=step, time, Res_Linf, Res_L2\n";
		fout.close();
	}
	else
	{
		std::ofstream fout(residual_file, std::ios::app);
		fout << current_iter << "\t\t" << GlobalData::GetDouble("currentTime") << "\t\t" << m_res_max << "\t\t" << m_res_ave << std::endl;
		fout.close();
	}
}

void FieldSimulation::SolveOneStep()
{
	for (size_t iter_field = 0; iter_field < m_field_manager->GetFieldNum(); iter_field++)
	{
		auto solver = m_field_manager->GetField(iter_field)->GetSolver();
		solver->Solve();
	}
}

void FieldSimulation::PreSolve()
{
	int currentIter = GlobalData::GetInt("currentIter");
	GlobalData::Update("currentIter", ++currentIter);
	for (size_t iter_field = 0; iter_field < m_field_manager->GetFieldNum(); iter_field++)
	{
		auto solver = m_field_manager->GetField(iter_field)->GetSolver();
		solver->Preprocess();
	}
}

void FieldSimulation::PostSolve()
{
	for (size_t iter_field = 0; iter_field < m_field_manager->GetFieldNum(); iter_field++)
	{
		auto solver = m_field_manager->GetField(iter_field)->GetSolver();
		solver->Postprocess();
	}
	CommFieldData();
	int currentIter = GlobalData::GetInt("currentIter");
	int calResidualIter = GlobalData::GetInt("calResidualIter");
	int writeFieldIter = GlobalData::GetInt("writeFieldIter");
	if (currentIter % calResidualIter == 0 && ContinueSolve())
	{
		CalcResidual();
		Log::info("iter= {}, dt={:E}, res_max= {:E}, res_ave= {:E}", GlobalData::GetInt("currentIter"), GlobalData::GetDouble("dt"), m_res_max, m_res_ave);
		SaveResidual();
	}
	if (currentIter % writeFieldIter == 0 && ContinueSolve())
	{
		Log::info("Backup field data...");
		Log::info("physical time: {:E}", GlobalData::GetDouble("currentTime"));
		SaveFieldData();
	}
}

void FieldSimulation::CommFieldData()
{
	for (size_t iter_field = 0; iter_field < m_field_manager->GetFieldNum(); iter_field++)
	{
		auto field = m_field_manager->GetField(iter_field);
		auto field_data = field->GetData();
		auto comm_info = m_field_manager->GetFieldDataCommInfo(iter_field);
		if (comm_info == nullptr)
		{
			continue;
		}
		auto& recv_data_name = comm_info->GetRecvDataName();
		int recv_node_num = comm_info->GetRecvNodeNum();
		for (int i_recv_node = 0; i_recv_node < recv_node_num; i_recv_node++)
		{
			auto send_field = m_field_manager->GetField(comm_info->GetIdxSendField()[i_recv_node]);
			auto send_field_data = send_field->GetData();
			int idx_send_field = comm_info->GetIdxSendField()[i_recv_node];
			int idx_send_node = comm_info->GetIdxSendNode()[i_recv_node];
			for (size_t i_recv_name = 0; i_recv_name < recv_data_name.size(); i_recv_name++)
			{
				std::string data_name = recv_data_name[i_recv_name];
				double data = send_field_data->GetData(data_name, idx_send_node);
				comm_info->SetRecvDataCache(i_recv_name, i_recv_node, data);
			}
		}
	}
	for (size_t iter_field = 0; iter_field < m_field_manager->GetFieldNum(); iter_field++)
	{
		auto field = m_field_manager->GetField(iter_field);
		auto field_data = field->GetData();
		auto comm_info = m_field_manager->GetFieldDataCommInfo(iter_field);
		if (comm_info == nullptr)
		{
			continue;
		}
		auto& recv_data_name = comm_info->GetRecvDataName();
		int recv_node_num = comm_info->GetRecvNodeNum();
		auto recv_node = comm_info->GetIdxRecvNode();
		for (int i_recv_node = 0; i_recv_node < recv_node_num; i_recv_node++)
		{
			for (size_t i_recv_name = 0; i_recv_name < recv_data_name.size(); i_recv_name++)
			{
				std::string data_name = recv_data_name[i_recv_name];
				field_data->GetData(data_name, recv_node[i_recv_node]) = comm_info->GetRecvDataCache(i_recv_name, i_recv_node);
			}
		}
	}
}
