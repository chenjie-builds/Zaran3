#include "FieldSimulation.h"
#include "File.h"
#include "FlowSolver.h"
#include "Log.h"
#include "MathBasic.h"
#include "NSFieldFN.h"
#include <fstream>
using namespace zaran;

NSFieldSimulation::NSFieldSimulation(const shared_ptr<FieldManager> &field_manager)
{
    m_field_manager = field_manager;
    m_res_flag = false;
}

NSFieldSimulation::~NSFieldSimulation() = default;

void NSFieldSimulation::Initialize() const
{
    if (!GlobalData::IsExist("iteration.current_iter"))
    {
        GlobalData::Update("iteration.current_iter", 0);
    }
    if (!GlobalData::IsExist("iteration.current_time"))
    {
        double startTime = GlobalData::GetDouble("iteration.start_time");
        GlobalData::Update("iteration.current_time", startTime);
    }
    for (size_t iter_field = 0; iter_field < m_field_manager->GetFieldNum(); iter_field++)
    {
        m_field_manager->GetField(iter_field)->GetSolver()->Init();
    }
}

void NSFieldSimulation::SaveDataTecplot() const
{
     //m_visual->WriteTecASCII(m_field_manager);
    m_visual->WriteTecplotBinary(m_field_manager);
    // m_visual->WriteVtkASCII(m_field_manager);
    // m_visual->WriteCGNS(m_field_manager);
}

void NSFieldSimulation::SolveField()
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
        int currentIter = GlobalData::GetInt("iteration.current_iter");
        GlobalData::Update("iteration.current_iter", ++currentIter);
        PreSolve();
        SolveOneStep();
        PostSolve();
    }
    PostSolve();
    SaveFieldData();
}

void NSFieldSimulation::CalcResidual()
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

void NSFieldSimulation::SaveFieldData()
{
    SaveDataTecplot();
    std::string back_dir = GlobalData::GetString("init.backup_folder");
    int currentIter = GlobalData::GetInt("iteration.current_iter");
    back_dir += "/iter=" + std::to_string(currentIter);
    std::string work_dir = GlobalData::GetString("work_dir");
    back_dir = work_dir + "/" + back_dir;
    CreateFolder(back_dir);
    BackupFieldData(back_dir);
    BackupResidual(back_dir);
    BackupLog(back_dir);
    BackupGlobalData(back_dir);
}

void NSFieldSimulation::BackupFieldData(std::string &back_folder)
{
    for (size_t iter_field = 0; iter_field < m_field_manager->GetFieldNum(); iter_field++)
    {
        auto solver = m_field_manager->GetField(iter_field)->GetSolver();
        solver->BackupField(back_folder);
    }
}

void NSFieldSimulation::BackupResidual(std::string &back_folder)
{
    std::string residual_file = GlobalData::GetString("output.residual_file");
    std::string residual_file_back = back_folder + "/" + residual_file;
    if (IsFileExist(residual_file_back) == false)
    {
        return;
    }

    if (IsFileExist(residual_file_back) == true)
    {
        zaran::RemoveFile(residual_file_back);
    }
    CopySingleFile(residual_file, residual_file_back);
}

void NSFieldSimulation::BackupLog(std::string &back_folder)
{
    std::string log_file = "log.txt";
    std::string log_file_back = back_folder + "/" + log_file;
    if (IsFileExist(log_file_back) == false)
    {
        return;
    }

    if (IsFileExist(log_file_back) == true)
    {
        RemoveFile(log_file_back);
    }
    CopySingleFile(log_file, log_file_back);
}

void NSFieldSimulation::BackupGlobalData(std::string &back_folder)
{
    std::string global_file = "zaran.ini";
    std::string global_file_back = back_folder + "/" + global_file;
    if (IsFileExist(global_file_back) == false)
    {
        return;
    }
    if (IsFileExist(global_file_back) == true)
    {
        RemoveFile(global_file_back);
    }
    GlobalData::Backup(back_folder);
}

bool NSFieldSimulation::ContinueSolve()
{
    double end_time = GlobalData::GetDouble("iteration.end_time");
    int current_iter = GlobalData::GetInt("iteration.current_iter");
    int cal_res_iter = GlobalData::GetInt("iteration.residual_interval");
    int max_iter = GlobalData::GetInt("iteration.max_iter");
    int min_res = GlobalData::GetInt("iteration.min_residual_order");
    double currentTime = GlobalData::GetDouble("iteration.current_time");
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
void NSFieldSimulation::SaveResidual() const
{
    int current_iter = GlobalData::GetInt("iteration.current_iter");
    string residual_file = GlobalData::GetString("output.residual_file");
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
        fout << current_iter << "\t\t" << GlobalData::GetDouble("iteration.current_time") << "\t\t" << m_res_max << "\t\t"
             << m_res_ave << std::endl;
        fout.close();
    }
}

void NSFieldSimulation::SolveOneStep()
{
    for (size_t iter_field = 0; iter_field < m_field_manager->GetFieldNum(); iter_field++)
    {
        auto solver = m_field_manager->GetField(iter_field)->GetSolver();
        solver->Solve();
    }
}

void NSFieldSimulation::PreSolve()
{
    int currentIter = GlobalData::GetInt("iteration.current_iter");
    GlobalData::Update("iteration.current_iter", ++currentIter);
    double min_dt_global = LARGE_NUMBER;
    for (size_t iter_field = 0; iter_field < m_field_manager->GetFieldNum(); iter_field++)
    {
        auto solver = m_field_manager->GetField(iter_field)->GetSolver();
        solver->Preprocess();
    }
    CalcTimeStep();
}

void NSFieldSimulation::PostSolve()
{
    for (size_t iter_field = 0; iter_field < m_field_manager->GetFieldNum(); iter_field++)
    {
        auto solver = m_field_manager->GetField(iter_field)->GetSolver();
        solver->Postprocess();
    }
    CommFieldData();
    int currentIter = GlobalData::GetInt("iteration.current_iter");
    int calResidualIter = GlobalData::GetInt("iteration.residual_interval");
    int writeFieldIter = GlobalData::GetInt("iteration.write_interval");
    if (currentIter % calResidualIter == 0)
    {
        CalcResidual();
        Log::info("iter= {}, dt={:E}, res_max= {:E}, res_ave= {:E}", GlobalData::GetInt("iteration.current_iter"),
                  GlobalData::GetDouble("iteration.dt"), m_res_max, m_res_ave);
        SaveResidual();
    }
    if (currentIter % writeFieldIter == 0)
    {
        Log::info("Backup field data...");
        Log::info("physical time: {:E}", GlobalData::GetDouble("iteration.current_time"));
        SaveFieldData();
    }
}

void NSFieldSimulation::CommFieldData()
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
        auto &recv_data_name = comm_info->GetRecvDataName();
        int recv_node_num = comm_info->GetRecvNodeNum();
#ifdef USE_OMP
#pragma omp parallel for
#endif
        for (int i_recv_node = 0; i_recv_node < recv_node_num; i_recv_node++)
        {
            auto send_field = m_field_manager->GetField(comm_info->GetTgtFieldIdx()[i_recv_node]);
            auto send_field_data = send_field->GetData();
            int idx_send_field = comm_info->GetTgtFieldIdx()[i_recv_node];
            int idx_send_node = comm_info->GetTgtNodeIdx()[i_recv_node];
            for (size_t i_recv_name = 0; i_recv_name < recv_data_name.size(); i_recv_name++)
            {
                std::string data_name = recv_data_name[i_recv_name];
                double data = send_field_data->GetData(data_name, idx_send_node);
                comm_info->UpdateSrcDataBuffer(i_recv_name, i_recv_node, data);
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
        auto &recv_data_name = comm_info->GetRecvDataName();
        int recv_node_num = comm_info->GetRecvNodeNum();
        auto recv_node = comm_info->GetSrcDataIdx();
#ifdef USE_OMP
#pragma omp parallel for
#endif
        for (int i_recv_node = 0; i_recv_node < recv_node_num; i_recv_node++)
        {
            for (size_t i_recv_name = 0; i_recv_name < recv_data_name.size(); i_recv_name++)
            {
                std::string data_name = recv_data_name[i_recv_name];
                field_data->GetData(data_name, recv_node[i_recv_node]) =
                    comm_info->GetSrcDataBuffer(i_recv_name, i_recv_node);
            }
        }
    }
}

void NSFieldSimulation::CalcTimeStep()
{
    double min_dt_global = LARGE_NUMBER;
    for (size_t iter_field = 0; iter_field < m_field_manager->GetFieldNum(); iter_field++)
    {
        auto solver = m_field_manager->GetField(iter_field)->GetSolver();
        auto ns_solver = std::dynamic_pointer_cast<NSSolver>(solver);
        ns_solver->CalcTimeStepLocal();
        double min_dt_local = LARGE_NUMBER;
        ns_solver->CalcMinTimeStep(min_dt_local);
        min_dt_global = Min(min_dt_global, min_dt_local);
    }
    GlobalData::Update("iteration.dt", min_dt_global);
    double current_time = GlobalData::GetDouble("iteration.current_time");
    double end_time = GlobalData::GetDouble("iteration.end_time");
    if (current_time + min_dt_global > end_time)
    {
        min_dt_global = end_time - current_time;
        current_time = end_time;
    }
    else
    {
        current_time += min_dt_global;
    }
    GlobalData::Update("iteration.current_time", current_time);
    int isSteady = GlobalData::GetBool("task.is_steady") ? 1 : 0;
    if (isSteady == 0)
    {
        for (size_t iter_field = 0; iter_field < m_field_manager->GetFieldNum(); iter_field++)
        {
            auto solver = m_field_manager->GetField(iter_field)->GetSolver();
            auto ns_solver = std::dynamic_pointer_cast<NSSolver>(solver);
            ns_solver->ReduceTimeStep(min_dt_global);
        }
    }
}
