#include "DEMFieldSimulation.h"
#include "GlobalData.h"
#include "Log.h"
#include "File.h"
#include "Visual.h"

namespace zaran
{

DEMFieldSimulation::DEMFieldSimulation(shared_ptr<FieldManager> field_manager)
    : m_field_manager(std::move(field_manager))
{
    // 找到第一个 DEMField
    for (size_t i = 0; i < m_field_manager->GetFieldNum(); ++i)
    {
        auto f = std::dynamic_pointer_cast<DEMField>(m_field_manager->GetField(i));
        if (f)
        {
            m_dem_field = f;
            break;
        }
    }
}

void DEMFieldSimulation::Initialize()
{
    if (!m_dem_field)
    {
        Log::error("DEMFieldSimulation: no DEMField found in FieldManager!");
        return;
    }
    m_dem_field->GetSolver()->Init();
    Log::info("DEMFieldSimulation: initialized with {} particles",
              m_dem_field->GetDEMData()->GetParticleNum());
}

bool DEMFieldSimulation::ContinueSolve() const
{
    auto* para = m_dem_field->GetDEMSolverParam().get();
    int current_iter = GlobalData::IsExist("dem.current_iter")
                           ? GlobalData::GetInt("dem.current_iter") : 0;
    double current_time = GlobalData::IsExist("dem.current_time")
                              ? GlobalData::GetDouble("dem.current_time") : 0.0;

    if (current_iter >= para->GetMaxIter()) return false;
    if (current_time >= para->GetEndTime()) return false;
    return true;
}

void DEMFieldSimulation::SolveOneStep()
{
    for (size_t i = 0; i < m_field_manager->GetFieldNum(); ++i)
    {
        auto f = std::dynamic_pointer_cast<DEMField>(m_field_manager->GetField(i));
        if (f) f->GetSolver()->Solve();
    }
}

void DEMFieldSimulation::SaveFieldData(int iter) const
{
    std::string work_dir    = GlobalData::GetString("work_dir");
    std::string result_folder = GlobalData::IsExist("output.result_folder")
                                    ? GlobalData::GetString("output.result_folder")
                                    : "result";
    std::string back_folder   = GlobalData::IsExist("init.backup_folder")
                                    ? GlobalData::GetString("init.backup_folder")
                                    : "backup";

    // 写粒子 VTP
    Visual vis;
    vis.WriteParticleVTP(m_field_manager, iter);

    // 备份粒子 CSV
    std::string iter_folder = work_dir + "/" + back_folder + "/iter=" + std::to_string(iter);
    CreateFolder(iter_folder);
    auto* dem_solver = m_dem_field->GetDEMSolver().get();
    dem_solver->BackupField(iter_folder);
}

void DEMFieldSimulation::SolveField()
{
    Initialize();

    GlobalData::Update("dem.current_iter", 0);
    GlobalData::Update("dem.current_time", 0.0);

    auto* para = m_dem_field->GetDEMSolverParam().get();
    double dt = para->GetTimeStep();
    int output_iter = para->GetOutputIter();

    SaveFieldData(0);

    while (ContinueSolve())
    {
        int current_iter = GlobalData::GetInt("dem.current_iter") + 1;
        GlobalData::Update("dem.current_iter", current_iter);
        double current_time = GlobalData::GetDouble("dem.current_time") + dt;
        GlobalData::Update("dem.current_time", current_time);

        SolveOneStep();

        if (current_iter % output_iter == 0)
        {
            Log::info("DEM iter={}, time={:E}, particles={}",
                      current_iter, current_time,
                      m_dem_field->GetDEMData()->GetParticleNum());
            SaveFieldData(current_iter);
        }
    }

    // 最终输出（若最后一步恰为输出步则跳过，避免重复）
    int final_iter = GlobalData::GetInt("dem.current_iter");
    if (final_iter % output_iter != 0)
    {
        SaveFieldData(final_iter);
    }
    Log::info("DEM simulation finished. Total iter={}", final_iter);
}

} // namespace zaran
