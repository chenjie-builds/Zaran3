#include "DEMFieldSimulation.h"
#include "GlobalData.h"
#include "Log.h"
#include "File.h"
#include "Visual.h"
#include "DEMSolver.h"
#include "CavityGasModel.h"
#include <fstream>
#include <cstdio>
#include <filesystem>

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

    // 可视化格式：output.particle_format = tecplot | vtp | both | none
    //   tecplot（默认）：写 Tecplot ASCII（result/particles_<iter>.dat 与 bonds_<iter>.dat）
    //   vtp           ：写 VTK PolyData（result/particles_<iter>.vtp，供 ParaView）
    std::string format = GlobalData::IsExist("output.particle_format")
                             ? GlobalData::GetString("output.particle_format")
                             : "tecplot";
    if (format != "tecplot" && format != "vtp" && format != "both" && format != "none")
    {
        Log::warn("output.particle_format='{}' 无法识别（应为 tecplot|vtp|both|none），按 tecplot 处理",
                  format);
        format = "tecplot";
    }

    Visual vis;
    if (format == "tecplot" || format == "both")
    {
        // 键合文件是否写出：沿用 output.bond_details（默认 true）。
        // 该开关在 ParaView 模式下只影响数组多少，在 Tecplot 模式下决定是否
        // 生成 bonds_<iter>.dat（键合可达数万条，该文件通常是最大的一个）。
        const bool with_bonds = !GlobalData::IsExist("output.bond_details")
                                || GlobalData::GetBool("output.bond_details");
        const double solution_time = GlobalData::IsExist("dem.current_time")
                                         ? GlobalData::GetDouble("dem.current_time") : 0.0;
        vis.WriteParticleTecplotASCII(m_field_manager, iter, solution_time, with_bonds);
    }
    if (format == "vtp" || format == "both")
    {
        vis.WriteParticleVTP(m_field_manager, iter);
    }

    // 备份粒子 CSV
    std::string iter_folder = work_dir + "/" + back_folder + "/iter=" + std::to_string(iter);
    CreateFolder(iter_folder);
    auto* dem_solver = m_dem_field->GetDEMSolver().get();
    dem_solver->BackupField(iter_folder);

    // 气腔加载诊断：result/cavity_gas.csv（每帧追加一行）
    SaveCavityGasReport(iter, result_folder);
    // 破碎形态诊断：result/fragment_stats.csv（出现断键后每帧追加一行）
    SaveFragmentStats(iter, result_folder);
}

void DEMFieldSimulation::SaveFragmentStats(int iter, const std::string& result_folder) const
{
    auto* solver = m_dem_field->GetDEMSolver().get();
    if (solver == nullptr) return;
    const DEMSolver::EnergyBudget budget = solver->GetEnergyBudget();
    if (budget.broken_bonds == 0) return;   // 未起裂不写，避免无意义的空表

    const std::string work_dir = GlobalData::GetString("work_dir");
    const DEMSolver::FragmentStats st = solver->ComputeFragmentStats();
    const double time = GlobalData::IsExist("dem.current_time")
                            ? GlobalData::GetDouble("dem.current_time") : 0.0;

    std::string path = work_dir + "/" + result_folder + "/fragment_stats.csv";
    // 注意：本文件只在"已出现断键"时才写，所以 iter == 0 时往往还没有文件，
    // 用 iter==0 判断表头会永远不成立（踩过一次）。改为看文件是否存在/为空。
    const bool need_header = !std::filesystem::exists(path)
        || std::filesystem::file_size(path) == 0;
    std::ofstream fout(path, std::ios::app);
    if (!fout.is_open())
    {
        Log::warn("DEMFieldSimulation: 无法写入 '{}'", path);
        return;
    }
    if (need_header)
    {
        fout << "iter,time,fragments,isolated_particles,largest_fragment,"
                "largest_fragment_mass_fraction,fragments_ge_20,"
                "fragments_ge_20_mass_fraction,largest_crack_band,"
                "active_bonds,broken_bonds,broken_fraction,max_radius_m";
        for (int i = 0; i < 12; ++i)
        {
            fout << ",radial_broken_" << i;
        }
        fout << "\n";
    }
    char line[1024];
    const long long total_bonds = static_cast<long long>(
        budget.active_bonds + budget.broken_bonds);
    std::snprintf(line, sizeof(line),
        "%d,%.12e,%lld,%lld,%lld,%.12e,%lld,%.12e,%lld,%lld,%lld,%.12e,%.12e",
        iter, time,
        static_cast<long long>(st.fragments),
        static_cast<long long>(st.isolated_particles),
        static_cast<long long>(st.largest_fragment),
        st.largest_fragment_mass_fraction,
        static_cast<long long>(st.fragments_ge_20),
        st.fragments_ge_20_mass_fraction,
        static_cast<long long>(st.largest_crack_band),
        static_cast<long long>(budget.active_bonds),
        static_cast<long long>(budget.broken_bonds),
        total_bonds > 0 ? static_cast<double>(budget.broken_bonds) / total_bonds : 0.0,
        st.max_radius);
    fout << line;
    for (int i = 0; i < 12; ++i)
    {
        std::snprintf(line, sizeof(line), ",%.9f", st.radial_broken_fraction[i]);
        fout << line;
    }
    fout << "\n";
}

void DEMFieldSimulation::SaveCavityGasReport(int iter, const std::string& result_folder) const
{
    auto* solver = m_dem_field->GetDEMSolver().get();
    if (solver == nullptr || !solver->GetCavityGas().IsEnabled()) return;

    const std::string work_dir = GlobalData::GetString("work_dir");
    const CavityGasModel& gas = solver->GetCavityGas();
    const DEMSolver::EnergyBudget budget = solver->GetEnergyBudget();
    const double time = GlobalData::IsExist("dem.current_time")
                            ? GlobalData::GetDouble("dem.current_time") : 0.0;

    // 能量账本：U0 + W_ext = U + KE + ΣE_bond + ΣE_fracture + ΣE_dissipation
    //                        + ΣE_contact_elastic + ΣE_rebound_protection
    const double source = gas.GetInitialInternalEnergy() + gas.GetExternalWork();
    const double stored = gas.GetInternalEnergy()
        + budget.kinetic_translation + budget.kinetic_rotation
        + budget.bond_elastic + budget.bond_fracture + budget.dissipation
        + budget.contact_elastic + budget.rebound_loss;
    const double residual = source > 0.0 ? (stored - source) / source : 0.0;

    std::string path = work_dir + "/" + result_folder + "/cavity_gas.csv";
    std::ofstream fout(path, std::ios::app);
    if (!fout.is_open())
    {
        Log::warn("DEMFieldSimulation: 无法写入 '{}'", path);
        return;
    }
    if (iter == 0)
    {
        fout << "iter,time,cavity_pressure_Pa,cavity_volume_m3,cavity_area_m2,"
                "cavity_area_geom_m2,cavity_radius_eff_m,gas_internal_energy_J,"
                "gas_work_J,external_work_J,kinetic_translation_J,kinetic_rotation_J,"
                "bond_elastic_J,bond_fracture_J,dissipation_J,contact_elastic_J,"
                "rebound_loss_J,energy_residual_rel,"
                "max_speed_m_s,max_fragment_distance_m,active_bonds,broken_bonds,"
                "ring_particles,vented,gas_force_sum_N,ring_max_gap_rad\n";
    }
    char line[1400];
    std::snprintf(line, sizeof(line),
        "%d,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,"
        "%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%lld,%lld,%lld,%d,%.12e,%.12e\n",
        iter, time,
        gas.GetPressure(), gas.GetVolume(), gas.GetArea(), gas.GetAreaGeometric(),
        gas.GetEffectiveRadius(), gas.GetInternalEnergy(), gas.GetGasWork(),
        gas.GetExternalWork(),
        budget.kinetic_translation, budget.kinetic_rotation,
        budget.bond_elastic, budget.bond_fracture, budget.dissipation,
        budget.contact_elastic, budget.rebound_loss, residual,
        budget.max_speed, budget.max_fragment_distance,
        static_cast<long long>(budget.active_bonds),
        static_cast<long long>(budget.broken_bonds),
        static_cast<long long>(gas.GetRingCount()), gas.IsVented() ? 1 : 0,
        gas.GetAppliedForceSum(), gas.GetRingMaxGap());
    fout << line;
}

void DEMFieldSimulation::SolveField()
{
    // 诊断 CSV 是"每次运行一份"的追加型文件：上一次运行的残留会被当成本次数据
    // （表现为表头出现两次、或校验脚本读到两轮拼起来的行）。开跑前先删除。
    {
        const std::string work_dir = GlobalData::GetString("work_dir");
        const std::string folder = GlobalData::IsExist("output.result_folder")
            ? GlobalData::GetString("output.result_folder") : "result";
        for (const char* name : { "cavity_gas.csv", "fragment_stats.csv" })
        {
            const std::string path = work_dir + "/" + folder + "/" + name;
            std::error_code ec;
            std::filesystem::remove(path, ec);
        }
    }

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
