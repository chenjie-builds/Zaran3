#include "DEMSolverParam.h"
#include "Log.h"

namespace zaran
{

void DEMSolverParam::Init()
{
    SolverParam::Init();

    if (GlobalData::IsExist("dem.time_step"))
        m_dt = GlobalData::GetDouble("dem.time_step");

    if (GlobalData::IsExist("dem.end_time"))
        m_end_time = GlobalData::GetDouble("dem.end_time");

    if (GlobalData::IsExist("dem.max_iter"))
        m_max_iter = GlobalData::GetInt("dem.max_iter");

    if (GlobalData::IsExist("dem.output_iter"))
        m_output_iter = GlobalData::GetInt("dem.output_iter");

    if (GlobalData::IsExist("dem.contact_model"))
        m_contact_model = GlobalData::GetString("dem.contact_model");

    if (GlobalData::IsExist("dem.particle_file"))
        m_particle_file = GlobalData::GetString("dem.particle_file");

    if (GlobalData::IsExist("dem.young_modulus"))
        m_young_modulus = GlobalData::GetDouble("dem.young_modulus");

    if (GlobalData::IsExist("dem.poisson_ratio"))
        m_poisson_ratio = GlobalData::GetDouble("dem.poisson_ratio");

    if (GlobalData::IsExist("dem.friction_coeff"))
        m_friction_coeff = GlobalData::GetDouble("dem.friction_coeff");

    if (GlobalData::IsExist("dem.restitution_coeff"))
        m_restitution_coeff = GlobalData::GetDouble("dem.restitution_coeff");

    // 重力向量（分量分别读取）
    if (GlobalData::IsExist("dem.gravity_x"))
        m_gravity.x() = GlobalData::GetDouble("dem.gravity_x");
    if (GlobalData::IsExist("dem.gravity_y"))
        m_gravity.y() = GlobalData::GetDouble("dem.gravity_y");
    if (GlobalData::IsExist("dem.gravity_z"))
        m_gravity.z() = GlobalData::GetDouble("dem.gravity_z");

    Log::info("DEM SolverParam Init: dt={:E}, contact_model={}, particle_file={}",
              m_dt, m_contact_model, m_particle_file);
}

} // namespace zaran
