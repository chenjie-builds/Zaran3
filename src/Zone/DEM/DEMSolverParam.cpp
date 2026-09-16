#include "DEMSolverParam.h"
#include "Log.h"
#include "ZaranError.h"

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

    if (GlobalData::IsExist("dem.density"))
        m_density = GlobalData::GetDouble("dem.density");
    if (GlobalData::IsExist("dem.reaction_enabled"))
        m_reaction_enabled = GlobalData::GetBool("dem.reaction_enabled");
    if (GlobalData::IsExist("dem.ignition_temperature"))
        m_ignition_temperature = GlobalData::GetDouble("dem.ignition_temperature");
    if (GlobalData::IsExist("dem.mixed_phase_threshold"))
        m_mixed_phase_threshold = GlobalData::GetDouble("dem.mixed_phase_threshold");
    if (GlobalData::IsExist("dem.burnout_threshold"))
        m_burnout_threshold = GlobalData::GetDouble("dem.burnout_threshold");
    if (GlobalData::IsExist("dem.max_temperature"))
        m_max_temperature = GlobalData::GetDouble("dem.max_temperature");
    if (GlobalData::IsExist("dem.jwl_a")) m_jwl_a = GlobalData::GetDouble("dem.jwl_a");
    if (GlobalData::IsExist("dem.jwl_b")) m_jwl_b = GlobalData::GetDouble("dem.jwl_b");
    if (GlobalData::IsExist("dem.jwl_r1")) m_jwl_r1 = GlobalData::GetDouble("dem.jwl_r1");
    if (GlobalData::IsExist("dem.jwl_r2")) m_jwl_r2 = GlobalData::GetDouble("dem.jwl_r2");
    if (GlobalData::IsExist("dem.jwl_omega")) m_jwl_omega = GlobalData::GetDouble("dem.jwl_omega");
    if (GlobalData::IsExist("dem.burn_a")) m_burn_a = GlobalData::GetDouble("dem.burn_a");
    if (GlobalData::IsExist("dem.burn_b")) m_burn_b = GlobalData::GetDouble("dem.burn_b");
    if (GlobalData::IsExist("dem.burn_pressure_threshold"))
        m_burn_pressure_threshold = GlobalData::GetDouble("dem.burn_pressure_threshold");
    if (GlobalData::IsExist("dem.burn_temperature_threshold"))
        m_burn_temperature_threshold = GlobalData::GetDouble("dem.burn_temperature_threshold");
    if (GlobalData::IsExist("dem.gas_temperature_reference"))
        m_gas_temperature_reference = GlobalData::GetDouble("dem.gas_temperature_reference");
    if (GlobalData::IsExist("dem.pressure_cap")) m_pressure_cap = GlobalData::GetDouble("dem.pressure_cap");
    if (GlobalData::IsExist("dem.bond_break_strain"))
        m_bond_break_strain = GlobalData::GetDouble("dem.bond_break_strain");
    if (GlobalData::IsExist("dem.bond_peak_strain"))
        m_bond_peak_strain = GlobalData::GetDouble("dem.bond_peak_strain");
    if (GlobalData::IsExist("dem.reaction_weakening"))
        m_reaction_weakening = GlobalData::GetDouble("dem.reaction_weakening");
    if (GlobalData::IsExist("dem.ignition_center_x"))
        m_ignition_center.x() = GlobalData::GetDouble("dem.ignition_center_x");
    if (GlobalData::IsExist("dem.ignition_center_y"))
        m_ignition_center.y() = GlobalData::GetDouble("dem.ignition_center_y");
    if (GlobalData::IsExist("dem.ignition_center_z"))
        m_ignition_center.z() = GlobalData::GetDouble("dem.ignition_center_z");
    if (GlobalData::IsExist("dem.ignition_radius"))
        m_ignition_radius = GlobalData::GetDouble("dem.ignition_radius");
    if (GlobalData::IsExist("dem.ignition_center_temperature"))
        m_ignition_center_temperature = GlobalData::GetDouble("dem.ignition_center_temperature");
    if (GlobalData::IsExist("dem.lattice_thickness"))
        m_lattice_thickness = GlobalData::GetDouble("dem.lattice_thickness");
    if (GlobalData::IsExist("dem.internal_heat_limit"))
        m_internal_heat_limit = GlobalData::GetDouble("dem.internal_heat_limit");
    if (GlobalData::IsExist("dem.voronoi_enabled"))
        m_voronoi_enabled = GlobalData::GetBool("dem.voronoi_enabled");
    if (GlobalData::IsExist("dem.voronoi_update_interval"))
        m_voronoi_update_interval = GlobalData::GetInt("dem.voronoi_update_interval");
    if (GlobalData::IsExist("dem.voronoi_max_volume_ratio"))
        m_voronoi_max_volume_ratio = GlobalData::GetDouble("dem.voronoi_max_volume_ratio");
    if (GlobalData::IsExist("dem.box.x_min")) m_box_x_min = GlobalData::GetDouble("dem.box.x_min");
    if (GlobalData::IsExist("dem.box.x_max")) m_box_x_max = GlobalData::GetDouble("dem.box.x_max");
    if (GlobalData::IsExist("dem.box.y_min")) m_box_y_min = GlobalData::GetDouble("dem.box.y_min");
    if (GlobalData::IsExist("dem.box.y_max")) m_box_y_max = GlobalData::GetDouble("dem.box.y_max");
    if (GlobalData::IsExist("dem.surface_energy"))
        m_surface_energy = GlobalData::GetDouble("dem.surface_energy");
    if (GlobalData::IsExist("dem.grain_boundary_energy"))
        m_grain_boundary_energy = GlobalData::GetDouble("dem.grain_boundary_energy");
    if (GlobalData::IsExist("dem.high_pressure_rnn"))
        m_high_pressure_rnn = GlobalData::GetDouble("dem.high_pressure_rnn");
    if (GlobalData::IsExist("dem.high_pressure_knn"))
        m_high_pressure_knn = GlobalData::GetDouble("dem.high_pressure_knn");
    if (GlobalData::IsExist("dem.high_pressure_exponent"))
        m_high_pressure_exponent = GlobalData::GetDouble("dem.high_pressure_exponent");
    if (GlobalData::IsExist("dem.normal_viscosity"))
        m_normal_viscosity = GlobalData::GetDouble("dem.normal_viscosity");
    if (GlobalData::IsExist("dem.tangential_viscosity"))
        m_tangential_viscosity = GlobalData::GetDouble("dem.tangential_viscosity");

    // --- 初始弹簧连接网络（可选）---
    if (GlobalData::IsExist("dem.spring_network_enabled"))
        m_spring_network_enabled = GlobalData::GetBool("dem.spring_network_enabled");
    if (GlobalData::IsExist("dem.spring_network_gap"))
        m_spring_network_gap = GlobalData::GetDouble("dem.spring_network_gap");
    if (GlobalData::IsExist("dem.spring_network_stiffness"))
        m_spring_network_stiffness = GlobalData::GetDouble("dem.spring_network_stiffness");
    if (GlobalData::IsExist("dem.spring_network_tangential_stiffness"))
        m_spring_network_tangential_stiffness = GlobalData::GetDouble("dem.spring_network_tangential_stiffness");
    if (GlobalData::IsExist("dem.spring_network_fracture_strain"))
        m_spring_network_fracture_strain = GlobalData::GetDouble("dem.spring_network_fracture_strain");

    // --- 刚性边界（可选）---
    if (GlobalData::IsExist("dem.rigid_boundary"))
        m_rigid_boundary_enabled = GlobalData::GetBool("dem.rigid_boundary");

    // 重力向量（分量分别读取）
    if (GlobalData::IsExist("dem.gravity_x"))
        m_gravity.x() = GlobalData::GetDouble("dem.gravity_x");
    if (GlobalData::IsExist("dem.gravity_y"))
        m_gravity.y() = GlobalData::GetDouble("dem.gravity_y");
    if (GlobalData::IsExist("dem.gravity_z"))
        m_gravity.z() = GlobalData::GetDouble("dem.gravity_z");

    if (m_dt <= 0.0)
        throw ZaranError("dem.time_step must be positive");
    if (m_end_time < 0.0)
        throw ZaranError("dem.end_time must be non-negative");
    if (m_max_iter < 0)
        throw ZaranError("dem.max_iter must be non-negative");
    if (m_output_iter <= 0)
        throw ZaranError("dem.output_iter must be positive");
    if (m_young_modulus <= 0.0)
        throw ZaranError("dem.young_modulus must be positive");
    if (m_poisson_ratio <= -1.0 || m_poisson_ratio >= 0.5)
        throw ZaranError("dem.poisson_ratio must be in (-1, 0.5)");
    if (m_friction_coeff < 0.0)
        throw ZaranError("dem.friction_coeff must be non-negative");
    if (m_restitution_coeff <= 0.0 || m_restitution_coeff > 1.0)
        throw ZaranError("dem.restitution_coeff must be in (0, 1]");
    if (m_density <= 0.0)
        throw ZaranError("dem.density must be positive");
    if (m_ignition_temperature <= 0.0 || m_max_temperature <= 0.0)
        throw ZaranError("DEM reaction temperatures must be positive");
    if (m_mixed_phase_threshold < 0.0 || m_mixed_phase_threshold >= 1.0
        || m_burnout_threshold <= m_mixed_phase_threshold || m_burnout_threshold > 1.0)
        throw ZaranError("DEM reaction phase thresholds are invalid");
    if (m_jwl_a < 0.0 || m_jwl_b < 0.0 || m_jwl_r1 <= 0.0 || m_jwl_r2 <= 0.0
        || m_jwl_omega < 0.0 || m_burn_a < 0.0 || m_burn_b < 0.0
        || m_burn_pressure_threshold < 0.0 || m_burn_temperature_threshold < 0.0
        || m_gas_temperature_reference <= 0.0 || m_pressure_cap <= 0.0
        || m_bond_break_strain <= 0.0 || m_bond_peak_strain <= 0.0
        || m_bond_peak_strain > m_bond_break_strain || m_reaction_weakening < 0.0
        || m_reaction_weakening > 1.0 || m_ignition_radius < 0.0
        || m_ignition_center_temperature < 0.0 || m_lattice_thickness <= 0.0
        || m_internal_heat_limit < 0.0 || m_internal_heat_limit > 1.0)
        throw ZaranError("DEM gas/reaction parameters are invalid");
    if (m_voronoi_update_interval <= 0 || m_voronoi_max_volume_ratio < 1.0
        || m_box_x_min >= m_box_x_max || m_box_y_min >= m_box_y_max
        || m_surface_energy < 0.0 || m_grain_boundary_energy < 0.0
        || m_grain_boundary_energy > m_surface_energy
        || m_high_pressure_rnn < 0.0 || m_high_pressure_rnn >= 1.0
        || m_high_pressure_knn < 0.0 || m_high_pressure_exponent < 0.0
        || m_normal_viscosity < 0.0 || m_tangential_viscosity < 0.0)
        throw ZaranError("DEM Voronoi/material parameters are invalid");
    if (m_spring_network_gap < 0.0 || m_spring_network_stiffness < 0.0
        || m_spring_network_tangential_stiffness < 0.0
        || m_spring_network_fracture_strain < 0.0)
        throw ZaranError("DEM spring network parameters are invalid");

    Log::info("DEM SolverParam Init: dt={:E}, contact_model={}, particle_file={}, spring_network={}",
              m_dt, m_contact_model, m_particle_file, m_spring_network_enabled ? "on" : "off");
}

} // namespace zaran
