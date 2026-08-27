#include "FlowSolverPara.h"
#include "GlobalData.h"
#include "PerfectGas.h"
#include "Log.h"
#include "ZaranError.h"
#include <iostream>
using namespace zaran;
FlowSolverParam::FlowSolverParam()
{
}

FlowSolverParam::~FlowSolverParam()
{
	if (m_gas)
	{
		delete m_gas;
		m_gas = nullptr;
	}
	if (m_cfl)
	{
		delete m_cfl;
		m_cfl = nullptr;
	}
}

void FlowSolverParam::Init()
{
	m_equ_num = 5;
	double inflow_Ma = GlobalData::GetDouble("freestream.mach_number");
	double inflow_gamma = GlobalData::GetDouble("freestream.gamma");
	double inflow_density = GlobalData::GetDouble("freestream.density");
	double inflow_temperature = GlobalData::GetDouble("freestream.temperature");
	double inflow_Mw = GlobalData::GetDouble("freestream.mw");
	double ref_length = GlobalData::GetDouble("freestream.ref_length");
	m_dimensionless.SetRefValue(inflow_density, inflow_gamma, inflow_Mw, ref_length, inflow_temperature);
	double inflow_attack_angle = GlobalData::GetDouble("freestream.attack_angle");
	double inflow_slide_angle = GlobalData::GetDouble("freestream.slide_angle");
	m_gas = new PerfectGas(inflow_Mw, inflow_gamma, m_dimensionless);
	double inflow_sonic_speed = m_gas->CalcSoundSpeed(m_dimensionless.GetTempDL(inflow_temperature));
	m_inflow_velocity_x = inflow_Ma * inflow_sonic_speed * cos(inflow_attack_angle) * cos(inflow_slide_angle);
	m_inflow_velocity_y = inflow_Ma * inflow_sonic_speed * sin(inflow_attack_angle) * cos(inflow_slide_angle);
	m_inflow_velocity_z = inflow_Ma * inflow_sonic_speed * sin(inflow_slide_angle);
	m_inflow_density = m_dimensionless.GetDensityDL(inflow_density);
	m_inflow_pressure = m_gas->CalcPressure(m_inflow_density, m_dimensionless.GetTempDL(inflow_temperature));
	string inflow_type = GlobalData::GetString("init.field_type");
	if (inflow_type == "FarFlow")
	{
		m_init_field_type = InitFieldType::FarFlow;
	}
	else if (inflow_type == "FarFieldNoVelocity")
	{
		m_init_field_type = InitFieldType::FarFlowNoVelocity;
	}
	else if (inflow_type == "Backup")
	{
		m_init_field_type = InitFieldType::Backup;
	}
	else if (inflow_type == "Explosion")
	{
		m_init_field_type = InitFieldType::Explosion;
	}
	else if (inflow_type == "Vortex")
	{
		m_init_field_type = InitFieldType::Vortex;
	}
	else
	{
        Log::warn("Unsupported Init Field Type: {}, Please Check Control File!", inflow_type);
		m_init_field_type = InitFieldType::FarFlow;
	}
	m_is_viscous = GlobalData::GetBool("task.is_viscous") ? 1 : 0;
	Log::info("Inflow Parameters: Density: {}, Velocity: ({},{},{}), Pressure: {}", m_inflow_density, m_inflow_velocity_x, m_inflow_velocity_y, m_inflow_velocity_z, m_inflow_pressure);
	int rkStage = GlobalData::GetInt("task.rk_stage");
	if (rkStage == 1)
	{
		m_rk_coef.resize(1);
		m_rk_coef[0] = { 1.0,0.0,1.0 };
	}
	else if (rkStage == 2)
	{
		m_rk_coef.resize(2);
		m_rk_coef[0] = { 1.0,0.0,1.0 };
		m_rk_coef[1] = { 0.5,0.5,0.5 };
	}
	else if (rkStage == 3)
	{
		m_rk_coef.resize(3);
		m_rk_coef[0] = { 1.0,0.0,1.0 };
		m_rk_coef[1] = { 3.0 / 4.0,1.0 / 4.0,1.0 / 4.0 };
		m_rk_coef[2] = { 1.0 / 3.0,2.0 / 3.0,2.0 / 3.0 };
	}
	else
	{
		Log::warn("Wrong Runge-Kutta Stage: {}, Please Check Control File!", rkStage);
		throw ZaranError("Wrong Runge-Kutta Stage");
	}
	int gradScheme = GlobalData::GetInt("task.grad_scheme");
	if (gradScheme == 0)
	{
		m_grad_scheme = GradScheme::WLS;
	}
	else if (gradScheme == 1)
	{
		m_grad_scheme = GradScheme::UFDM;
	}
	else if (gradScheme == -1)
	{
		m_grad_scheme = GradScheme::NONE;
	}
	else
	{
		Log::warn("Unsupported Gradient Scheme parameter: {}, Please Check Control File!", gradScheme);
		throw ZaranError("Unsupported Gradient Scheme");
	}

	m_backup_field_file_name = GlobalData::GetString("init.backup_file");
	std::string riemann_solver_type = GlobalData::GetString("flow.riemann_solver");
	if (riemann_solver_type == "HLLC")
	{
		m_riemann_solver_type = RiemannSolverType::HLLC;
	}
	else if (riemann_solver_type == "VanLeer")
	{
		m_riemann_solver_type = RiemannSolverType::VanLeer;
	}
	else if (riemann_solver_type == "Roe")
	{
		m_riemann_solver_type = RiemannSolverType::Roe;
	}
	else if (riemann_solver_type == "StegerWarming")
	{
		m_riemann_solver_type = RiemannSolverType::StegerWarming;
	}
	else if (riemann_solver_type == "Ausmpw")
	{
		m_riemann_solver_type = RiemannSolverType::Ausmpw;
	}
	else
	{
		Log::warn("Unsupported Riemann Solver Type: {}, Please Check Control File!", riemann_solver_type);
		throw ZaranError("Unsupported Riemann Solver Type");
	}
	InitCflNumber();
	InitLimiter();
}

void zaran::FlowSolverParam::InitCflNumber()
{
	double cfl_min = GlobalData::GetDouble("cfl.min");
	double cfl_max = GlobalData::GetDouble("cfl.max");
	int start_step = GlobalData::GetInt("iteration.current_iter");
	int grow_step = GlobalData::GetInt("cfl.grow_step");
	double reduce_factor = GlobalData::GetDouble("cfl.reduce_factor");
	m_cfl = new CFL(cfl_min, cfl_max, start_step, grow_step, reduce_factor);
}

void zaran::FlowSolverParam::InitLimiter()
{

	std::string limiter_type = GlobalData::GetString("space.limiter");
	if (limiter_type == "noLimiter")
	{
		n_limiter_type = LimiterType::none;
	}
	else if (limiter_type == "barth")
	{
		n_limiter_type = LimiterType::barth;
	}
	else if (limiter_type == "vk")
	{
		n_limiter_type = LimiterType::vk;
	}
	else if (limiter_type == "1st-order")
	{
		n_limiter_type = LimiterType::first_order;
	}
	else
	{
		Log::warn("Unsupportted Limiter: {}, Please Check Control File!", limiter_type);
		throw ZaranError("Unsupported Limiter");
	}
	m_max_limit = GlobalData::GetDouble("space.max_limiter_slope");
	if (m_max_limit < 0.0 || m_max_limit > 1.0)
	{
		Log::warn("Wrong Max Limiter Slope: {}, Please Check Control File!", m_max_limit);
		throw ZaranError("Wrong Max Limiter Slope");
	}
}

const InitFieldType& FlowSolverParam::GetInitFieldType() const
{
	return m_init_field_type;
}

int FlowSolverParam::GetIsViscous() const
{
	return m_is_viscous;
}

int zaran::FlowSolverParam::GetCurrentStep() const
{
	return m_current_step;
}

double zaran::FlowSolverParam::GetCflNumber() const
{
	return m_cfl->GetCFL(m_current_step);
}

void zaran::FlowSolverParam::ReduceCflNumber()
{
	m_cfl->ReduceCFL();
}
int zaran::FlowSolverParam::GetRkStep() const
{
	return m_rk_coef.size();
}
const dynamic_array<double>& FlowSolverParam::GetRkCoef(int iStage)const 
{
	return m_rk_coef[iStage];
}

void FlowSolverParam::SetInitFieldType(InitFieldType& initflowType)
{
	m_init_field_type = initflowType;
}

void FlowSolverParam::SetIsViscous(const int& isViscous)
{
	m_is_viscous = isViscous;
}


const GradScheme& FlowSolverParam::GetGradScheme() const
{
	return m_grad_scheme;
}

const LimiterType& FlowSolverParam::GetLimiterType() const
{
	return n_limiter_type;
}

const std::string& zaran::FlowSolverParam::GetBackupFieldFileName() const
{
	return m_backup_field_file_name;
}

const Dimensionless& zaran::FlowSolverParam::GetDimensionless() const
{
	return m_dimensionless;
}

const RiemannSolverType& zaran::FlowSolverParam::GetRiemannSolverType() const
{
	return m_riemann_solver_type;
}

const double zaran::FlowSolverParam::GetMaxLimit() const
{
	return m_max_limit;
}

void FlowSolverParam::SetGradScheme(const GradScheme& gradScheme)
{
	m_grad_scheme = gradScheme;
}

void FlowSolverParam::SetLimiterType(const LimiterType& limiterType)
{
	n_limiter_type = limiterType;
}

void zaran::FlowSolverParam::SetBackupFieldFileName(const std::string& backupFieldFileName)
{
	m_backup_field_file_name = backupFieldFileName;
}

void zaran::FlowSolverParam::SetCurrentStep(const int& currentStep)
{
	m_current_step = currentStep;
}

const int& zaran::FlowSolverParam::GetEqNum() const
{
	return m_equ_num;
}

const double& zaran::FlowSolverParam::GetInflowDensity() const
{
	return m_inflow_density;
}

const double& zaran::FlowSolverParam::GetInflowVelocityX() const
{
	return m_inflow_velocity_x;
}

const double& zaran::FlowSolverParam::GetInflowVelocityY() const
{
	return m_inflow_velocity_y;
}

const double& zaran::FlowSolverParam::GetInflowVelocityZ() const
{
	return m_inflow_velocity_z;
}

const double& zaran::FlowSolverParam::GetInflowPressure() const
{
	return m_inflow_pressure;
}

const double& zaran::FlowSolverParam::GetInflowTemperature() const
{
	return m_inflow_temperature;
}
