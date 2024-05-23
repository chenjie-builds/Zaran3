#include"flowsolverpara.h"
#include "GlobalData.h"
#include "log.h"
#include "PerfectGas.h"
using namespace zaran;
FlowSolverPara::FlowSolverPara()
{

}

FlowSolverPara::~FlowSolverPara()
{
	if (m_gas)
	{
		delete m_gas;
		m_gas = nullptr;
	}
	Log::info("FlowSolverPara is destroyed!");
	exit(0);
}

void FlowSolverPara::Init()
{
	m_equ_num = 5;
	double inflow_Ma = GlobalData::GetDouble("inflowMachNumber");
	double inflow_gamma = GlobalData::GetDouble("inflowGamma");
	double inflow_density = GlobalData::GetDouble("inflowDensity");
	double inflow_temperature = GlobalData::GetDouble("inflowTemperature");
	double inflow_Mw = GlobalData::GetDouble("inflowMw");
	double ref_length = GlobalData::GetDouble("refLength");
	m_dimensionless.SetRefValue(inflow_density, inflow_gamma, inflow_Mw, ref_length, inflow_temperature);
	double inflow_attack_angle = GlobalData::GetDouble("inflowAttackAngle");
	double inflow_slide_angle = GlobalData::GetDouble("inflowSlideAngle");
	m_gas = new PerfectGas(inflow_Mw, inflow_gamma, m_dimensionless);
	double inflow_sonic_speed = m_gas->CalcSoundSpeed(m_dimensionless.GetTempDL(inflow_temperature));
	m_inflow_velocity_x = inflow_Ma * inflow_sonic_speed * cos(inflow_attack_angle) * cos(inflow_slide_angle);
	m_inflow_velocity_y = inflow_Ma * inflow_sonic_speed * sin(inflow_attack_angle) * cos(inflow_slide_angle);
	m_inflow_velocity_z = inflow_Ma * inflow_sonic_speed * sin(inflow_slide_angle);
	m_inflow_density = m_dimensionless.GetDensityDL(inflow_density);
	m_inflow_pressure = m_gas->CalcPressure(m_inflow_density, m_dimensionless.GetTempDL(inflow_temperature));
	string inflow_type = GlobalData::GetString("initFieldType");
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
	else
	{
		m_init_field_type = InitFieldType::FarFlow;
	}
	m_is_viscous = GlobalData::GetInt("isViscous");
	m_cfl = GlobalData::GetDouble("cflNumber");
	Log::info("Inflow Parameters: Density: {}, Velocity: ({},{},{}), Pressure: {}", m_inflow_density, m_inflow_velocity_x, m_inflow_velocity_y, m_inflow_velocity_z, m_inflow_pressure);
	int rkStage = GlobalData::GetInt("rkStage");
	if (rkStage == 1)
	{
		m_rk_coef = { 1.0 };
	}
	else if (rkStage == 3)
	{
		m_rk_coef = { 0.1918,0.4929,1.0 };
	}
	else if (rkStage == 4)
	{
		m_rk_coef = { 0.1084,0.2602,0.5052,1.0 };
	}
	else if (rkStage == 5)
	{
		m_rk_coef = { 0.0695,0.1602,0.2898,0.5060 };
	}
	else
	{
		Log::warn("Wrong Runge-Kutta Stage: {}, Please Check Control File!", rkStage);
		system("pause");
	}
	int gradScheme = GlobalData::GetInt("gradScheme");
	if (gradScheme == 0)
	{
		m_grad_scheme = GradScheme::wls;
	}
	else if (gradScheme == 1)
	{
		m_grad_scheme = GradScheme::ufdm;
	}
	else if (gradScheme == -1)
	{
		m_grad_scheme = GradScheme::noGrad;
	}
	else
	{
		Log::warn("Unsupported Gradient Scheme parameter: {}, Please Check Control File!", gradScheme);
		system("pause");
	}

	std::string limiter_type = GlobalData::GetString("limiterType");
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
		system("pause");
	}
	m_backup_field_file_name = GlobalData::GetString("backupFieldFileName");

}

const InitFieldType& FlowSolverPara::GetInitFieldType() const
{
	return m_init_field_type;
}

const int& FlowSolverPara::GetIsViscous() const
{
	return m_is_viscous;
}

const double& FlowSolverPara::GetCflNumber() const
{
	return m_cfl;
}

const DArray& FlowSolverPara::GetRKCoef() const
{
	return m_rk_coef;
}



void FlowSolverPara::SetInitFieldType(InitFieldType& initflowType)
{
	m_init_field_type = initflowType;
}

void FlowSolverPara::SetIsViscous(const int& isViscous)
{
	m_is_viscous = isViscous;
}

void FlowSolverPara::SetCflNumber(const double& cfl)
{
	m_cfl = cfl;
}

void FlowSolverPara::SetRKCoef(const DArray& rkCoef)
{
	m_rk_coef = rkCoef;
}

const GradScheme& FlowSolverPara::GetGradScheme() const
{
	return m_grad_scheme;
}

const LimiterType& FlowSolverPara::GetLimiterType() const
{
	return n_limiter_type;
}

const std::string& zaran::FlowSolverPara::GetBackupFieldFileName() const
{
	return m_backup_field_file_name;
}

const Dimensionless& zaran::FlowSolverPara::GetDimensionless() const
{
	return m_dimensionless;
}

void FlowSolverPara::SetGradScheme(const GradScheme& gradScheme)
{
	m_grad_scheme = gradScheme;
}

void FlowSolverPara::SetLimiterType(const LimiterType& limiterType)
{
	n_limiter_type = limiterType;
}

void zaran::FlowSolverPara::SetBackupFieldFileName(const std::string& backupFieldFileName)
{
	m_backup_field_file_name = backupFieldFileName;
}

const int& zaran::FlowSolverPara::GetEquNum() const
{
	return m_equ_num;

}

const double& zaran::FlowSolverPara::GetInflowDensity() const
{
	return m_inflow_density;
}

const double& zaran::FlowSolverPara::GetInflowVelocityX() const
{
	return m_inflow_velocity_x;
}

const double& zaran::FlowSolverPara::GetInflowVelocityY() const
{
	return m_inflow_velocity_y;
}

const double& zaran::FlowSolverPara::GetInflowVelocityZ() const
{
	return m_inflow_velocity_z;
}

const double& zaran::FlowSolverPara::GetInflowPressure() const
{
	return m_inflow_pressure;
}

const double& zaran::FlowSolverPara::GetInflowTemperature() const
{
	return m_inflow_temperature;
}
