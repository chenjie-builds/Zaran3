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

}

void FlowSolverPara::Init()
{
	double inflow_Ma = GlobalData::GetDouble("inflowMachNumber");
	double inflow_gamma = GlobalData::GetDouble("inflowGamma");
	double inflow_density = GlobalData::GetDouble("inflowDensity");
	double inflow_temperature = GlobalData::GetDouble("inflowTemperature");
	double inflow_Mw = GlobalData::GetDouble("inflowMw");
	double ref_length = GlobalData::GetDouble("refLength");
	m_dimensionless.SetRefValue(inflow_density, inflow_gamma, inflow_Mw, 1.0, inflow_temperature);
	double inflow_attack_angle = GlobalData::GetDouble("inflowAttackAngle");
	double inflow_slide_angle = GlobalData::GetDouble("inflowSlideAngle");
	PerfectGas gas(inflow_Mw, inflow_gamma, m_dimensionless);
	double inflow_sonice_speed = gas.GetSonicSpeed(m_dimensionless.GetTempDL(inflow_temperature));
	m_inflow_velocity_x = inflow_Ma * inflow_sonice_speed * cos(inflow_attack_angle) * cos(inflow_slide_angle);
	m_inflow_velocity_y = inflow_Ma * inflow_sonice_speed * sin(inflow_attack_angle) * cos(inflow_slide_angle);
	m_inflow_velocity_z = inflow_Ma * inflow_sonice_speed * sin(inflow_slide_angle);
	m_inflow_density =m_dimensionless.GetDensityDL(inflow_density);
	m_inflow_pressure =gas.GetPressureFromDensityAndTemperature(m_inflow_density, m_dimensionless.GetTempDL(inflow_temperature));
	initFieldType_ = GlobalData::GetInt("initFieldType");
	isViscous_ = GlobalData::GetInt("isViscous");
	cflNumber_ = GlobalData::GetDouble("cflNumber");
	int rkStage = GlobalData::GetInt("rkStage");
	if (rkStage == 1)
	{
		rkCoef_ = { 1.0 };
	}
	else if (rkStage == 3)
	{
		rkCoef_ = { 0.1918,0.4929,1.0 };
	}
	else if (rkStage == 4)
	{
		rkCoef_ = { 0.1084,0.2602,0.5052,1.0 };
	}
	else if (rkStage == 5)
	{
		rkCoef_ = { 0.0695,0.1602,0.2898,0.5060 };
	}
	else
	{
		ZaranLog::warn("Wrong Runge-Kutta Stage: {}, Please Check Control File!", rkStage);
		system("pause");
	}
	int gradScheme = GlobalData::GetInt("gradScheme");
	if (gradScheme == 0)
	{
		gradSchem_ = GradScheme::wls;
	}
	else if (gradScheme == 1)
	{
		gradSchem_ = GradScheme::ufdm;
	}
	else if (gradScheme == -1)
	{
		gradSchem_ = GradScheme::noGrad;
	}
	else
	{
		ZaranLog::warn("Unsupported Gradient Scheme parameter: {}, Please Check Control File!", gradScheme);
		system("pause");
	}

	std::string limiterType = GlobalData::GetString("limiterType");
	if (limiterType == "noLimiter")
	{
		limiterType_ = LimiterType::nolimit;
	}
	else if (limiterType == "barth")
	{
		limiterType_ = LimiterType::barth;
	}
	else if (limiterType == "vk")
	{
		limiterType_ = LimiterType::vk;
	}
	else
	{
		ZaranLog::warn("Unsupportted Limiter: {}, Please Check Control File!", limiterType);
		system("pause");
	}

}

const int& FlowSolverPara::GetInitFieldType() const
{
	return initFieldType_;
}

const int& FlowSolverPara::GetIsViscous() const
{
	return isViscous_;
}

const double& FlowSolverPara::GetCflNumber() const
{
	return cflNumber_;
}

const DArray& FlowSolverPara::GetRKCoef() const
{
	return rkCoef_;
}



void FlowSolverPara::SetInitFieldType(const int& initflowType)
{
	initFieldType_ = initflowType;
}

void FlowSolverPara::SetIsViscous(const int& isViscous)
{
	isViscous_ = isViscous;
}

void FlowSolverPara::SetCflNumber(const double& cfl)
{
	cflNumber_ = cfl;
}

void FlowSolverPara::SetRKCoef(const DArray& rkCoef)
{
	rkCoef_ = rkCoef;
}

const GradScheme& FlowSolverPara::GetGradScheme() const
{
	return gradSchem_;
}

const LimiterType& FlowSolverPara::GetLimiterType() const
{
	return limiterType_;
}

void FlowSolverPara::SetGradScheme(const GradScheme& gradScheme)
{
	gradSchem_ = gradScheme;
}

void FlowSolverPara::SetLimiterType(const LimiterType& limiterType)
{
	limiterType_ = limiterType;
}

const int& zaran::FlowSolverPara::GetInflowDensity() const
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
