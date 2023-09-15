#include"flowsolverpara.h"
#include "GlobalData.h"
#include "log.h"
using namespace zaran;
FlowSolverPara::FlowSolverPara()
{

}

FlowSolverPara::~FlowSolverPara()
{

}

void FlowSolverPara::Init()
{
	refMachNumber_ = GlobalData::GetDouble("refMachNumber");
	refGamma_ = GlobalData::GetDouble("refGamma");
	refDensity_ = GlobalData::GetDouble("refDensity");
	refSonicSpeed_ = GlobalData::GetDouble("refSonicSpeed");
	refPressure_ = GlobalData::GetDouble("refPressure");
	refTemperature_ = GlobalData::GetDouble("refTemperature");
	initFieldType_ = GlobalData::GetInt("initFieldType");
	isViscous_ = GlobalData::GetInt("isViscous");
	cflNumber_ = GlobalData::GetDouble("cflNumber");
	primInflow_.resize(5);
	primInflow_ << 1.4, 3, 0, 0, 1.0;
	primInflow_[1] = refMachNumber_;
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
	if (limiterType == "minmod")
	{
		limiterType_ = LimiterType::minmod;
	}
	else if (limiterType == "vanleer")
	{
		limiterType_ = LimiterType::vanleer;
	}
	else if (limiterType == "vanalbada")
	{
		limiterType_ = LimiterType::vanalbada;
	}
	else if (limiterType == "nolimit")
	{
		limiterType_ = LimiterType::nolimit;
	}
	else if (limiterType == "oneorder")
	{
		limiterType_ = LimiterType::oneorder;
	}
	else if (limiterType == "barth")
	{
		limiterType_ = LimiterType::barth;
	}
	else
	{
		ZaranLog::warn("Unsupportted Limiter: {}, Please Check Control File!", limiterType);
		system("pause");
	}

}

const double& FlowSolverPara::GetRefMachNumber() const
{
	return refMachNumber_;
}

const double& FlowSolverPara::GetRefGamma() const
{
	return refGamma_;
}

const double& FlowSolverPara::GetRefDensity() const
{
	return refDensity_;
}

const double& FlowSolverPara::GetRefSonicSpeed() const
{
	return refSonicSpeed_;
}

const double& FlowSolverPara::GetRefPressure() const
{
	return refPressure_;
}

const double& FlowSolverPara::GetRefTemperature() const
{
	return refTemperature_;
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

const DVector& FlowSolverPara::GetPrimitiveInflow() const
{
	return primInflow_;
}

const DArray& FlowSolverPara::GetRKCoef() const
{
	return rkCoef_;
}

void FlowSolverPara::SetRefMachNumber(const int& refMach)
{
	refMachNumber_ = refMach;

}

void FlowSolverPara::SetRefGamma(const int& refGamma)
{
	refGamma_ = refGamma;
}

void FlowSolverPara::SetRefDensity(const int& refDensity)
{
	refDensity_ = refDensity;
}

void FlowSolverPara::SetRefSonicSpeed(const int& refSonicSpped)
{
	refSonicSpeed_ = refSonicSpped;
}

void FlowSolverPara::SetRefTemprature(const int& refTemperature)
{
	refTemperature_ = refTemperature;
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

void FlowSolverPara::SetPrimitiveInflow(const DVector& primInflow)
{
	primInflow_ = primInflow;
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
