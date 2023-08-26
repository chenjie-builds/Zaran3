#include"PerfectGas.h"
#include"MathBasic.h"
PerfectGas::PerfectGas(const double& Mw, const double& gamma, const RefValue& refValue) :Gas{ Mw, gamma, refValue }
{
	gamma_ = 1.4;
	double Rm = GetRm();
	Cv_ = Rm / (gamma_ - 1);
	Cp_ = gamma_ * Rm / (gamma_ - 1);
}
double PerfectGas::GetTemp(const double& density, const double& p)
{
	return p / (density * GetRm());
}

double PerfectGas::GetMul(const double& T)
{
	return Southerland(T);
}

double PerfectGas::GetMut(const double& T)
{
	return 0.0;//TODO：μt由湍流模型给定，层流等于0
}

double PerfectGas::GetMu(const double& T)
{
	return GetMul(T) + GetMut(T);//μ=μl+μt
}

double PerfectGas::GetKl(const double& T)
{
	return GetMul(T) * Cp_ / Prl;
}

double PerfectGas::GetKt(const double& T)
{

	return GetMut(T) * Cp_ / Prt;
}

double PerfectGas::GetK(const double& T)
{
	return GetKl(T) + GetKt(T);
}
double PerfectGas::GetEnergy(const double& T, const double& velocity)
{
	return GetRm() * T / (gamma_ - 1) + velocity * velocity / 2.0;
}
double PerfectGas::GetEnergy(const double& density, const double& pressure, const double& velocity)
{
	return pressure / ((gamma_ - 1) * density) + velocity * velocity / 2.0;
}