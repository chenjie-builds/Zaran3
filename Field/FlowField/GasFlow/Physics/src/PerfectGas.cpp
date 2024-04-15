#include"PerfectGas.h"
#include"MathBasic.h"
using namespace zaran;
PerfectGas::PerfectGas(const double& Mw, const double& gamma, const Dimensionless& refValue) :Gas{ Mw, gamma, refValue }
{
	gamma_ = 1.4;
	double Rm = GetRm();
	Cv_ = Rm / (gamma_ - 1);
	Cp_ = gamma_ * Rm / (gamma_ - 1);
}
double PerfectGas::GetSonicSpeed(const double& T)
{
	return sqrt(gamma_ * T);
}
double PerfectGas::GetTemp(const double& density, const double& p)
{
	return gamma_ * p / density;
}

double PerfectGas::GetMul(const double& T)
{
	return 0.0;//TODO：μl由southerland公式给定，暂时等于0
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

double PerfectGas::GetPressureFromDensityAndTemperature(const double& density, const double& T)
{
	return density * T / gamma_;
}
