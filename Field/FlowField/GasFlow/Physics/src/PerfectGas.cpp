#include"PerfectGas.h"
#include"MathBasic.h"
using namespace zaran;
PerfectGas::PerfectGas(const double& Mw, const double& gamma, const Dimensionless& refValue) :Gas{ Mw, gamma, refValue }
{
	m_Mw = Mw;
	m_Rm = GAS_CONSTANT / m_Mw / refValue.GetRefRm();
	m_gamma = gamma;
	m_T0 = 273.15 / refValue.GetRefTemp();
	m_Ts = 110.4 / refValue.GetRefTemp();
	m_mu0 = 1.716e-5 / (refValue.GetRefDensity() * refValue.GetRefVelocity() * refValue.GetRefLength());
	m_Prl = 0.72;
	m_Prt = 0.9;
	m_cp = m_gamma / (m_gamma - 1.0) * m_Mw;
}
double PerfectGas::CalcSoundSpeed(const double& T)
{
	return sqrt(T);
}
double PerfectGas::CalcSoundSpeed(const double& density, const double& pressure)
{
	return sqrt(GetGamma() * pressure / density);
}
double PerfectGas::CalcTemperature(const double& density, const double& p)
{
	return GetGamma() * p / density;
}

double PerfectGas::CalcMul(const double& T)
{
	return m_mu0 * Southerland(T, m_T0, m_Ts);
}

double PerfectGas::CalcMut(const double& T)
{
	return 0.0;//TODO：μt由湍流模型给定，层流等于0
}

double PerfectGas::CalcMu(const double& T)
{
	return CalcMul(T) + CalcMut(T);//μ=μl+μt
}

double PerfectGas::CalcKl(const double& T)
{
	return CalcMul(T) * GetCp() / m_Prl;
}

double PerfectGas::CalcKt(const double& T)
{

	return CalcMut(T) * GetCp() / m_Prt;
}

double PerfectGas::CalcK(const double& T)
{
	return CalcKl(T) + CalcKt(T);
}
double PerfectGas::CalcEnergy(const double& T, const double& velocity)
{
	return GetRm() * T / (m_gamma - 1) + velocity * velocity / 2.0;
}
double PerfectGas::CalcEnergy(const double& density, const double& pressure, const double& velocity)
{
	return pressure / ((m_gamma - 1) * density) + velocity * velocity / 2.0;
}

double PerfectGas::CalcPressure(const double& density, const double& T)const
{
	return density * T / m_gamma;
}
