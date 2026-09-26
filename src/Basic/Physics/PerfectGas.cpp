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
	// 无量纲 Sutherland：μ*(T*) = μ0* · (T*/T0*)^1.5 · (T0*+Ts*)/(T*+Ts*)
	//
	// ⚠ 原来的写法是 `m_mu0 * Southerland(T, m_T0, m_Ts)`，但 4 参数 Southerland 的
	//   第 2 个形参是 **mu0** 而不是 T0（`Southerland(T, mu0, T0, Ts)`）。于是
	//   m_T0 = 273.16/T_ref 被当成 μ0、m_Ts = 110.4/T_ref 被当成 T0 用了，
	//   算出来的 μ 比真值大 5~6 倍（实测：空气 439 K 时给出 1.35e-4 而不是 2.43e-5）。
	//   该函数此前**没有任何调用点**（粘性残差在现有算例里都没开），所以从未暴露；
	//   阶段 3 的曳力需要 Re = ρ|u|d/μ，才把它逼出来。
	//   修好后 T* = 1 处应精确还原 1.716e-5 Pa·s（见 tests/demcfd_shock_tube/verify.py 的 V11）。
	return Southerland(T, m_mu0, m_T0, m_Ts);
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

void PerfectGas::Prim2Cons(const double* prim, double* cons)
{
	double v2 = prim[1] * prim[1] + prim[2] * prim[2] + prim[3] * prim[3];
	cons[0] = prim[0];
	cons[1] = prim[0] * prim[1];
	cons[2] = prim[0] * prim[2];
	cons[3] = prim[0] * prim[3];
	cons[4] = 0.5 * prim[0] * v2 + prim[4] / (m_gamma - 1);
}

void PerfectGas::Cons2Prim(const double* cons, double* prim)
{
	prim[0] = cons[0];
	prim[1] = cons[1] / cons[0];
	prim[2] = cons[2] / cons[0];
	prim[3] = cons[3] / cons[0];
	double v2 = prim[1] * prim[1] + prim[2] * prim[2] + prim[3] * prim[3];
	prim[4] = (cons[4] - 0.5 * cons[0] * v2) * (m_gamma - 1);
}
