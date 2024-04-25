#include "Dimensionless.h"
#include"CommonPara.h"
#include"MathBasic.h"
#include "Gas.h"
#include<cmath>
using namespace zaran;
Dimensionless::Dimensionless()
{
	m_gamma = 1.4;
	m_Mw = 0.02896;
	m_length = 1.0;
	m_density = 1.0;
	m_temperature = 2.488917144919069e+003;
	m_Rm = m_gamma * GAS_CONSTANT / m_Mw;
	m_velocity = sqrt(m_gamma * GAS_CONSTANT / m_Mw * m_temperature);//当地声速
	m_pressure = m_density * m_velocity * m_velocity;
	m_time = m_length / m_velocity;
	m_mass = m_density * m_length * m_length * m_length;
	m_Re = m_density * m_velocity * m_length / Southerland(m_temperature);
	double mu0 = 1.711e-5 / m_density / m_velocity / m_length;
	double T0 = 273.0 / m_temperature;
	double Ts = 110.4 / m_temperature;
	m_Re = 1.0 / Southerland(1.0, mu0, T0, Ts) * m_length;
}
Dimensionless::Dimensionless(const double& refDensity, const double& refGamma, const double& refMw, const double& refLength, const double& refTemp)
{
	SetRefValue(refDensity, refGamma, refMw, refLength, refTemp);

}

void Dimensionless::SetRefValue(const double& refDensity, const double& refGamma, const double& refMw, const double& refLength, const double& refTemp)
{
	m_gamma = refGamma;
	m_Mw = refMw;
	m_length = refLength;
	m_density = refDensity;
	m_temperature = refTemp;
	m_velocity = sqrt(m_gamma * GAS_CONSTANT / m_Mw * m_temperature);//当地声速
	m_Rm = m_velocity * m_velocity / m_temperature;
	m_pressure = m_density * m_velocity * m_velocity;
	m_time = m_length / m_velocity;
	m_mass = m_density * m_length * m_length * m_length;
	m_Re = m_density * m_velocity * m_length / Southerland(m_temperature);
	double mu0 = 1.711e-5 / m_density / m_velocity / m_length;
	double T0 = 273.0 / m_temperature;
	double Ts = 110.4 / m_temperature;
	m_Re = m_density * m_velocity * m_length / Southerland(1.0, mu0, T0, Ts);
}
