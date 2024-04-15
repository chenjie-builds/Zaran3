//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	RefValue.h															||
//*	@brief	参考量类																||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
class Dimensionless
{
public:
	Dimensionless();
	Dimensionless(const double& refDensity, const double& refGamma, const double& refMw, const double& refLength, const double& refTemp);
	void SetRefValue(const double& refDensity, const double& refGamma, const double& refMw, const double& refLength, const double& refTemp);
	// 获取无量纲的参数（无量纲）
	double GetMwDL(double Mw) const { return Mw / m_Mw; }
	double GetDensityDL(double density) const { return density / m_density; }
	double GetTempDL(double temp) const { return temp / m_temperature; }
	double GetRmDL(double Rm) const { return Rm / m_Rm; }
	double GetVelocityDL(double velocity) const { return velocity / m_velocity; }
	double GetPressureDL(double pressure) const { return pressure / m_pressure; }
	double GetTimeDL(double time) const { return time / m_time; }
	double GetMassDL(double mass) const { return mass / m_mass; }
	// 获取真实的参数（有量纲）
	double GetMwReal(double MwDL) const { return MwDL * m_Mw; }
	double GetDensityReal(double densityDL) const { return densityDL * m_density; }
	double GetTempReal(double tempDL) const { return tempDL * m_temperature; }
	double GetRmReal(double RmDL) const { return RmDL * m_Rm; }
	double GetVelocityReal(double velocityDL) const { return velocityDL * m_velocity; }
	double GetPressureReal(double pressureDL) const { return pressureDL * m_pressure; }
	double GetTimeReal(double timeDL) const { return timeDL * m_time; }
	double GetMassReal(double massDL) const { return massDL * m_mass; }

	//  获取参考量（有量纲）
	double GetRefGamma()const { return m_gamma; }
	double GetRefMw() const { return m_Mw; }
	double GetRefLength()const { return m_length; }
	double GetRefDensity()const { return m_density; }
	double GetRefTemp() const { return m_temperature; }
	double GetRefRm()const { return m_Rm; }
	double GetRefVelocity() const { return m_velocity; }
	double GetRefPressure() const { return m_pressure; }
	double GetRefTime() const { return m_time; }
	double GetRefMass() const { return m_mass; }

	double GetRe() const { return m_Re; }



private:
	double m_gamma;		//参考气体的比热比
	double m_Mw;			//参考气体的分子量(kg/mol)
	double m_length;		//参考长度(m)
	double m_density;	//参考密度(kg/m^3)
	double m_temperature;		//参考温度(K)
	double m_Rm;			//气体常数的参考值(J/kg/K,K*m^2/s^2)
	double m_velocity;	//参考速度(m/s)
	double m_pressure;	//参考压力(Pa)
	double m_time;		//参考时间(s)
	double m_mass;		//参考质量(kg)
	double m_Re;			//雷诺数(参考来流声速)
};