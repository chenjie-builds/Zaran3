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
class RefValue
{
public:
	RefValue();
	RefValue( const double&refDensity,const double &refGamma,const double&refMw,const double&refLength,const double&refTemp);
	void SetRefDensity(const double& refDensity) { density_ = refDensity; };
	void SetRefVelocity(const double& refVelocity) { velocity_ = refVelocity; };
	void SetRefPressure(const double& refPressure) { pressure_ = refPressure; };
	double GetRefGamma()const { return gama_; }
	double GetRefMw() const { return mw_; }
	double GetRefLength()const { return length_; }
	double GetRefDensity()const { return density_; }
	double GetRefTemp() const { return temp_; }
	double GetRefRm()const { return Rm_; }
	double GetRefVelocity() const { return velocity_; }
	double GetRefPressure() const { return pressure_; }
	double GetRefTime() const { return time_; }
	double GetRefMass() const { return mass_; }
	double GetRe() const { return Re_; }
private:
	double gama_;		//参考气体的比热比
	double mw_;			//参考气体的分子量(kg/mol)
	double length_;		//参考长度(m)
	double density_;	//参考密度(kg/m^3)
	double temp_;		//参考温度(K)
	double Rm_;			//气体常数的参考值(J/kg/K,K*m^2/s^2)
	double velocity_;	//参考速度(m/s)
	double pressure_;	//参考压力(Pa)
	double time_;		//参考时间(s)
	double mass_;		//参考质量(kg)
	double Re_;			//雷诺数(参考来流声速)
};