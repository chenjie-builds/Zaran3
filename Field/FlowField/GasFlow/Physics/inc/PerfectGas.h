//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	PerfectGas.h														||
//*	@brief	完全气体类					||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include"Gas.h"
//理想气体类
class PerfectGas : public Gas
{
public:
	PerfectGas(const double& Mw = 0.028964, const double& gamma = 1.4, const Dimensionless& refValue = Dimensionless());
	double GetGamma() const override { return m_gamma; }
	double GetMw() const override { return m_Mw; }
	double GetRm() const override { return m_Rm; }
	double GetCp() const override { return m_cp; }
	double GetCv() const override { return m_cp / m_gamma; }
	double CalcSoundSpeed(const double& T)override;  //根据温度计算声速
	double CalcSoundSpeed(const double& density, const double& pressure) override;
	double CalcTemperature(const double& density, const double& p)override;  //根据状态方程计算温度
	double CalcMul(const double& T)override;   //根据温度，利用southerland公式计算层流粘性系数μl
	double CalcMut(const double& T)override;   //根据温度，计算湍流粘性系数μt
	double CalcMu(const double& T)override;    //根据温度求解粘性系数μ
	double CalcKl(const double& T)override;    //根据温度，求出层流传热系数
	double CalcKt(const double& T)override;    //根据温度，求出湍流传热系数
	double CalcK(const double& T)override;     //根据温度，求出传热系数
	double CalcEnergy(const double& T, const double& velocity)override;//根据温度和速度求出总内能
	double CalcEnergy(const double& density, const double& pressure, const double& velocity)override;//根据密度压力和速度求出总内能
	//根据密度和温度求出压力
	double CalcPressure(const double& density, const double& T)const override;
	void Prim2Cons(const double* prim, double* cons)override;
	void Cons2Prim(const double* cons, double* prim)override;
private:
	double m_gamma;	//比热比 
	double m_Mw;     //摩尔质量 kg/mol
	double m_Rm;	//无量纲气体常数
	double m_T0;	//无量纲Southerland公式中的温度T0
	double m_Ts;	//无量纲Southerland公式中的常数Ts
	double m_mu0;	//无量纲Southerland公式中的粘性系数Mu0
	double m_Prl;	//层流Prandtl数
	double m_Prt;	//湍流Prandtl数
	double m_cp;	//定压比热
};