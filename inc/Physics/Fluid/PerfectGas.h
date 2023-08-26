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
//空气
const double Prl = 0.72;
const double Prt = 0.9;

//理想气体类
class PerfectGas : public Gas
{
public:
	PerfectGas(const double& Mw = 0.028964, const double& gamma = 1.4, const RefValue& refValue = RefValue());


	double GetGamma() const { return gamma_; }
	double GetTemp(const double& density, const double& p);  //根据状态方程计算温度
	double GetMul(const double& T);   //根据温度，利用southerland公式计算层流粘性系数μl
	double GetMut(const double& T);   //根据温度，计算湍流粘性系数μt
	double GetMu(const double& T);    //根据温度求解粘性系数μ
	double GetKl(const double& T);    //根据温度，求出层流传热系数
	double GetKt(const double& T);    //根据温度，求出湍流传热系数
	double GetK(const double& T);     //根据温度，求出传热系数
	double GetCp(const double& T) { return Cp_; };	//根据温度，求出定压比热Cp
	double GetCv(const double& T) { return Cv_; };	//根据温度，求出定容比热Cv
	virtual double GetEnergy(const double& T, const double& velocity);//根据温度和速度求出总内能
	virtual double GetEnergy(const double& density, const double& pressure, const double& velocity);//根据密度压力和速度求出总内能
private:
	double gamma_;	//比热比 
	double Cp_;		//定压比热	J(mol*k)
	double Cv_;		//定容比热	J(mol*k)
};