//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Gas.h																||
//*	@brief	气体基类																||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include"CommonPara.h"
#include"Dimensionless.h"
const double GAS_CONSTANT = 8.3143;//通用气体常数J/(mol*k)
enum class FluidType
{
	PerfectGas,
	RealGas,
};
//表示气体的抽象基类，当前只有基本功能
//后续还有粘性，分子量等功能加入
class Gas
{
public:
	Gas(const double& Mw=0.028964, const double& gamma=1.4, const Dimensionless& refValue = Dimensionless());

	double GetRm()const { return Rm_; }
	double GetMw() const { return mw_; }
	virtual double GetGamma() const = 0;
	virtual double GetSonicSpeed(const double& T)=0;  //根据温度计算声速
	virtual double GetTemp(const double& density, const double& p)=0;  //根据状态方程计算温度
	virtual double GetMul(const double& T)=0;   //根据温度，利用southerland公式计算层流粘性系数μl
	virtual double GetMut(const double& T)=0;   //根据温度，计算湍流粘性系数μt
	virtual double GetMu(const double& T)=0;    //根据温度求解粘性系数μ
	virtual double GetKl(const double& T)=0;    //根据温度，求出层流传热系数
	virtual double GetKt(const double& T)=0;    //根据温度，求出湍流传热系数
	virtual double GetK(const double& T)=0;     //根据温度，求出传热系数
	virtual double GetCp(const double& T) = 0;	//根据温度，求出定压比热Cp
	virtual double GetCv(const double& T) = 0;	//根据温度，求出定容比热Cv
	virtual double GetEnergy(const double& T, const double& velocity) = 0;//根据温度和速度求出总内能
	virtual double GetEnergy(const double& density, const double& pressure, const double& velocity) = 0;//根据密度压力和速度求出总内能
private:
	double mw_;     //摩尔质量 kg/mol
	double Rm_;		//气体常数
};
