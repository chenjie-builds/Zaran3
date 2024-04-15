//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	FlowSolverPara.h													||
//*	@brief	流场求解器参数类,定义了流场求解器正常进行的一些参数	 					||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "SolverPara.h"
#include "Limiter.h"
#include"BasicType.h"
#include"Dimensionless.h"
namespace zaran
{
	//梯度求解方法
	enum class GradScheme
	{
		ufdm,//非结构有限差分法
		wls,//最小二乘法
		noGrad//不计算梯度
	};
	class FlowSolverPara :public SolverPara
	{
	public:
		FlowSolverPara();
		virtual ~FlowSolverPara();
		void Init() override;
	public:
		void SetInitFieldType(const int& initflowType);
		void SetIsViscous(const int& isViscous);
		void SetCflNumber(const double& cfl);
		void SetRKCoef(const DArray& rkCoef);
		void SetGradScheme(const GradScheme& gradScheme);
		void SetLimiterType(const LimiterType& limiterType);
	public:
		const int& GetInflowDensity()const;
		const double& GetInflowVelocityX()const;
		const double& GetInflowVelocityY()const;
		const double& GetInflowVelocityZ()const;
		const double& GetInflowPressure()const;
		const double& GetInflowTemperature()const;
		const int& GetInitFieldType()const;
		const int& GetIsViscous()const;
		const double& GetCflNumber()const;
		const DArray& GetRKCoef()const;
		const GradScheme& GetGradScheme()const;
		const LimiterType& GetLimiterType()const;
	private:
		// 无量纲参数
		Dimensionless m_dimensionless;
		double m_inflow_density;
		double m_inflow_velocity_x;
		double m_inflow_velocity_y;
		double m_inflow_velocity_z;
		double m_inflow_pressure;
		double m_inflow_temperature;
		// 流场初始化方式
		// =0 为静止流场
		// =1 为来流
		int initFieldType_;
		// 是否为粘性
		int isViscous_;
		// CFL数, 决定计算步长
		double cflNumber_;
		// RK步数
		DArray rkCoef_;
		// 梯度方法
		GradScheme gradSchem_;
		// 限制器
		LimiterType limiterType_;
	};
}