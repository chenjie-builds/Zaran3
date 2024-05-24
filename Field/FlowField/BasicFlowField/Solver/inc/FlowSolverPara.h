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
#include"Gas.h"
#include "RiemannSolver.h"
namespace zaran
{
	//梯度求解方法
	enum class GradScheme
	{
		ufdm,//非结构有限差分法
		wls,//最小二乘法
		noGrad//不计算梯度
	};
	//流场初始化方式
	enum class InitFieldType
	{
		FarFlow,//远场
		FarFlowNoVelocity,//远场无速度
		Backup//备份文件
	};
	class FlowSolverPara :public SolverPara
	{
	public:
		FlowSolverPara();
		virtual ~FlowSolverPara();
		void Init() override;
	public:
		void SetInitFieldType(InitFieldType& initflowType);
		void SetIsViscous(const int& isViscous);
		void SetCflNumber(const double& cfl);
		void SetRKCoef(const DArray& rkCoef);
		void SetGradScheme(const GradScheme& gradScheme);
		void SetLimiterType(const LimiterType& limiterType);
		void SetBackupFieldFileName(const std::string& backupFieldFileName);
	public:
		const int& GetEquNum()const;
		const double& GetInflowDensity()const;
		const double& GetInflowVelocityX()const;
		const double& GetInflowVelocityY()const;
		const double& GetInflowVelocityZ()const;
		const double& GetInflowPressure()const;
		const double& GetInflowTemperature()const;
		const InitFieldType& GetInitFieldType()const;
		const int& GetIsViscous()const;
		const double& GetCflNumber()const;
		const DArray& GetRKCoef()const;
		const GradScheme& GetGradScheme()const;
		const LimiterType& GetLimiterType()const;
		const std::string& GetBackupFieldFileName()const;
		const Dimensionless& GetDimensionless()const;
		Gas* GetGas() { return m_gas; }
		const RiemannSolverType& GetRiemannSolverType()const;
	private:
		int m_equ_num;
		Gas* m_gas;
		// 无量纲参数
		Dimensionless m_dimensionless;
		double m_inflow_density;
		double m_inflow_velocity_x;
		double m_inflow_velocity_y;
		double m_inflow_velocity_z;
		double m_inflow_pressure;
		double m_inflow_temperature;
		// 流场初始化方式
		InitFieldType m_init_field_type;
		// 是否为粘性
		int m_is_viscous;
		// CFL数, 决定计算步长
		double m_cfl;
		// RK步数
		DArray m_rk_coef;
		// 梯度方法
		GradScheme m_grad_scheme;
		// 限制器
		LimiterType n_limiter_type;
		// 备份流场文件名，用于续算
		std::string m_backup_field_file_name;
		RiemannSolverType m_riemann_solver_type;
	};
}