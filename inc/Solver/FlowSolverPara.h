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
#include "solverpara.h"
#include "Limiter.h"
#include <vector>
#include <Eigen/Dense>
namespace zaran
{
	using namespace Eigen;
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
		void SetRefMachNumber(const int& refMach);
		void SetRefGamma(const int& refGamma);
		void SetRefDensity(const int& refDensity);
		void SetRefSonicSpeed(const int& refSonicSpped);
		void SetRefTemprature(const int& refTemperature);
		void SetInitFieldType(const int& initflowType);
		void SetIsViscous(const int& isViscous);
		void SetCflNumber(const int& cfl);
		void SetPrimitiveInflow(const VectorXd& primInflow);
		void SetRKCoef(const std::vector<double>& rkCoef);
		void SetGradScheme(const GradScheme& gradScheme);
		void SetLimiterType(const LimiterType& limiterType);
	public:
		const double& GetRefMachNumber()const;
		const double& GetRefGamma()const;
		const double& GetRefDensity()const;
		const double& GetRefSonicSpeed()const;
		const double& GetRefPressure()const;
		const double& GetRefTemperature()const;
		const int& GetInitFieldType()const;
		const int& GetIsViscous()const;
		const double& GetCflNumber()const;
		const VectorXd& GetPrimitiveInflow()const;
		const std::vector<double>& GetRKCoef()const;
		const GradScheme& GetGradScheme()const;
		const LimiterType& GetLimiterType()const;
	private:
		//referance flow parameter
		double refMachNumber_;
		double refGamma_;
		double refDensity_;
		double refSonicSpeed_;
		double refPressure_;
		double refTemperature_;
		// 流场初始化方式
		// =0 为静止流场
		// =1 为来流
		int initFieldType_;
		// 是否为粘性
		int isViscous_;
		// CFL数, 决定计算步长
		double cflNumber_;
		// 来流原始变量
		VectorXd primInflow_;
		// RK步数
		std::vector<double> rkCoef_;
		// 梯度方法
		GradScheme gradSchem_;
		// 限制器
		LimiterType limiterType_;
	};
}