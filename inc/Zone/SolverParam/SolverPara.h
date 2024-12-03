//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	SolverParam.h														||
//*	@brief	求解器参数基类														||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
namespace zaran
{
	/*
	求解器对应的参数类
	*/
	class SolverParam
	{
	public:
		SolverParam();
		virtual ~SolverParam();
		virtual void Init();
	private:
	};


}