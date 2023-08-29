//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	SolverPara.h														||
//*	@brief	求解器参数基类														||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
namespace zaran
{
	/*
	求解器对应的参数类
	*/
	class SolverPara
	{
	public:
		SolverPara();
		virtual ~SolverPara();
		virtual void Init();
	private:
	};


}