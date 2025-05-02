/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file SolverParam.h
 * \brief SolverParam class, used to store solver parameters.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
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