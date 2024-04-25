//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	Visual.h															||
//*	@brief	Á÷³¡ÏÔÊ¾																||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma	once
#include"FieldSolver.h"
namespace zaran
{
	class Visual
	{
	public:
		void WriteTecplot(FieldSolver* solver);
		void WriteTecplotBinary(FieldSolver* solver);
		void WriteTecplot2D(FieldSolver* solver);
		void WriteTecplotPoint(FieldSolver* solver);
		void WriteVTK(FieldSolver* solver);
		void WriteTecplotZaran3D(FieldSolver* solver);
		void WriteTecplotZaran3DBinary(FieldSolver* solver);
	};
}