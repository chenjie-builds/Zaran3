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
		void WriteTecplot(Ptr<FieldSolver>& solver);
		void WriteTecplotBinary(Ptr<FieldSolver>& solver);
		void WriteTecplot2D(Ptr<FieldSolver>& solver);
		void WriteTecplotPoint(Ptr<FieldSolver>& solver);
		void WriteVTK(Ptr<FieldSolver>& solver);
		void WriteTecplotZaran3D(Ptr<FieldSolver>& solver);
		void WriteTecplotZaran3DBinary(Ptr<FieldSolver>& solver);
	};
}