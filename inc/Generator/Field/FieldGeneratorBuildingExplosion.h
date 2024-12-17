//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software
//|| *	Copyright (C) ,Since 2020
//||
//*-----------------------------------------------------------------------------||
//*	License
//|| *	This file is part of ZaRan.
//|| *
//|| *	@file	FieldGeneratorBuildingExplosion.h
//|| *	@brief	求解器工厂类, 生成求解器
//|| *	@author	Chen Jie.
//||
//==============================================================================||
#pragma once
#include "FieldGenerator.h"
#include "GridBlock.h"
#include "GridGeneratorBlock.h"
#include "NSFieldBlockExplosion.h"
namespace zaran
{
	class FieldGeneratorBuildingExplosion : public FieldGenerator
	{
	public:
		FieldGeneratorBuildingExplosion(GridType grid_type, FieldSolverType solver_type, Dimension dim);
		shared_ptr<FieldManager>  Create() override;
	private:
		void CreateGrid();
		void TagGrid();
		void CreateField();
		void OutTest();
	private:
		shared_ptr<GridBlock> m_grid;
		shared_ptr<FieldManager> m_field_manager;
	};


} // namespace zaran