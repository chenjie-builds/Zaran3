//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	FieldSolver.h														||
//*	@brief	流场求解器															||
//*	@author	Chen Jie.															||
//==============================================================================||

#pragma once
#include"Solver.h"
#include"FieldData.h"
namespace zaran
{
	enum class FieldSolverType
	{
		NS_Struct,
		NS_FNFDM,
	};
	/// @brief 场求解器基类
	class FieldSolver :public Solver
	{
	public:
		FieldSolver(int index, string name, SolverPara* para, GridBase* grid, FieldData* fieldData) ;
		virtual ~FieldSolver();
	public:
		// 初始化流场数据
		virtual void InitField() = 0;
		// 边界处理
		virtual void BoundaryCondition() = 0;
		// 更新场数据, 根据迭代中间结果更新场数据
		virtual void UpdateField() = 0;
		// 备份场数据
		virtual void BackupField(std::string&back_folder) = 0;
		// 获取求解器求解方程的个数
		int GetEquNum() { return m_equ_num; }
		void SetFieldData(FieldData* fieldData) ;
		// 获取Field Data
		FieldData* GetFieldData() { return m_field_data; }
	protected:
		void SetEquNum(int n) { m_equ_num = n; }
	private:
		// 求解器求解方程的个数
		int m_equ_num;
		// 场数据
		FieldData* m_field_data;
	};
}
