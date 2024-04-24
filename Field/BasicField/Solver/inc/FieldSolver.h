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
		NS_2D,
		NS_2D_Struct,
		NS_3D,
		NS_3D_Struct,
		NS_ZaRan_2D,
		NS_ZaRan_3D
	};
	/// @brief 场求解器基类
	class FieldSolver :public Solver
	{
	public:
		FieldSolver() {}
		virtual ~FieldSolver() {}
	public:
		// 初始化流场数据
		virtual void InitField() = 0;
		// 边界处理
		virtual void BoundaryCondition() = 0;
		// 生成场数据, 开辟内存
		virtual void CreateData() = 0;
		// 注册场数据, 根据求解器类型对场数据分类注册
		virtual void RegisterFieldData() = 0;
		// 更新场数据, 根据迭代中间结果更新场数据
		virtual void UpdateField() = 0;
		// 备份场数据
		virtual void BackupField(std::string&back_folder) = 0;

		// 获取求解器求解方程的个数
		int GetNumberOfEquations() { return m_NumberOfEquations; }
		void SetFieldData(Ptr<FieldData>& fieldData) { m_FieldData = fieldData; }
		// 获取Field Data
		Ptr<FieldData>& GetFieldData() { return m_FieldData; }
	protected:
		void SetNumberOfEquations(int n) { m_NumberOfEquations = n; }
	private:
		// 求解器求解方程的个数
		int m_NumberOfEquations;
		Ptr<FieldData> m_FieldData;
	};
}
