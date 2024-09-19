//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	STL.h																||
//*	@brief	STL 模型																||
//*	@author	Chen Jie.															||
//==============================================================================||

#pragma once
#include "Model.h"
#include <vtkExtractCells.h>
#include <vtkIdList.h>
#include <vtkLineSource.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkCellLocator.h>
#include <vtkOBBTree.h>
#include <vtkPoints.h>
#include <vtkSphereSource.h>
#include<vtkImplicitPolyDataDistance.h>
namespace zaran
{
	class PolyDataModel : public Model
	{
	public:
		PolyDataModel() : m_tolerance(1e-6) {}
		void SetPolyData(vtkSmartPointer<vtkPolyData> &polyData, double tolerance);
		virtual bool InModel(const double *pt) const override;										// 判断点pt是否在
		virtual void GetClosestPoint(const double *point_input, double *point_find) const override; // 求出点pt在模型离散点数组中最近的点
		virtual double GetClosestDistance(const double *pt) const override;
		// 判断是否为封闭模型
		bool IsClosed();
		// 设置容忍度
		void SetTolerance(double tol);
		vtkSmartPointer<vtkPolyData> &GetPolyData() { return m_polyData; }

	private:
		// 容忍度
		double m_tolerance;
		vtkSmartPointer<vtkImplicitPolyDataDistance> m_distance;
		vtkSmartPointer<vtkOBBTree> m_cell_locator;
		vtkSmartPointer<vtkCellLocator> m_cellLocator;
		
		// polyData
		vtkSmartPointer<vtkPolyData> m_polyData;
	};
}