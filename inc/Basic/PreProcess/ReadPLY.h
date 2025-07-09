/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file ReadPLY.h
 * \brief PLYReader class, used to read PLY Point Cloud file.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include <vtkSTLReader.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include"BasicType.h"
namespace zaran
{

	class PlyReader
	{
	public:
		PlyReader() {};
		~PlyReader() {};
		void ReadPlyFile(const char* filename);
		dynamic_array<dynamic_array<double>>& GetPoints()  { return m_points; } // 获取点云数据
	private:
		dynamic_array<dynamic_array<double>> m_points; // 存储点云数据
	};
}