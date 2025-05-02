/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file STLReader.h
 * \brief STLReader class, used to read STL file.
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

	class STLReader
	{
	public:
		STLReader() {};
		~STLReader() {};
		void ReadSTLFile(const char* filename);
		vtkSmartPointer<vtkPolyData>&  GetPolyData() { return m_poly_data; }
	private:
		vtkSmartPointer<vtkPolyData> m_poly_data;
	};
}