/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file DataManagerNSStruct.h
 * \brief DataManagerNSStruct class, used to manage the data of the NS equation in a structured grid.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once 
#include"DataManagerNS.h"

namespace zaran
{
	class DataManagerNSStruct :public DataManagerNS
	{
	public:
		DataManagerNSStruct(shared_ptr<FieldData> fieldData, int ni, int nj, int nk);
		~DataManagerNSStruct();
	public:
		void CreateData()override;
		void RegisterData()override;
	public:
		//! @brief 设置半点左侧基本变量
		/// @param iEqu 变量索引
		/// @param iDim 方向索引
		/// @param iNode 节点索引
		/// @param value 输入的值
		void SetMidNodePrimLeft(int iEqu, int iDim, int iNode, double value);
		/// @brief  设置半点左侧基本变量
		/// @param eq_index 变量索引
		/// @param dim 半点的方向索引
		/// @param data_index 数据索引
		/// @param data_value基本变量值
		void SetMidNodePrimLeft(index_type eq_index, dimension_type dim, index_type data_index, double data_value);
		//! @brief 设置半点右侧基本变量
		/// @param iEqu 变量索引
		/// @param iDim 方向索引
		/// @param iNode 节点索引
		/// @param value 输入的值
		void SetMidNodePrimRight(int iEqu, int iDim, int iNode, double value);
		void SetMidNodePrim(int iEqu, int iDim, int iNode, double value_left, double value_right);
		//! @brief 设置半点数值通量
		/// @param iEqu 变量索引
		/// @param iDim 方向索引
		/// @param iNode 节点索引
		/// @param value 输入的值
		void SetMidNodeFlux(int iEqu, int iDim, int iNode, double value);
		//! @brief 获取半点左侧基本变量
		/// @param iEqu 变量索引
		/// @param iDim 方向索引
		/// @param iNode 节点索引
		double GetMidNodePrimLeft(int iEqu, int iDim, int iNode);
		/// @brief 获取半点右侧基本变量
		/// @param iEqu 变量索引
		/// @param iDim 方向索引
		/// @param iNode 节点索引
		double GetMidNodePrimRight(int iEqu, int iDim, int iNode);
		/// @brief 获取半点数值通量
		/// @param iEqu 变量索引
		/// @param iDim 方向索引
		/// @param iNode 节点索引
		double GetMidNodeFlux(int iEqu, int iDim, int iNode);
	private:
		/// @brief 半点左侧的基本变量，1D变量索引，2D方向索引(0:i+1/2,1:j+1/2,2:k+1/2)，3D节点索引
		double*** m_midnode_prim_left;
		/// @brief 半点右侧的基本变量，1D变量索引，2D方向索引(1:i+1/2,2:j+1/2,3:k+1/2)，3D节点索引
		double*** m_midnode_prim_right;
		/// @brief 半点的数值通量，1D变量索引，2D方向索引(1:i+1/2,2:j+1/2,3:k+1/2)，3D节点索引
		double*** m_midnode_flux;
	};
}