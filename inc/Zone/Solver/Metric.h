/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file Metric.h
 * \brief Metric class, used to calculate the metric of the grid.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "GridFNFDM.h"
namespace zaran
{
	/// @brief 网格度量类
	/// @details 用于计算网格度量
	class Metric
	{
	public:
		Metric(count_type node_num);
		~Metric();
		/// @brief 获取逆度量系数
		/// @param dim 维度索引, 0:x, 1:y, 2:z, 3:t
		/// @param iNode 节点索引
		/// @return 逆度量系数，数组大小为4，分别为xi,eta,zeta,tau
		double* GetInvMetrics(int dim, int iNode);
		/// @brief 获取度量系数
		/// @param dim 维度索引, 0:xi, 1:eta, 2:zeta, 3:tau
		/// @param iNode 节点索引
		/// @return 度量系数，数组大小为4，分别为x,y,z,t
		double* GetMetrics(int dim, int iNode);
		double* GetX(int iNode);
		double* GetY(int iNode);
		double* GetZ(int iNode);
		double* GetT(int iNode);
		//获取Xi 方向度量
		double* GetXi(int iNode);
		//获取Eta 方向度量
		double* GetEta(int iNode);
		//获取Zeta 方向度量
		double* GetZeta(int iNode);
		//获取Tau度量
		double* GetTau(int iNode);
		//获取Jacobian
		double& GetJacobian(int iNode);
	public:
		void CalcMetric(int iNode, const double* xRight, const double* xLeft, const double* yRight, const double* yLeft, const double* zRight = nullptr, const double* zLeft = nullptr);
	private:
		const count_type m_metric_num = 33;
		dynamic_array<double> m_metric;
	};
}