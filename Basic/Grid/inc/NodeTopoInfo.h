//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	NodeTopoInfo.h														||
//*	@brief	节点的拓扑结构														||
//*			包含节点类型，更新标记，邻居模板，坐标变换系数，面元单元连接信息			||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "CoordTrans.h"
#include"BasicType.h"
namespace zaran
{
	enum class NodeType//节点类型
	{
		inner,//内部点，不是边界
		inlet,//入口
		outlet,//出口
		slipWall,//滑移壁面
		noSlipWall,//无滑移壁面
		hole,//洞边界
		symmetry,//对称面
		farfield,//远场
		corner,//角点
		ghost,//ghonst point
		userDefined,//用户自定义
		undefined//未定义
	};
	class NodeTopo
	{
	public:
		NodeTopo();
		virtual~NodeTopo();
		// 设置节点坐标
		void SetCoordinate(Array<DVector3D>& coordinate) { m_coordinate = coordinate; }
		// 设置节点类型
		void SetType(Array<NodeType>& type) { m_type = type; }
		// 设置I方向的节点模板
		void SetTemplateI(Array<IArray>& template_i) { m_template_i = template_i; }
		// 设置J方向的节点模板
		void SetTemplateJ(Array<IArray>& template_j) { m_template_j = template_j; }
		// 设置K方向的节点模板
		void SetTemplateK(Array<IArray>& template_k) { m_template_k = template_k; }
		// 设置邻居点云
		void SetNeighborCloud(Array<IArray>& neighborNodeIndex) { m_neighbor_node_index = neighborNodeIndex; }
	public:
		// 返回节点类型
		Array<NodeType>& GetType() { return m_type; }
		NodeType& GetType(const int& index) { return m_type[index]; }
		// 返回I方向的节点模板
		Array<IArray>& GetTemplateI() { return m_template_i; }
		IArray& GetTemplateI(const int& index) { return m_template_i[index]; }
		// 返回J方向的节点模板
		Array<IArray>& GetTemplateJ() { return m_template_j; }
		IArray& GetTemplateJ(const int& index) { return m_template_j[index]; }
		// 返回K方向的节点模板
		Array<IArray>& GetTemplateK() { return m_template_k; }
		IArray& GetTemplateK(const int& index) { return m_template_k[index]; }
		// 返回节点坐标
		Array<DVector3D>& GetCoordinate() { return m_coordinate; }
		DVector3D& GetCoordinate(const int& index) { return m_coordinate[index]; }
		// 返回邻居点云
		Array<IArray>& GetNeighborCloud() { return m_neighbor_node_index; }
		IArray& GetNeighborCloud(const int& index) { return m_neighbor_node_index[index]; }
	private:
		//节点坐标
		Array<DVector3D> m_coordinate;
		//节点类型
		Array<NodeType> m_type;
		//邻居节点模板，用于计算,分别为i,j,k三个方向的模板
		//如一阶迎风记录如下
		//i-1,i,i+1
		//其中i,j,k均为自身
		Array<IArray> m_template_i;
		Array<IArray> m_template_j;
		Array<IArray> m_template_k;
		//邻居节点点云
		//对于常规节点，为其网格线连接的邻居
		//用于计算梯度等
		Array<IArray> m_neighbor_node_index;
		//对应的面元编号
		Array<IArray> m_face_index;
		//对应的单元编号
		Array<IArray> m_cell_index;
	};
}