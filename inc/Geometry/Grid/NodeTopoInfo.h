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
		undefined//未定义
	};
	class NodeTopoInfo
	{
	public:
		NodeTopoInfo();
		virtual~NodeTopoInfo();
		// 设置节点坐标
		void SetCoordinate(const double& x, const double& y, const double& z = 0);
		// 设置节点类型
		void SetType(const NodeType& type);
		// 设置节点标记
		void SetTag(const int& tag);
		// 设置I方向的节点模板
		void SetNeighborTemplateI(IArray& neighborTemplateI);
		// 设置J方向的节点模板
		void SetNeighborTemplateJ(IArray& neighborTemplateJ);
		// 设置K方向的节点模板
		void SetNeighborTemplateK(IArray& neighborTemplateK);
		// 设置邻居点云
		void SetNeighborCloud(const IArray& neighborCloud);
	public:
		// 返回节点类型
		const NodeType& GetType()const;
		// 返回节点标记
		const int& GetTag()const;
		// 返回I方向的节点模板
		IArray& GetNeighborTemplateI();
		// 返回J方向的节点模板
		IArray& GetNeighborTemplateJ();
		// 返回K方向的节点模板
		IArray& GetNeighborTemplateK();
		// 返回节点坐标
		DVector3D& GetCoordinate();
		// 返回邻居点云
		IArray& GetNeighborCloud();
	private:
		//节点坐标
		DVector3D coordinate_;
		//节点类型
		NodeType type_;
		// 节点标记，决定是否更新
		//=1 更新, =0 不更新
		int tag_;
		//邻居节点模板，用于计算,分别为i,j,k三个方向的模板
		//如一阶迎风记录如下
		//i-1,i,i+1
		//其中i,j,k均为自身
		IArray neighborTemplateI_;
		IArray neighborTemplateJ_;
		IArray neighborTemplateK_;
		//邻居节点点云
		//对于常规节点，为其网格线连接的邻居
		//用于计算梯度等
		IArray neighborCloud_;
		//对应的面元编号
		IArray face_;
		//对应的单元编号
		IArray cell_;
	};
}