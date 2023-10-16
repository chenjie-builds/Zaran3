//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	CoordTrans.h														||
//*	@brief	坐标变换系数类, 计算节点的坐标变换系数									||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include"BasicType.h"
//坐标变换
//对坐标进行标记的right 和left表示在坐标轴上的方向，如：x轴向，位于坐标轴右侧的为xright，左侧的为xleft

namespace zaran
{
	class CoordTrans
	{
	public:
		using Coordinate = DVector3D;
		CoordTrans();
		CoordTrans(const CoordTrans& ct);
		//直接输入坐标变换所有的点，注意点的数量和维数对应
		void CalcCoordTrans(int dim, double dt, Coordinate tRight, Coordinate tLeft, Coordinate xRight, Coordinate xLeft, Coordinate yRight, Coordinate yLeft);
		void CalcCoordTrans(int dim, double dt, Coordinate tRight, Coordinate tLeft, Coordinate xRight, Coordinate xLeft, Coordinate yRight, Coordinate yLeft, Coordinate zRight, Coordinate zLeft);
		void CalcCoordTrans(int dim, Coordinate xRight, Coordinate xLeft, Coordinate yRight, Coordinate yLeft);
		void CalcCoordTrans(int dim, Coordinate xRight, Coordinate xLeft, Coordinate yRight, Coordinate yLeft, Coordinate zRight, Coordinate zLeft);
		//输入坐标变换所需要点的vector进行坐标变换初始化,注意vector大小和维数对应
		void CalcCoordTrans(int dim, const std::vector<Coordinate>& coord_neib);
		//返回坐标变换系数
		const DArray& GetX() const { return x_; }
		const DArray& GetY()const { return y_; }
		const DArray& GetZ()const { return z_; }
		const DArray& GetT()const { return t_; }
		const DArray& GetXi()const { return xi_; }
		const DArray& GetEta()const { return eta_; }
		const DArray& GetZeta()const { return zeta_; }
		const DArray& GetTau() const { return tau_; }
		const double& J()const { return jacob_; }
		const double& dim()const { return dimension_; }
	private:
		void CoordTransTime2D(double dt, Coordinate tRight, Coordinate tLeft, Coordinate xRight, Coordinate xLeft, Coordinate yRight, Coordinate yLeft);
		void CoordTransTime3D(double dt, Coordinate tRight, Coordinate tLeft, Coordinate xRight, Coordinate xLeft, Coordinate yRight, Coordinate yLeft, Coordinate zRight, Coordinate zLeft);
		void CoordTransNoTime2D(Coordinate xRight, Coordinate xLeft, Coordinate yRight, Coordinate yLeft);
		void CoordTransNoTime3D(Coordinate xRight, Coordinate xLeft, Coordinate yRight, Coordinate yLeft, Coordinate zRight, Coordinate zLeft);
	private:
		//xξ,xη,xζ,xτ
		DArray x_;
		//yξ,yη,yζ,yτ
		DArray y_;
		//zξ,zη,zζ,zτ
		DArray z_;
		//tξ,tη,tζ,tτ
		DArray t_;
		//ξx,ξy,ξz,ξt
		DArray xi_;
		//ηx,ηy,ηz,ηt
		DArray eta_;
		//ζx,ζy,ζz,ζt
		DArray zeta_;
		//τx,τy,τz,τt
		DArray tau_;
		double jacob_;
		int dimension_;//维数
	};
}