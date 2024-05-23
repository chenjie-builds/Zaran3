#include<iostream>
#include "CoordTrans.h"
#include "log.h"
using namespace zaran;
void CoordTrans::CoordTransNoTime2D(const double* xRight, const double* xLeft, const double* yRight, const double* yLeft)
{
	//xξ yξ zξ tξ
	x_[0] = 0.5 * (xRight[0] - xLeft[0]);
	y_[0] = 0.5 * (xRight[1] - xLeft[1]);
	z_[0] = 0;
	t_[0] = 0;

	//xη yη zη tη,
	x_[1] = 0.5 * (yRight[0] - yLeft[0]);
	y_[1] = 0.5 * (yRight[1] - yLeft[1]);
	z_[1] = 0;
	t_[1] = 0;
	//    |x_xi   x_eta|
	// J= |            |
	//    |y_xi   y_eta|
	jacob_ = x_[0] * y_[1] - x_[1] * y_[0];
	jacob_ = 1.0 / jacob_;

	//ξx ηx ζx τx,
	xi_[0] = jacob_ * y_[1];
	eta_[0] = -jacob_ * y_[0];
	zeta_[0] = 0;
	tau_[0] = 0;

	//ξy ηy ζy τy,
	xi_[1] = -jacob_ * x_[1];
	eta_[1] = jacob_ * x_[0];
	zeta_[1] = 0;
	tau_[1] = 0;

	//ξz ηz ζz τz,
	xi_[2] = 0;
	eta_[2] = 0;
	zeta_[2] = 0;
	tau_[2] = 0;

	//xτ yτ zτ tτ,
	x_[3] = 0;
	y_[3] = 0;
	z_[3] = 0;
	t_[3] = 1;
	//ξt ηt ζt τt,
	xi_[3] = 0;
	eta_[3] = 0;
	zeta_[3] = 0;
	tau_[3] = 1;

}
void CoordTrans::CoordTransNoTime3D(const double* xRight, const double* xLeft, const double* yRight, const double* yLeft, const double* zRight, const double* zLeft)
{
	//xξ yξ zξ tξ
	x_[0] = 0.5 * (xRight[0] - xLeft[0]);
	y_[0] = 0.5 * (xRight[1] - xLeft[1]);
	z_[0] = 0.5 * (xRight[2] - xLeft[2]);
	t_[0] = 0;

	//xη yη zη tη,
	x_[1] = 0.5 * (yRight[0] - yLeft[0]);
	y_[1] = 0.5 * (yRight[1] - yLeft[1]);
	z_[1] = 0.5 * (yRight[2] - yLeft[2]);
	t_[1] = 0;

	//xζ yζ zζ tζ,
	x_[2] = 0.5 * (zRight[0] - zLeft[0]);
	y_[2] = 0.5 * (zRight[1] - zLeft[1]);
	z_[2] = 0.5 * (zRight[2] - zLeft[2]);
	t_[2] = 0;

	//J
	jacob_ = x_[0] * (y_[1] * z_[2] - y_[2] * z_[1])
		- x_[1] * (y_[0] * z_[2] - y_[2] * z_[0])
		+ x_[2] * (y_[0] * z_[1] - y_[1] * z_[0]);
	jacob_ = 1.0 / jacob_;

	//ξx ηx ζx τx,
	xi_[0] = jacob_ * (y_[1] * z_[2] - y_[2] * z_[1]);
	xi_[1] = jacob_ * (x_[2] * z_[1] - x_[1] * z_[2]);
	xi_[2] = jacob_ * (x_[1] * y_[2] - x_[2] * y_[1]);

	eta_[0] = jacob_ * (y_[2] * z_[0] - y_[0] * z_[2]);
	eta_[1] = jacob_ * (x_[0] * z_[2] - x_[2] * z_[0]);
	eta_[2] = jacob_ * (x_[2] * y_[0] - x_[0] * y_[2]);
	tau_[0] = 0;

	//ξy ηy ζy τy,
	tau_[1] = 0;

	//ξz ηz ζz τz,
	zeta_[0] = jacob_ * (y_[0] * z_[1] - y_[1] * z_[0]);
	zeta_[1] = jacob_ * (x_[1] * z_[0] - x_[0] * z_[1]);
	zeta_[2] = jacob_ * (x_[0] * y_[1] - x_[1] * y_[0]);
	tau_[2] = 0;

	//xτ yτ zτ tτ,
	x_[3] = 0;
	y_[3] = 0;
	z_[3] = 0;
	t_[3] = 1;
	//ξt ηt ζt τt,
	xi_[3] = -(x_[3] * xi_[0] + y_[3] * xi_[1] + z_[3] * xi_[2]);
	eta_[3] = -(x_[3] * eta_[0] + y_[3] * eta_[1] + z_[3] * eta_[2]);
	zeta_[3] = -(x_[3] * zeta_[0] + y_[3] * zeta_[1] + z_[3] * zeta_[2]);
	tau_[3] = 1;

}
void CoordTrans::CoordTransTime2D(double dt, const double* tRight, const double* tLeft, const double* xRight, const double* xLeft, const double* yRight, const double* yLeft)
{
	CoordTransNoTime2D(xRight, xLeft, yRight, yLeft);
	//xτ yτ zτ tτ,
	x_[3] = (tRight[0] - tLeft[0]) / dt;
	y_[3] = (tRight[1] - tLeft[1]) / dt;
	z_[3] = 0;
	t_[3] = 1;
	//ξt ηt ζt τt,
	xi_[3] = -(x_[3] * xi_[0] + y_[3] * xi_[1] + z_[3] * xi_[2]);
	eta_[3] = -(x_[3] * eta_[0] + y_[3] * eta_[1] + z_[3] * eta_[2]);
	zeta_[3] = -(x_[3] * zeta_[0] + y_[3] * zeta_[1] + z_[3] * zeta_[2]);
	tau_[3] = 1;

}
void CoordTrans::CoordTransTime3D(double dt, const double* tRight, const double* tLeft, const double* xRight, const double* xLeft, const double* yRight, const double* yLeft, const double* zRight, const double* zLeft)
{
	CoordTransNoTime3D(xRight, xLeft, yRight, yLeft, zRight, zLeft);
	//xτ yτ zτ tτ,
	x_[3] = (tRight[0] - tLeft[0]) / dt;
	y_[3] = (tRight[1] - tLeft[1]) / dt;
	z_[3] = (tRight[2] - tLeft[2]) / dt;
	t_[3] = 1;
	//ξt ηt ζt τt,
	xi_[3] = -(x_[3] * xi_[0] + y_[3] * xi_[1] + z_[3] * xi_[2]);
	eta_[3] = -(x_[3] * eta_[0] + y_[3] * eta_[1] + z_[3] * eta_[2]);
	zeta_[3] = -(x_[3] * zeta_[0] + y_[3] * zeta_[1] + z_[3] * zeta_[2]);
	tau_[3] = 1;

}


void CoordTrans::CalcCoordTrans(int dim, double dt, const double* tRight, const double* tLeft, const double* xRight, const double* xLeft, const double* yRight, const double* yLeft)
{
	CoordTrans();
	dimension_ = dim;
	if (dim == 2)
		CoordTransTime2D(dimension_, tRight, tLeft, xRight, xLeft, yRight, yLeft);
	else
		exit(0);
}

void CoordTrans::CalcCoordTrans(int dim, double dt, const double* tRight, const double* tLeft, const double* xRight, const double* xLeft, const double* yRight, const double* yLeft, const double* zRight, const double* zLeft)
{  //未加动网格 2020.07.31
	CoordTrans();
	dimension_ = dim;
	if (dim == 3)
		CoordTransTime3D(dt, tRight, tLeft, xRight, xLeft, yRight, yLeft, zRight, zLeft);
	else
		exit(0);
}

void CoordTrans::CalcCoordTrans(int dim, const double* xRight, const double* xLeft, const double* yRight, const double* yLeft)
{
	x_.resize(4);
	y_.resize(4);
	z_.resize(4);
	t_.resize(4);
	xi_.resize(4);
	eta_.resize(4);
	zeta_.resize(4);
	tau_.resize(4);
	jacob_ = 0;
	dimension_ = 0;
	dimension_ = dim;
	if (dim == 2)
		CoordTransNoTime2D(xRight, xLeft, yRight, yLeft);
	else
		exit(0);

}

void CoordTrans::CalcCoordTrans(int dim, const double* xRight, const double* xLeft, const double* yRight, const double* yLeft, const double* zRight, const double* zLeft)
{
	x_.resize(4);
	y_.resize(4);
	z_.resize(4);
	t_.resize(4);
	xi_.resize(4);
	eta_.resize(4);
	zeta_.resize(4);
	tau_.resize(4);
	jacob_ = 0;
	dimension_ = 0;
	dimension_ = dim;
	if (dim == 3)
		CoordTransNoTime3D(xRight, xLeft, yRight, yLeft, zRight, zLeft);
	else
		exit(0);

}

void CoordTrans::CalcCoordTrans(int dim, const std::vector<const double*>& coord_neib)
{
	CoordTrans();
	dimension_ = dim;
	if (dimension_ == 2 && coord_neib.size() == 4)//二维
		CoordTransNoTime2D(coord_neib[0], coord_neib[2], coord_neib[1], coord_neib[3]);
	else if (dimension_ == 3 && coord_neib.size() == 6)//三维
	{
		CoordTransNoTime3D(coord_neib[0], coord_neib[2], coord_neib[1], coord_neib[3], coord_neib[4], coord_neib[5]);
	}
	else//意外情况
	{
		Log::error("wrong neibor num ={}, dim={}", coord_neib.size(), dim);
		Log::error("EXIT NOW!");
		exit(0);
	}
}

CoordTrans::CoordTrans()
{
	x_.resize(4);
	y_.resize(4);
	z_.resize(4);
	t_.resize(4);
	xi_.resize(4);
	eta_.resize(4);
	zeta_.resize(4);
	tau_.resize(4);
	jacob_ = 0;
	dimension_ = 0;
}

CoordTrans::CoordTrans(const CoordTrans& ct)
{
	x_ = ct.xi_;
	y_ = ct.y_;
	z_ = ct.z_;
	t_ = ct.t_;
	xi_ = ct.xi_;
	eta_ = ct.eta_;
	zeta_ = ct.zeta_;
	tau_ = ct.tau_;
	jacob_ = ct.jacob_;
	dimension_ = ct.dimension_;
}
