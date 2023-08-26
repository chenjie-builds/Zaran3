#include "SpecialField.h"

void CalcSupersonicVortex(const double& x, const double& y, const double& x0, const double y0, const double& ci, const double& Mi, const double& rhoi, Eigen::VectorXd& primitive)
{
	if (primitive.size() != 5)
	{
		primitive.resize(5);
	}
	double gama = 1.4;
	double ri = 1.0;
	double dx = x - x0;
	double dy = y - y0;
	double r = sqrt(dx * dx + dy * dy);
	double U = ci * Mi / r;
	double u = y * U / r;
	double v = -x * U / r;
	double rho = rhoi * pow(1 + (gama - 1) / 2.0 * Mi * Mi * (1 - pow(ri, 2.0) / pow(r, 2.0)), 1.0 / (gama - 1));
	double p = pow(rho, gama) / gama;
	primitive[0] = rho;
	primitive[1] = u;
	primitive[2] = v;
	primitive[3] = 0;
	primitive[4] = p;
}

void CalcIsentropicVortex(const double& x, const double& y, const double& beta, Eigen::VectorXd& primitive)
{
	double xv, yv;
	xv = yv = 0;
	double x1 = x - xv;
	double y1 = y - yv;
	double r = sqrt(x1 * x1 + y1 * y1);
	if (primitive.size() != 5)
	{
		primitive.resize(5);
	}
	double gamma = 1.4;
	double pi = 3.141592653;
	double density, u, v, w, p;
	density = pow(1 - (gamma - 1) * beta * beta * exp(1 - r * r) / (8 * gamma * pi * pi), 1.0 / (gamma - 1));
	u = beta * exp(0.5 * (1 - r * r)) * (-y1) / (2 * pi);
	v = beta * exp(0.5 * (1 - r * r)) * (x1) / (2 * pi);
	w = 0;
	p = pow(density, gamma);

	//density = r * r;

	primitive[0] = density;
	primitive[1] = u;
	primitive[2] = v;
	primitive[3] = w;
	primitive[4] = p;
	//ZaranLog::info("rho,u,v,p={},{},{},{}", density, u, v, p);
}