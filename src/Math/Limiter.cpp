#include "Limiter.h"
#include "MathBasic.h"
#include<cmath>
namespace zaran
{
	double NoLimiter(const double& x, const double& y)
	{
		return 1.0;
	}

	double MinMod(const double& x, const double& y)
	{
		return 0.5 * (Sign(1.0, x) + Sign(1.0, y)) * Min(abs(x), abs(y));
	}

	double VanLeer(const double& x, const double& y)
	{
		double eps = 1.0e-12;
		return (Sign(1.0, x) + Sign(1.0, y)) * x * y / (abs(x + y) + eps);
	}

	double MixMinModVanLeer(const double& x, const double& y)
	{
		double eps = 1.0e-15;
		double z, z2, z3, xy2, absx, absy, min_van;

		z = 0.5 * (Sign(1.0, x) + Sign(1.0, y));
		absx = abs(x);
		absy = abs(y);
		xy2 = (x + x) * y;
		z2 = xy2 / (x * x + y * y + eps);
		z3 = xy2 / (absx + absy + eps);
		return z * ((1.0 - z2) * Min(absx, absy) + z2 * z3);
	}
	double VanAlbada(const double& x, const double& y)
	{
		double eps = 1.0e-6;
		return 0.5 * (Sign(1.0, x) + Sign(1.0, y)) * ((x * y) * abs(x + y)) / (x * x + y * y + eps);
	}

	double OneOrder(const double& x, const double& y)
	{
		return 0;
	}

}