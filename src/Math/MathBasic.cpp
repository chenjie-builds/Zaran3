#include<cmath>
#include"MathBasic.h"
#include "CommonPara.h"
double zaran::DistanceOfTwoPoints(const double* pt1, const double* pt2)
{
	return sqrt((pt1[0] - pt2[0]) * (pt1[0] - pt2[0]) +
		(pt1[1] - pt2[1]) * (pt1[1] - pt2[1]) +
		(pt1[2] - pt2[2]) * (pt1[2] - pt2[2]));
}

template < typename T >
T Abs(const T& a)
{
	return (a < 0) ? -a : a;
}
template < typename T >
bool zaran::IsErrorData(T data)
{
#ifdef WIN32
	if (!_isnan(data) || Abs(data) < SMALL_NUMBER)
	{
		return false;
	}
#else
	if (Abs(data) < LARGE_NUMBER && (std::isnormal(data) || Abs(data) < SMALL_NUMBER))
	{
		return false;
	}
#endif
	return true;
}
double zaran::GetRand(const double& range_low, const double& range_high)
{
	return rand() / double(RAND_MAX) * (range_high - range_low) + range_low;
}
//southerland¹«Ê½
double zaran::Southerland(double T, double mu0, double T0, double Ts)
{
	return mu0 * pow(T / T0, 1.5) * (T0 + Ts) / (T + Ts);
}

void zaran::CircleFrom3Point(const double* pt1, const double* pt2, const double* pt3, double& radius, double* center)
{
	double x1 = pt1[0], x2 = pt2[0], x3 = pt3[0];
	double y1 = pt1[1], y2 = pt2[1], y3 = pt3[1];
	double a = x1 - x2;
	double b = y1 - y2;
	double c = x1 - x3;
	double d = y1 - y3;
	double e = ((x1 * x1 - x2 * x2) + (y1 * y1 - y2 * y2)) / 2.0;
	double f = ((x1 * x1 - x3 * x3) + (y1 * y1 - y3 * y3)) / 2.0;
	double det = b * c - a * d;
	if (abs(det) < SMALL_NUMBER)
	{
		radius = LARGE_NUMBER;
		center[0] = center[1] = 0;
		return;
	}
	double x0 = -(d * e - b * f) / det;
	double y0 = -(a * f - c * e) / det;
	center[0] = x0;
	center[1] = y0;

	radius = zaran::DistanceOfTwoPoints(pt1, center);
}
double zaran::AngleOfTwoArray(const double& x1, const double& y1, const double& x2, const double& y2)
{
	double dot = x2 * x1 + y2 * y1;
	double cross = x2 * y1 - y2 * x1;
	double angle = atan2(cross, dot);
	if (angle < 0)
		angle += 2 * PI;
	return angle;
}

double zaran::test1()
{
	double a = 1.0;
	return a;
}
