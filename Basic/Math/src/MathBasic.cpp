#include<cmath>
#include"MathBasic.h"
#include "CommonPara.h"
namespace zaran
{

	double DistanceOfTwoPoints(const double* pt1, const double* pt2)
	{
		return sqrt((pt1[0] - pt2[0]) * (pt1[0] - pt2[0]) +
			(pt1[1] - pt2[1]) * (pt1[1] - pt2[1]) +
			(pt1[2] - pt2[2]) * (pt1[2] - pt2[2]));
	}

	template < typename T >
	bool IsErrorData(T data)
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
	double GetRand(const double& range_low, const double& range_high)
	{
		return rand() / double(RAND_MAX) * (range_high - range_low) + range_low;
	}
	//southerland formula
	double Southerland(double T, double mu0, double T0, double Ts)
	{
		return mu0 * pow(T / T0, 1.5) * (T0 + Ts) / (T + Ts);
	}

	void CircleFrom3Point(const double* pt1, const double* pt2, const double* pt3, double& radius, double* center)
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
	double AngleOfTwoArray(const double& x1, const double& y1, const double& x2, const double& y2)
	{
		double dot = x2 * x1 + y2 * y1;
		double cross = x2 * y1 - y2 * x1;
		double angle = atan2(cross, dot);
		if (angle < 0)
			angle += 2 * PI;
		return angle;
	}
	double AngleOfTwoArray3D(const double* A, const double* B)
	{
		// 计算点积
		double dot = A[0] * B[0] + A[1] * B[1] + A[2] * B[2];

		// 计算向量 A 的模长
		double normA = sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2]);

		// 计算向量 B 的模长
		double normB = sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);

		// 计算夹角的余弦值
		double cosTheta = dot / (normA * normB);

		// 计算并返回夹角（弧度）
		return acos(cosTheta);
	}
	double TriangleArea(const double* pt1, const double* pt2, const double* pt3)
	{
		double a = DistanceOfTwoPoints(pt1, pt2);
		double b = DistanceOfTwoPoints(pt2, pt3);
		double c = DistanceOfTwoPoints(pt3, pt1);
		double p = (a + b + c) / 2;
		return sqrt(p * (p - a) * (p - b) * (p - c));
	}
	//计算四边形面积
	//可能存在四个点不在一个平面上的情况
	//先计算四边形的中心点，然后分割成四个三角形进行计算
	double QuadrangleArea(const double* pt1, const double* pt2, const double* pt3, const double* pt4)
	{
		double center[3] = { 0 };
		for (int i = 0; i < 3; i++)
		{
			center[i] = (pt1[i] + pt2[i] + pt3[i] + pt4[i]) / 4;
		}
		double area = TriangleArea(pt1, pt2, center) + TriangleArea(pt2, pt3, center) +
			TriangleArea(pt3, pt4, center) + TriangleArea(pt4, pt1, center);
		return area;
	}

	void CrossProduct(const double* A, const double* B, double* C)
	{
		C[0] = A[1] * B[2] - A[2] * B[1];
		C[1] = A[2] * B[0] - A[0] * B[2];
		C[2] = A[0] * B[1] - A[1] * B[0];
	}
    double DotProduct(const double *A, const double *B)
    {
        return A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
    }
}