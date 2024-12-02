#include <cmath>
#include "MathBasic.h"
#include "CommonPara.h"
#include "Eigen/Dense"
namespace zaran
{

	double DistanceOfTwoPoints(const double *pt1, const double *pt2)
	{
		return sqrt((pt1[0] - pt2[0]) * (pt1[0] - pt2[0]) +
					(pt1[1] - pt2[1]) * (pt1[1] - pt2[1]) +
					(pt1[2] - pt2[2]) * (pt1[2] - pt2[2]));
	}

	template <typename T>
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
	double GetRand(const double &range_low, const double &range_high)
	{
		return rand() / double(RAND_MAX) * (range_high - range_low) + range_low;
	}
	// southerland formula
	double Southerland(double T, double mu0, double T0, double Ts)
	{
		return mu0 * pow(T / T0, 1.5) * (T0 + Ts) / (T + Ts);
	}

	void CircleFrom3Point(const double *pt1, const double *pt2, const double *pt3, double &radius, double *center)
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
	double AngleOfTwoArray(const double &x1, const double &y1, const double &x2, const double &y2)
	{
		double dot = x2 * x1 + y2 * y1;
		double cross = x2 * y1 - y2 * x1;
		double angle = atan2(cross, dot);
		if (angle < 0)
			angle += 2 * PI;
		return angle;
	}
	double AngleOfTwoArray3D(const double *A, const double *B)
	{
		// 计算点积
		double dot = A[0] * B[0] + A[1] * B[1] + A[2] * B[2];

		// 计算向量 A 的模长
		double normA = sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2]);

		// 计算向量 B 的模长
		double normB = sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);

		// 计算夹角的余弦值
		double cosTheta = dot / (normA * normB);
		cosTheta = Min(1.0 - SMALL_NUMBER, cosTheta);
		cosTheta = Max(-1.0 + SMALL_NUMBER, cosTheta);

		// 计算并返回夹角（弧度）
		return acos(cosTheta);
	}
	double TriangleArea(const double *pt1, const double *pt2, const double *pt3)
	{
		double a = DistanceOfTwoPoints(pt1, pt2);
		double b = DistanceOfTwoPoints(pt2, pt3);
		double c = DistanceOfTwoPoints(pt3, pt1);
		double p = (a + b + c) / 2;
		return sqrt(p * (p - a) * (p - b) * (p - c));
	}
	// 计算四边形面积
	// 可能存在四个点不在一个平面上的情况
	// 先计算四边形的中心点，然后分割成四个三角形进行计算
	double QuadrangleArea(const double *pt1, const double *pt2, const double *pt3, const double *pt4)
	{
		double center[3] = {0};
		for (int i = 0; i < 3; i++)
		{
			center[i] = (pt1[i] + pt2[i] + pt3[i] + pt4[i]) / 4;
		}
		double area = TriangleArea(pt1, pt2, center) + TriangleArea(pt2, pt3, center) +
					  TriangleArea(pt3, pt4, center) + TriangleArea(pt4, pt1, center);
		return area;
	}

	void CrossProduct(const double *A, const double *B, double *C)
	{
		C[0] = A[1] * B[2] - A[2] * B[1];
		C[1] = A[2] * B[0] - A[0] * B[2];
		C[2] = A[0] * B[1] - A[1] * B[0];
	}
	double DotProduct(const double *A, const double *B)
	{
		return A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
	}
	void LineFit2D(const double **coord, int point_num, double &a, double &b, double &c)
	{
		if (point_num < 2)
		{
			a = b = c = 0;
			return;
		}
		bool is_same_x = true;
		for (int i = 1; i < point_num; i++)
		{
			if (abs(coord[i][0] - coord[0][0]) > EPSILON_NUMBER)
			{
				is_same_x = false;
				break;
			}
		}
		if (is_same_x)
		{
			a = 1;
			b = 0;
			c = -coord[0][0];
			return;
		}
		bool is_same_y = true;
		for (int i = 1; i < point_num; i++)
		{
			if (abs(coord[i][1] - coord[0][1]) > EPSILON_NUMBER)
			{
				is_same_y = false;
				break;
			}
		}
		if (is_same_y)
		{
			a = 0;
			b = 1;
			c = -coord[0][1];
			return;
		}
		else
		{
			double sum_x = 0.0, sum_y = 0.0;
			double sum_xy = 0.0, sum_xx = 0.0;
			for (int i = 0; i < point_num; ++i)
			{
				double x = coord[i][0];
				double y = coord[i][1];
				sum_x += x;
				sum_y += y;
				sum_xy += x * y;
				sum_xx += x * x;
			}
			double denominator = point_num * sum_xx - sum_x * sum_x;
			if (abs(denominator) < SMALL_NUMBER)
			{
				a = 1;
				b = 0;
				c = -sum_x / point_num;
				return;
			}
			a = (point_num * sum_xy - sum_x * sum_y) / denominator;
			b = (sum_xx * sum_y - sum_x * sum_xy) / denominator;
			c = b;
			b = -1;
		}
	}
	void PlaneFit3D(const double **coord, int point_num, double &A, double &B, double &C, double &D)
	{
		if (point_num < 3)
		{
			A = B = C = D = 0;
			return;
		}
		//检查是否所有点x坐标都相同
		bool is_same_x = true;
		for (int i = 1; i < point_num; i++)
		{
			if (abs(coord[i][0] - coord[0][0]) > EPSILON_NUMBER)
			{
				is_same_x = false;
				break;
			}
		}
		//检查是否所有点y坐标都相同
		bool is_same_y = true;
		for (int i = 1; i < point_num; i++)
		{
			if (abs(coord[i][1] - coord[0][1]) > EPSILON_NUMBER)
			{
				is_same_y = false;
				break;
			}
		}
		//检查是否所有点z坐标都相同
		bool is_same_z = true;
		for (int i = 1; i < point_num; i++)
		{
			if (abs(coord[i][2] - coord[0][2]) > EPSILON_NUMBER)
			{
				is_same_z = false;
				break;
			}
		}
		if(is_same_x)
		{
			A = 1;
			B = C = 0;
			D = -coord[0][0]; 
			return;
		}
		if(is_same_y)
		{
			A = B = 0;
			C = 1;
			D = -coord[0][1]; // 使用y的平均值
			return;
		}
		if(is_same_z)
		{
			A = B = 0;
			C = 1;
			D = -coord[0][2]; // 使用z的平均值
			return;
		}
		double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
		double sum_xy = 0.0, sum_xz = 0.0, sum_yz = 0.0;
		double sum_xx = 0.0, sum_yy = 0.0, sum_zz = 0.0;

		// 计算必要的求和
		for (int i = 0; i < point_num; ++i)
		{
			double x = coord[i][0];
			double y = coord[i][1];
			double z = coord[i][2];
			sum_x += x;
			sum_y += y;
			sum_z += z;
			sum_xy += x * y;
			sum_xz += x * z;
			sum_yz += y * z;
			sum_xx += x * x;
			sum_yy += y * y;
			sum_zz += z * z;
		}

		// 构建矩阵 A
		Eigen::Matrix3d A_matrix;
		A_matrix << sum_xx, sum_xy, sum_xz,
					sum_xy, sum_yy, sum_yz,
					sum_xz, sum_yz, sum_zz;

		Eigen::Vector3d B_vector(sum_x, sum_y, sum_z);

		Eigen::Vector3d plane_params = A_matrix.colPivHouseholderQr().solve(B_vector);

		A = plane_params[0];
		B = plane_params[1];
		C = plane_params[2];

		// 计算常数项 D
		D = -1;
}
} // namespace zaran