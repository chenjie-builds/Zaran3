#include "MathBasic.h"
#include "CommonPara.h"
#include "Eigen/Dense"
#include <cmath>
namespace zaran
{

double DistanceOfTwoPoints(const double *pt1, const double *pt2)
{
    return sqrt((pt1[0] - pt2[0]) * (pt1[0] - pt2[0]) + (pt1[1] - pt2[1]) * (pt1[1] - pt2[1]) +
                (pt1[2] - pt2[2]) * (pt1[2] - pt2[2]));
}

template <typename T> bool IsErrorData(T data)
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
    double area = TriangleArea(pt1, pt2, center) + TriangleArea(pt2, pt3, center) + TriangleArea(pt3, pt4, center) +
                  TriangleArea(pt4, pt1, center);
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
        // 1. 计算点的均值
        double mean_x = 0.0, mean_y = 0.0;
        for (int i = 0; i < point_num; i++)
        {
            mean_x += coord[i][0];
            mean_y += coord[i][1];
        }
        mean_x /= point_num;
        mean_y /= point_num;

        // 2. 计算协方差矩阵元素
        double Sxx = 0.0, Syy = 0.0, Sxy = 0.0;
        for (int i = 0; i < point_num; i++)
        {
            double dx = coord[i][0] - mean_x;
            double dy = coord[i][1] - mean_y;
            Sxx += dx * dx;
            Syy += dy * dy;
            Sxy += dx * dy;
        }

        // 3. 解协方差矩阵的特征向量
        // |Sxx - λ, Sxy      |
        // |Sxy,     Syy - λ  |

        double trace = Sxx + Syy;
        double det = Sxx * Syy - Sxy * Sxy;
        double temp = sqrt(trace * trace - 4 * det);

        // 两个特征值
        double lambda1 = 0.5 * (trace + temp);
        double lambda2 = 0.5 * (trace - temp);

        // 取较小的特征值对应的特征向量作为法向量
        double lambda = (lambda1 < lambda2 ? lambda1 : lambda2);

        if (fabs(Sxy) > 1e-12)
        {
            a = lambda - Syy;
            b = Sxy;
        }
        else
        {
            // 对角矩阵的情况
            if (Sxx < Syy)
            {
                a = 1.0;
                b = 0.0;
            }
            else
            {
                a = 0.0;
                b = 1.0;
            }
        }

        // 归一化法向量
        double norm = sqrt(a * a + b * b);
        a /= norm;
        b /= norm;

        // 4. 计算 c，使直线通过均值点
        c = -(a * mean_x + b * mean_y);
    }
}
void PlaneFit3D(const double **coord, int point_num, double &A, double &B, double &C, double &D)
{
    if (point_num < 3)
    {
        A = B = C = D = 0;
        return;
    }
    // Create matrix for points
    Eigen::MatrixXd points(point_num, 3);
    for (int i = 0; i < point_num; ++i)
    {
        points(i, 0) = coord[i][0];
        points(i, 1) = coord[i][1];
        points(i, 2) = coord[i][2];
    }

    // Compute centroid
    Eigen::Vector3d centroid = points.colwise().mean();

    // Center points
    Eigen::MatrixXd centered = points.rowwise() - centroid.transpose();

    // Compute SVD
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(centered, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::Vector3d normal = svd.matrixV().col(2); // Normal is the last column (smallest singular value)

    // Plane equation: Ax + By + Cz + D = 0
    A = normal(0);
    B = normal(1);
    C = normal(2);
    D = -normal.dot(centroid);
}
} // namespace zaran