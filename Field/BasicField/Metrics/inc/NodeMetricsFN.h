#pragma once
#include "GridFNFDM.h"
namespace zaran
{
    /// @brief 网格度量类
    /// @details 用于计算网格度量
    class Metrics
    {
    public:
        Metrics(int node_num);
        ~Metrics();
        double* GetX(int iNode);
        double* GetY(int iNode);
        double* GetZ(int iNode);
        double* GetT(int iNode);
        //获取Xi 方向度量
        double* GetXi(int iNode);
        //获取Eta 方向度量
        double* GetEta(int iNode);
        //获取Zeta 方向度量
        double* GetZeta(int iNode);
        //获取Tau度量
        double* GetTau(int iNode);
        //获取Jacobian
        double& GetJacobian(int iNode);
    public:
        void CalcMetric(int iNode, const double* xRight, const double* xLeft, const double* yRight, const double* yLeft, const double* zRight = nullptr, const double* zLeft = nullptr);
    private:
        const int m_metric_num = 33;
        double* m_metric;
    };
}