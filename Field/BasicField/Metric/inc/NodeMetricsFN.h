#pragma once
#include "GridFNFDM.h"
namespace zaran
{
    class NodeMetric
    {
    public:
        NodeMetric(int node_num);
        ~NodeMetric();
        //计算度量
        // void CalcMetric();
        //检查Jacobian是否存在错误
        // void CheckJacobian();
        //对存在错误的Jacobian进行修正
        // void FixJacobian();
        //获取Xi 方向度量
        const  double* GetMetricXi(int iNode)const;
        //获取Eta 方向度量
        const double* GetMetricEta(int iNode)const;
        //获取Zeta 方向度量
        const double* GetMetricZeta(int iNode)const;
        //获取Tau度量
        const double* GetMetricTau(int iNode)const;
        //获取Jacobian
        double GetJacobian(int iNode)const;
    public :
        void CalcMetric(int iNode, const double* xRight, const double* xLeft, const double* yRight, const double* yLeft, const double* zRight = nullptr, const double* zLeft = nullptr);
    private:
        const int m_metric_num = 17;
        double* m_metric;
    };
}