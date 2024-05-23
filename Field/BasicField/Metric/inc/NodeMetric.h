#pragma once
#include "GridFNFDM.h"
namespace zaran
{
    class NodeMetric
    {
    public:
        NodeMetric(GridFN* grid);
        ~NodeMetric();
        //计算度量
        void CalcMetric();
        //检查Jacobian是否存在错误
        void CheckJacobian();
        //对存在错误的Jacobian进行修正
        void FixJacobian();
        //获取Xi 方向度量
        const  double* GetMetricXi(int iNode)const;
        //获取Eta 方向度量
        const double* GetMetricEta(int iNode)const;
        //获取Zeta 方向度量
        const double* GetMetricZeta(int iNode)const;
        //获取Tau度量
        const double* GetMetricTau(int iNode)const;
        //获取Jacobian
        const double GetJacobian(int iNode)const;
        const double GetMaxJacobian()const;
        const double GetMinJacobian()const;
        const int GetMaxJacobianNode()const;
        const int GetMinJacobianNode()const;   
    private:
        void CalcMetric(double* metric, const double* xRight, const double* xLeft, const double* yRight, const double* yLeft, const double* zRight = nullptr, const double* zLeft = nullptr);
        bool ErrorJacobi(double jacobi)const;
    private:
        const int m_metric_num = 17;
        GridFN* m_grid;
        double* m_metric;
        double max_jacobian;
        double min_jacobian;
        int max_jacobian_node;
        int min_jacobian_node;
        bool m_exist_error_jacobi;
    };
}