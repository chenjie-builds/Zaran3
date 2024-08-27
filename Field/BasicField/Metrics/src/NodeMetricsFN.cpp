#include "NodeMetricsFN.h"
#include "Log.h"
namespace zaran
{
    Metrics::Metrics(int node_num)
    {
        m_metric = new double[node_num * m_metric_num];
    }

    Metrics::~Metrics()
    {
        delete[] m_metric;
    }
    double *Metrics::GetX(int iNode)
    {
        return m_metric + iNode * m_metric_num;
    }
    double *Metrics::GetY(int iNode)
    {
        return m_metric + iNode * m_metric_num + 4;
    }
    double *Metrics::GetZ(int iNode)
    {
        return m_metric + iNode * m_metric_num + 8;
    }
    double *Metrics::GetT(int iNode)
    {
        return m_metric + iNode * m_metric_num + 12;
    }
    double *Metrics::GetXi(int iNode)
    {
        return m_metric + iNode * m_metric_num + 16;
    }

    double *Metrics::GetEta(int iNode)
    {
        return m_metric + iNode * m_metric_num + 20;
    }

    double *Metrics::GetZeta(int iNode)
    {
        return m_metric + iNode * m_metric_num + 24;
    }

    double *Metrics::GetTau(int iNode)
    {
        return m_metric + iNode * m_metric_num + 28;
    }
    double &Metrics::GetJacobian(int iNode)
    {
        return m_metric[iNode * m_metric_num + 32];
    }

    void Metrics::CalcMetric(int iNode, const double *iRight, const double *iLeft, const double *jRight, const double *jLeft, const double *kRight, const double *kLeft)
    {
        double coef_x[4], coef_y[4], coef_z[4], coef_t[4];
        double *coef_xi = GetXi(iNode);
        double *coef_eta = GetEta(iNode);
        double *coef_zeta = GetZeta(iNode);
        double *coef_tau = GetTau(iNode);
        double &jacobian = GetJacobian(iNode);
        coef_x[0] = 0.5 * (iRight[0] - iLeft[0]);
        coef_x[1] = 0.5 * (jRight[0] - jLeft[0]);
        coef_x[3] = 0.0;
        coef_y[0] = 0.5 * (iRight[1] - iLeft[1]);
        coef_y[1] = 0.5 * (jRight[1] - jLeft[1]);
        coef_y[3] = 0.0;
        coef_t[0] = 0.0;
        coef_t[1] = 0.0;
        coef_t[2] = 0.0;
        coef_t[3] = 1.0;
        if (kLeft == nullptr && kRight == nullptr)
        {
            coef_x[2] = 0.0;
            coef_y[2] = 0.0;
            coef_z[0] = 0.0;
            coef_z[1] = 0.0;
            coef_z[2] = 1.0;
            coef_z[3] = 0.0;
        }
        else
        {
            coef_x[2] = 0.5 * (kRight[0] - kLeft[0]);
            coef_y[2] = 0.5 * (kRight[1] - kLeft[1]);
            coef_z[0] = 0.5 * (iRight[2] - iLeft[2]);
            coef_z[1] = 0.5 * (jRight[2] - jLeft[2]);
            coef_z[2] = 0.5 * (kRight[2] - kLeft[2]);
            coef_z[3] = 0.0;
        }
        jacobian = coef_x[0] * (coef_y[1] * coef_z[2] - coef_y[2] * coef_z[1]) - coef_x[1] * (coef_y[0] * coef_z[2] - coef_y[2] * coef_z[0]) + coef_x[2] * (coef_y[0] * coef_z[1] - coef_y[1] * coef_z[0]);
        jacobian = 1.0 / jacobian;

        coef_xi[0] = (coef_y[1] * coef_z[2] - coef_y[2] * coef_z[1]);
        coef_xi[1] = -(coef_x[1] * coef_z[2] - coef_x[2] * coef_z[1]);
        coef_xi[2] = (coef_x[1] * coef_y[2] - coef_x[2] * coef_y[1]);
        coef_xi[3] = -(coef_x[3] * coef_xi[0] + coef_x[3] * coef_xi[1] + coef_x[3] * coef_xi[2]);

        coef_eta[0] = -(coef_y[0] * coef_z[2] - coef_y[2] * coef_z[0]);
        coef_eta[1] = (coef_x[0] * coef_z[2] - coef_x[2] * coef_z[0]);
        coef_eta[2] = -(coef_x[0] * coef_y[2] - coef_x[2] * coef_y[0]);
        coef_eta[3] = -(coef_y[3] * coef_eta[0] + coef_y[3] * coef_eta[1] + coef_y[3] * coef_eta[2]);

        coef_zeta[0] = (coef_y[0] * coef_z[1] - coef_y[1] * coef_z[0]);
        coef_zeta[1] = -(coef_x[0] * coef_z[1] - coef_x[1] * coef_z[0]);
        coef_zeta[2] = (coef_x[0] * coef_y[1] - coef_x[1] * coef_y[0]);
        coef_zeta[3] = -(coef_z[3] * coef_zeta[0] + coef_z[3] * coef_zeta[1] + coef_z[3] * coef_zeta[2]);

        coef_tau[0] = 0.0;
        coef_tau[1] = 0.0;
        coef_tau[2] = 0.0;
        coef_tau[3] = 1.0;
    }

}