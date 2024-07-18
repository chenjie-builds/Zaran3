#include "NodeMetricsFN.h"
#include "Log.h"
namespace zaran
{
    Metrics::Metrics(int node_num)
    {
        m_metric = new double[node_num * (4 * 4 + 1)];
    }

    Metrics::~Metrics()
    {
        delete[] m_metric;
    }

    // void NodeMetric::CalcMetric()
    // {
    //     auto node = m_grid->GetNode();
    //     int node_num = node->GetNodeNum();
    //     const double* xRight, * xLeft, * yRight, * yLeft, * zRight, * zLeft;
    //     xRight = xLeft = yRight = yLeft = zRight = zLeft = nullptr;
    //     max_jacobian = -LARGE_NUMBER;
    //     min_jacobian = LARGE_NUMBER;
    //     for (int iNode = 0;iNode < node_num;++iNode)
    //     {
    //         if (node->GetType(iNode) != NodeType::inner)
    //         {
    //             continue;
    //         }
    //         auto neighbors = node->GetNeighborNode(iNode);
    //         xLeft = node->GetCoord(neighbors[0]);
    //         xRight = node->GetCoord(neighbors[1]);
    //         yLeft = node->GetCoord(neighbors[2]);
    //         yRight = node->GetCoord(neighbors[3]);
    //         if (m_grid->GetDimension() == 3)
    //         {
    //             zLeft = node->GetCoord(neighbors[4]);
    //             zRight = node->GetCoord(neighbors[5]);
    //         }
    //         CalcMetric(m_metric + iNode * m_metric_num, xRight, xLeft, yRight, yLeft, zRight, zLeft);
    //         if (GetJacobian(iNode) > max_jacobian)
    //         {
    //             max_jacobian = GetJacobian(iNode);
    //             max_jacobian_node = iNode;
    //         }
    //         if (GetJacobian(iNode) < min_jacobian)
    //         {
    //             min_jacobian = GetJacobian(iNode);
    //             min_jacobian_node = iNode;
    //         }
    //     }
    // }

    // void NodeMetric::CheckJacobian()
    // {
    //     auto node = m_grid->GetNode();
    //     int node_num = node->GetNodeNum();
    //     for (int iNode = 0;iNode < node_num;++iNode)
    //     {
    //         if (node->GetType(iNode) != NodeType::inner)
    //         {
    //             continue;
    //         }
    //         if (ErrorJacobi(GetJacobian(iNode)))
    //         {
    //             m_exist_error_jacobi = true;
    //             return;
    //         }
    //     }
    // }
    // void NodeMetric::FixJacobian()
    // {
    //     if (!m_exist_error_jacobi)
    //         return;
    //     auto node = m_grid->GetNode();
    //     int node_num = node->GetNodeNum();
    //     const double* xRight, * xLeft, * yRight, * yLeft, * zRight, * zLeft;
    //     xRight = xLeft = yRight = yLeft = zRight = zLeft = nullptr;
    //     std::vector<int> temp_neibor;
    //     for (int iNode = 0;iNode < node_num;++iNode)
    //     {
    //         if (node->GetType(iNode) != NodeType::inner)
    //         {
    //             continue;
    //         }
    //         auto neighbors = node->GetNeighborNode(iNode);
    //         int neighbor_num = node->GetNeighborNodeNum(iNode);
    //         temp_neibor.resize(neighbor_num);
    //         for (int i = 0;i < neighbor_num;++i)
    //         {
    //             temp_neibor[i] = neighbors[i];
    //         }
    //         // 如果雅可比行列式为0或者无穷大, 可能是因为坐标轴平行
    //         if (abs(GetJacobian(iNode) < SMALL_NUMBER) || isnan(GetJacobian(iNode)) || isinf(GetJacobian(iNode)))
    //         {
    //             std::swap(temp_neibor[3], temp_neibor[4]);
    //             xLeft = node->GetCoord(neighbors[0]);
    //             xRight = node->GetCoord(neighbors[1]);
    //             yLeft = node->GetCoord(neighbors[2]);
    //             yRight = node->GetCoord(neighbors[3]);
    //             if (m_grid->GetDimension() == 3)
    //             {
    //                 zLeft = node->GetCoord(neighbors[4]);
    //                 zRight = node->GetCoord(neighbors[5]);
    //             }
    //             CalcMetric(m_metric + iNode * m_metric_num, xRight, xLeft, yRight, yLeft, zRight, zLeft);
    //             node->SetNeighborNode(iNode, neighbor_num, temp_neibor.data());
    //         }
    //         //如果雅可比行列式为负, 交换K方向的两个邻居节点
    //         if (GetJacobian(iNode) < 0)
    //         {
    //             std::swap(temp_neibor[4], temp_neibor[5]);
    //             xLeft = node->GetCoord(neighbors[0]);
    //             xRight = node->GetCoord(neighbors[1]);
    //             yLeft = node->GetCoord(neighbors[2]);
    //             yRight = node->GetCoord(neighbors[3]);
    //             if (m_grid->GetDimension() == 3)
    //             {
    //                 zLeft = node->GetCoord(neighbors[4]);
    //                 zRight = node->GetCoord(neighbors[5]);
    //             }
    //             CalcMetric(m_metric + iNode * m_metric_num, xRight, xLeft, yRight, yLeft, zRight, zLeft);
    //             node->SetNeighborNode(iNode, neighbor_num, temp_neibor.data());
    //         }
    //         // 如果修正后的雅可比行列式仍然为错误, 则报错
    //         if (ErrorJacobi(GetJacobian(iNode)))
    //         {
    //             Log::error("NodeMetric::FixJacobian: fix jacobian error failed");
    //             Log::error("NodeMetric::FixJacobian: node id = {}, Coordinate = ({},{},{}), jacobian = {}", iNode, node->GetCoord(iNode)[0], node->GetCoord(iNode)[1], node->GetCoord(iNode)[2], GetJacobian(iNode));
    //             Log::error("Exit...");
    //             exit(1);
    //         }
    //     }
    //     // 修复完成
    //     m_exist_error_jacobi = false;
    // }
    const double* Metrics::GetXi(int iNode) const
    {
        return m_metric + iNode * 17;
    }

    const double* Metrics::GetEta(int iNode) const
    {
        return m_metric + iNode * 17 + 4;
    }

    const double* Metrics::GetZeta(int iNode) const
    {
        return m_metric + iNode * 17 + 8;
    }

    const double* Metrics::GetMetricTau(int iNode) const
    {
        return m_metric + iNode * 17 + 12;
    }
    double Metrics::GetJacobian(int iNode) const
    {
        return m_metric[iNode * 17 + 16];
    }



    void Metrics::CalcMetric(int iNode, const double* iRight, const double* iLeft, const double* jRight, const double* jLeft, const double* kRight, const double* kLeft)
    {
        double coef_x[4], coef_y[4], coef_z[4], coef_t[4];
        double* coef_xi = m_metric + iNode * m_metric_num;
        double* coef_eta = coef_xi + 4;
        double* coef_zeta = coef_eta + 4;
        double* coef_tau = coef_zeta + 4;
        double& jacobian = *(coef_tau + 4);
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
        jacobian = coef_x[0] * (coef_y[1] * coef_z[2] - coef_y[2] * coef_z[1])
            - coef_x[1] * (coef_y[0] * coef_z[2] - coef_y[2] * coef_z[0])
            + coef_x[2] * (coef_y[0] * coef_z[1] - coef_y[1] * coef_z[0]);
        jacobian = 1.0 / jacobian;

        coef_xi[0] = jacobian * (coef_y[1] * coef_z[2] - coef_y[2] * coef_z[1]);
        coef_xi[1] = -jacobian * (coef_x[1] * coef_z[2] - coef_x[2] * coef_z[1]);
        coef_xi[2] = jacobian * (coef_x[1] * coef_y[2] - coef_x[2] * coef_y[1]);
        coef_xi[3] = -(coef_x[3] * coef_xi[0] + coef_x[3] * coef_xi[1] + coef_x[3] * coef_xi[2]);

        coef_eta[0] = -jacobian * (coef_y[0] * coef_z[2] - coef_y[2] * coef_z[0]);
        coef_eta[1] = jacobian * (coef_x[0] * coef_z[2] - coef_x[2] * coef_z[0]);
        coef_eta[2] = -jacobian * (coef_x[0] * coef_y[2] - coef_x[2] * coef_y[0]);
        coef_eta[3] = -(coef_y[3] * coef_eta[0] + coef_y[3] * coef_eta[1] + coef_y[3] * coef_eta[2]);

        coef_zeta[0] = jacobian * (coef_y[0] * coef_z[1] - coef_y[1] * coef_z[0]);
        coef_zeta[1] = -jacobian * (coef_x[0] * coef_z[1] - coef_x[1] * coef_z[0]);
        coef_zeta[2] = jacobian * (coef_x[0] * coef_y[1] - coef_x[1] * coef_y[0]);
        coef_zeta[3] = -(coef_z[3] * coef_zeta[0] + coef_z[3] * coef_zeta[1] + coef_z[3] * coef_zeta[2]);

        coef_tau[0] = 0.0;
        coef_tau[1] = 0.0;
        coef_tau[2] = 0.0;
        coef_tau[3] = 1.0;
    }


}