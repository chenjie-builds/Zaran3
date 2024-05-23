#include "GradWLSQ.h"
#include"Log.h"
using namespace zaran;
zaran::GradWLSQ::GradWLSQ(GridFN* grid):m_grid(grid)
{
        m_grid = grid;
    auto node = grid->GetNodeTopo();
    int node_num = node->GetNodeNum();
    m_omega = new double** [node_num];
    for (int iNode = 0;iNode < node_num;iNode++)
    {
        auto neighbors = node->GetNeighborNode(iNode);
        int neighbor_num = node->GetNeighborNodeNum(iNode);
        m_omega[iNode] = new double* [neighbor_num];
        for (int iNeigh = 0;iNeigh < neighbor_num;iNeigh++)
        {
            m_omega[iNode][iNeigh] = new double[3];
            for (int iDim = 0;iDim < 3;iDim++)
            {
                m_omega[iNode][iNeigh][iDim] = 0.0;
            }
        }
    }
}
GradWLSQ::~GradWLSQ()
{
    auto node = m_grid->GetNodeTopo();
    int node_num = node->GetNodeNum();
    for (int iNode = 0;iNode < node_num;iNode++)
    {
        auto neighbors = node->GetNeighborNode(iNode);
        int neighbor_num = node->GetNeighborNodeNum(iNode);
        for (int iNeigh = 0;iNeigh < neighbor_num;iNeigh++)
        {
            delete[] m_omega[iNode][iNeigh];
        }
        delete[] m_omega[iNode];
    }
    delete[] m_omega;
}

void GradWLSQ::CalcGradient(GridFN* grid, const double* data, double* grad_x, double* grad_y, double* grad_z)
{
    if (grid != m_grid)
    {
        Log::error("Gradient::CalcGradient: grid is not matched");
        return;
    }
    auto node = grid->GetNodeTopo();
    int node_num = node->GetNodeNum();
    for (int iNode = 0;iNode < node_num;iNode++)
    {
        auto neighbors = node->GetNeighborNode(iNode);
        int neighbor_num = node->GetNeighborNodeNum(iNode);
        grad_x[iNode] = grad_y[iNode] = grad_z[iNode] = 0.0;
        for (int iNeigh = 0;iNeigh < neighbor_num;iNeigh++)
        {
            int iNeighNode = neighbors[iNeigh];
            grad_x[iNode] += data[iNeighNode] * m_omega[iNode][iNeigh][0];
            grad_y[iNode] += data[iNeighNode] * m_omega[iNode][iNeigh][1];
            grad_z[iNode] += data[iNeighNode] * m_omega[iNode][iNeigh][2];
        }
    }
}