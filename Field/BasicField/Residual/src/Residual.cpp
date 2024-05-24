#include "Residual.h"
#include"Log.h"
using namespace zaran;
ResAnalyzerFN::ResAnalyzerFN(GridFN* grid, double** res, int var_num) :m_grid(grid), m_res(res), m_var_num(var_num)
{
    m_max_residual = new double[m_var_num];
    m_ave_residual = new double[m_var_num];
    m_max_residual_x = new double[m_var_num];
    m_max_residual_y = new double[m_var_num];
    m_max_residual_z = new double[m_var_num];
    m_max_residual_node = new int[m_var_num];
}
ResAnalyzerFN::~ResAnalyzerFN()
{
    delete[] m_max_residual;
    delete[] m_ave_residual;
    delete[] m_max_residual_x;
    delete[] m_max_residual_y;
    delete[] m_max_residual_z;
    delete[] m_max_residual_node;
    delete[] m_res;
}

void ResAnalyzerFN::Analyze()
{
    auto node = m_grid->GetNode();
    int node_num = node->GetNodeNum();
    for (int iVar = 0; iVar < m_var_num; iVar++)
    {
        m_max_residual[iVar] = -LARGE_NUMBER;
        m_ave_residual[iVar] = 0;
    }
    double max_res = -LARGE_NUMBER;
    double ave_res = 0;
    int max_res_node = 0;
    double max_res_x = 0;
    double max_res_y = 0;
    double max_res_z = 0;
    for (int iVar = 0;iVar < m_var_num;iVar++)
    {
        max_res = -LARGE_NUMBER;
        ave_res = 0;
        max_res_node = -1;
#pragma omp parallel for reduction(max:max_res) reduction(+:ave_res)
        for (int iNode = 0; iNode < node_num; iNode++)
        {
            if (node->GetType(iNode) != NodeType::inner)
                continue;
            if (abs(m_res[iVar][iNode]) > max_res)
            {
                max_res = abs(m_res[iVar][iNode]);
                max_res_node = iNode;
                max_res_x = node->GetCoord(iNode)[0];
                max_res_y = node->GetCoord(iNode)[1];
                if (m_grid->GetDimension() == 3)
                {
                    max_res_z = node->GetCoord(iNode)[2];
                }
            }
            ave_res += m_res[iVar][iNode] * m_res[iVar][iNode];
        }
        ave_res = sqrt(ave_res / node_num);
        m_max_residual[iVar] = max_res;
        m_ave_residual[iVar] = ave_res;
        m_max_residual_node[iVar] = max_res_node;
        m_max_residual_x[iVar] = max_res_x;
        m_max_residual_y[iVar] = max_res_y;
        m_max_residual_z[iVar] = max_res_z;
    }

}
double ResAnalyzerFN::GetMaxResidual(int iVar)
{
    return m_max_residual[iVar];
}

double ResAnalyzerFN::GetAveResidual(int iVar)
{
    return m_ave_residual[iVar];
}

double ResAnalyzerFN::GetMaxResidualX(int iVar)
{
    return m_max_residual_x[iVar];
}

double ResAnalyzerFN::GetMaxResidualY(int iVar)
{
    return m_max_residual_y[iVar];
}

double ResAnalyzerFN::GetMaxResidualZ(int iVar)
{
    return m_max_residual_z[iVar];
}
