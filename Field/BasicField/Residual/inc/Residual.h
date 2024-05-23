#pragma once
#include"GridFNFDM.h"
namespace zaran
{
    class ResAnalyzerFN
    {
    public:
        ResAnalyzerFN(GridFN* grid, double** res, int var_num);
        ~ResAnalyzerFN();
        void Analyze();
        double GetMaxResidual(int iVar);
        double GetAveResidual(int iVar);
        double GetMaxResidualX(int iVar);
        double GetMaxResidualY(int iVar);
        double GetMaxResidualZ(int iVar);
    private:
        int m_var_num;
        GridFN* m_grid;
        double** m_res;
        double* m_max_residual;
        double* m_ave_residual;
        double* m_max_residual_x;
        double* m_max_residual_y;
        double* m_max_residual_z;
        int* m_max_residual_node;
    };
}