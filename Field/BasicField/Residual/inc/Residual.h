#pragma once
#include"GridFNFDM.h"
namespace zaran
{
    class ResInfo
    {
    public:
        ResInfo(int var_num);
        ~ResInfo();
        double GetMaxRes(int iVar);
        double GetAveRes(int iVar);
        const double* GetMaxResCoord(int iVar)const;
        void SetMaxRes(int iVar, double value);
        void SetAveRes(int iVar, double value);
        void SetMaxResCoord(int iVar, const double* coord);
    private:
        int m_var_num;
        // 最大残差
        double* m_max_res;
        // 平均残差
        double* m_ave_res;
        // 最大残差坐标
        double* m_max_res_coord;
        // 最大残差编号
        int* m_max_res_idx;
    };
}