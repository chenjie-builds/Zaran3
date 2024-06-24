#include "Residual.h"
#include"Log.h"
namespace zaran
{

    ResInfo::ResInfo(int var_num) : m_var_num(var_num)
    {
        m_max_res = new double[m_var_num];
        m_ave_res = new double[m_var_num];
        m_max_res_coord = new double[m_var_num*3];
        m_max_res_idx = new int[m_var_num];
    }
    ResInfo::~ResInfo()
    {
        delete[] m_max_res;
        delete[] m_ave_res;
        delete[] m_max_res_coord;
        delete[] m_max_res_idx;
    }

    double ResInfo::GetMaxRes(int iVar)
    {
        return m_max_res[iVar];
    }

    double ResInfo::GetAveRes(int iVar)
    {
        return m_ave_res[iVar];
    }

    const double* ResInfo::GetMaxResCoord(int iVar) const
    {
        return m_max_res_coord + 3 * iVar;
    }

    void ResInfo::SetMaxRes(int iVar, double value)
    {
        m_max_res[iVar] = value;
    }

    void ResInfo::SetAveRes(int iVar, double value)
    {
        m_ave_res[iVar] = value;
    }

    void ResInfo::SetMaxResCoord(int iVar, const double* coord)
    {
        for (int i_dim = 0; i_dim < 3; i_dim++)
        {
            m_max_res_coord[3 * iVar + i_dim] = coord[i_dim];
        }
    }

}