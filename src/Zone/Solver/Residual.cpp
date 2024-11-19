#include "Residual.h"
#include"Log.h"
namespace zaran
{

    ResInfo::ResInfo(int equ_num) : m_equ_num(equ_num)
    {
        m_res_inf = new double[m_equ_num];
        m_ref_L2 = new double[m_equ_num];
        m_res_inf_coord = new double[m_equ_num*3];
        m_res_inf_idx = new int[m_equ_num];
    }
    ResInfo::~ResInfo()
    {
        delete[] m_res_inf;
        delete[] m_ref_L2;
        delete[] m_res_inf_coord;
        delete[] m_res_inf_idx;
    }

    double ResInfo::GetInfNorm(int iVar)
    {
        return m_res_inf[iVar];
    }

    double ResInfo::GetL2Norm(int iVar)
    {
        return m_ref_L2[iVar];
    }

    const double* ResInfo::GetInfNormCoord(int iVar) const
    {
        return m_res_inf_coord + 3 * iVar;
    }

    int ResInfo::GetInfNormIdx(int iVar) const
    {
        return m_res_inf_idx[iVar];
    }

    void ResInfo::SetInfNorm(int iVar, double value)
    {
        m_res_inf[iVar] = value;
    }

    void ResInfo::SetL2Norm(int iVar, double value)
    {
        m_ref_L2[iVar] = value;
    }

    void ResInfo::SetInfNormCoord(int iVar, const double* coord)
    {
        for (int i_dim = 0; i_dim < 3; i_dim++)
        {
            m_res_inf_coord[3 * iVar + i_dim] = coord[i_dim];
        }
    }

    void ResInfo::SetInfNormIdx(int iVar, int idx)
    {
        m_res_inf_idx[iVar] = idx;
    }

}