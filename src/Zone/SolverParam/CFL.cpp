#include "CFL.h"
#include "Log.h"
#include "ZaranError.h"
#include <iostream>
namespace zaran
{
    CFL::CFL(double cfl_min, double cfl_max, int start_step, int grow_step, double reduce_factor)
    {
        SetCFL(cfl_min, cfl_max, start_step, grow_step, reduce_factor);
    }

    CFL::~CFL()
    {
    }
    void CFL::SetCFL(double cfl_min, double cfl_max, int start_step, int grow_step, double reduce_factor)
    {
        if (cfl_min < 0 || cfl_max < 0 || start_step < 0 || grow_step < 0 || reduce_factor < 0)
        {
            Log::error("CFL number must be positive, please check the input parameters");
            Log::error("cfl_min = {}, cfl_max = {}, start_step = {}, grow_step = {}, reduce_factor = {}", cfl_min, cfl_max, start_step, grow_step, reduce_factor);
            throw ZaranError("CFL parameters must be positive");
        }
        if(cfl_min > cfl_max)
        {
            Log::error("cfl_min must be less than cfl_max, please check the input parameters");
            Log::error("cfl_min = {}, cfl_max = {}", cfl_min, cfl_max);
            throw ZaranError("cfl_min must be less than cfl_max");
        }
        if(reduce_factor > 1)
        {
            Log::error("reduce_factor must be less than 1, please check the input parameters");
            Log::error("reduce_factor = {}", reduce_factor);
            throw ZaranError("reduce_factor must be less than 1");
        }
        m_cfl_min = cfl_min;
        m_cfl_max = cfl_max;
        m_start_step = start_step;
        m_grow_step = grow_step;
        m_reduce_factor = reduce_factor;
    }
    void CFL::ReduceCFL()
    {
        m_cfl_max *= m_reduce_factor;
        if (m_cfl_max < m_cfl_min)
        {
            m_cfl_max = m_cfl_min;
        }
    }
    double CFL::GetCFL(int current_step)
    {
        if (current_step < m_start_step)
        {
            return m_cfl_min;
        }
        else if (current_step >= m_start_step + m_grow_step)
        {
            return m_cfl_max;
        }
        else
        {
            return m_cfl_min + (m_cfl_max - m_cfl_min) * (current_step - m_start_step) / m_grow_step;
        }
    }
} // namespace zaran