#pragma once
namespace zaran
{

    /// @brief CFL数计算类
    class CFL
    {
    public:
        /// @brief 构造函数
        /// @param cfl_min 最小CFL数
        /// @param cfl_max 最大CFL数
        /// @param start_step 起始步数
        /// @param grow_step 增长步数，即CFL数从最小值增长到最大值所需的步数
        /// @param reduce_factor 缩减因子,用于减小CFL数
        CFL(double cfl_min, double cfl_max, int start_step, int grow_step , double reduce_factor );
        ~CFL();
       void SetCFL(double cfl_min, double cfl_max, int start_step, int grow_step, double reduce_factor);
       /// @brief 减小CFL数
       /// @details 将CFL数减小到原来的reduce_factor倍
       void ReduceCFL();
       /// @brief 获取CFL数
       /// @param current_step 当前步数
        double GetCFL(int current_step);
    private:
        /// @brief CFL最小值
        double m_cfl_min;
        /// @brief CFL最大值
        double m_cfl_max;
        /// @brief 起始步数
        int m_start_step;
        /// @brief 增长步数
        int m_grow_step;
        /// @brief 缩减因子,用于减小CFL数
        double m_reduce_factor;
    };
}