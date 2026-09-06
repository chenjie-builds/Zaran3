/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file ReadDEMParticle.h
 * \brief Read DEM particle initial state from a CSV file.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "DEMParticle.h"
#include <string>
#include <vector>
namespace zaran
{
    /// @brief 从 CSV 文件读取粒子初始状态
    ///
    /// CSV 文件格式（首行为表头，忽略或兼容）：
    ///   id, group, radius, mass, px, py, pz, vx, vy, vz, ox, oy, oz
    /// 可选列：motion_type, young_modulus, poisson_ratio, friction_coeff, restitution_coeff
    /// motion_type: 0=固定，1=动力学，2=按输入速度运动的运动学粒子。
    ///
    /// mass 为 0 时，求解器根据 dem.density 计算质量。
    class ReadDEMParticle
    {
    public:
        ReadDEMParticle() = default;
        ~ReadDEMParticle() = default;

        /// @brief 读取 CSV 文件，填充粒子列表
        /// @param filename 文件路径（绝对路径）
        /// @param particles 输出粒子数组
        void ReadCSV(const std::string& filename, std::vector<DEMParticle>& particles) const;
    };
} // namespace zaran
