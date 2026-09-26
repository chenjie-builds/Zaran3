/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file CavityGasModel.cpp
 * \brief 高压气腔加载模型的实现。
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#include "CavityGasModel.h"
#include "Log.h"
#include "CommonPara.h"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace zaran
{
    namespace
    {
        /// @brief 把角度差折到 (−π, π]
        double WrapToPi(double angle)
        {
            while (angle > PI)  angle -= 2.0 * PI;
            while (angle <= -PI) angle += 2.0 * PI;
            return angle;
        }
    }

    void CavityGasModel::Initialize(const dynamic_array<DEMParticle>& particles,
                                    const CavityGasOptions& options,
                                    double lattice_spacing)
    {
        m_options = options;
        m_thickness = options.thickness > 0.0 ? options.thickness : 1.0;
        m_ring.clear();
        m_theta.clear();
        m_dtheta.clear();
        m_rho.clear();
        m_chord.clear();
        m_normal.clear();
        m_vented = false;
        m_vent_time = -1.0;
        m_gas_work = 0.0;
        m_external_work = 0.0;
        m_force_sum = 0.0;

        if (!m_options.enabled)
        {
            return;
        }
        if (!(m_options.cavity_radius > 0.0) || !(m_options.pressure_initial > 0.0))
        {
            Log::warn("CavityGasModel: gas_cavity_radius={} 或 gas_pressure_initial={} 非正，"
                "气腔加载被关闭", m_options.cavity_radius, m_options.pressure_initial);
            m_options.enabled = false;
            return;
        }

        const double a0 = m_options.cavity_radius;
        const double cx = m_options.center.x();
        const double cy = m_options.center.y();

        // --- 键合格距：优先用参数，否则退化为 2·中位粒子半径 ---
        double spacing = lattice_spacing;
        if (!(spacing > 0.0))
        {
            std::vector<double> radii;
            radii.reserve(particles.size());
            for (const auto& p : particles)
            {
                if (p.radius > 0.0) radii.push_back(p.radius);
            }
            if (!radii.empty())
            {
                std::nth_element(radii.begin(), radii.begin() + radii.size() / 2, radii.end());
                spacing = 2.0 * radii[radii.size() / 2];
            }
        }
        if (!(spacing > 0.0))
        {
            spacing = a0 * 0.1;
        }

        // --- 识别边界环：d_i ≤ a0 + r_p + tol·L0，且 d_i ≤ ratio·a0 ---
        const double shell = a0 + m_options.shell_tolerance * spacing;
        const double upper = m_options.max_ring_radius_ratio > 0.0
            ? m_options.max_ring_radius_ratio * a0 : 1.0e30;
        for (index_type i = 0; i < particles.size(); ++i)
        {
            const DEMParticle& p = particles[i];
            if (!p.active) continue;
            const double dx = p.pos.x() - cx;
            const double dy = p.pos.y() - cy;
            const double d = std::sqrt(dx * dx + dy * dy);
            if (d <= a0 + p.radius + m_options.shell_tolerance * spacing
                && d <= upper)
            {
                m_ring.push_back(i);
            }
        }
        if (m_ring.size() < 3)
        {
            Log::warn("CavityGasModel: 只识别到 {} 个边界粒子，无法构成闭合气腔，加载被关闭",
                m_ring.size());
            m_options.enabled = false;
            return;
        }

        // --- 按角度排序（环序一旦确定，之后固定不变）---
        std::sort(m_ring.begin(), m_ring.end(),
            [&](index_type ia, index_type ib)
            {
                const double aa = std::atan2(particles[ia].pos.y() - cy,
                                             particles[ia].pos.x() - cx);
                const double ab = std::atan2(particles[ib].pos.y() - cy,
                                             particles[ib].pos.x() - cx);
                return aa < ab;
            });

        const index_type n = static_cast<index_type>(m_ring.size());
        m_theta.assign(n, 0.0);
        m_dtheta.assign(n, 0.0);
        m_rho.assign(n, 0.0);
        m_chord.assign(n, 0.0);
        m_normal.assign(3 * n, 0.0);

        double min_rho = 1.0e30;
        double max_rho = -1.0e30;
        for (index_type i = 0; i < n; ++i)
        {
            const DEMParticle& p = particles[m_ring[i]];
            const double dx = p.pos.x() - cx;
            const double dy = p.pos.y() - cy;
            const double d = std::sqrt(dx * dx + dy * dy);
            m_rho[i] = std::max(0.0, d - p.radius);
            m_theta[i] = std::atan2(dy, dx);
            min_rho = std::min(min_rho, m_rho[i]);
            max_rho = std::max(max_rho, m_rho[i]);
        }

        // --- 初始几何面积：内切点多边形的有向面积绝对值 ---
        double twice_area = 0.0;
        for (index_type i = 0; i < n; ++i)
        {
            const index_type j = (i + 1) % n;
            const double xi = m_rho[i] * std::cos(m_theta[i]);
            const double yi = m_rho[i] * std::sin(m_theta[i]);
            const double xj = m_rho[j] * std::cos(m_theta[j]);
            const double yj = m_rho[j] * std::sin(m_theta[j]);
            twice_area += xi * yj - xj * yi;
        }
        const double area0 = 0.5 * std::fabs(twice_area);

        m_volume_0 = area0 * m_thickness;
        m_volume = m_volume_0;
        const double gamma = std::max(1.0 + 1.0e-6, m_options.polytropic_index);
        m_internal_energy_0 = m_options.pressure_initial * m_volume_0 / (gamma - 1.0);
        m_internal_energy = m_internal_energy_0;
        m_pressure = m_options.pressure_initial;
        m_area_geom = area0;

        Log::info("CavityGasModel: 气腔边界 {} 个粒子，初始面积 {:.6e} m^2，"
            "等效半径 {:.6e} m（输入 a0 = {:.6e} m），ρ ∈ [{:.6e}, {:.6e}] m",
            n, area0, std::sqrt(area0 / PI), a0, min_rho, max_rho);
        Log::info("CavityGasModel: p0 = {:.6e} Pa，gamma = {:.4f}，"
            "初始气体内能 U0 = {:.6e} J，V0 = {:.6e} m^3",
            m_options.pressure_initial, gamma, m_internal_energy_0, m_volume_0);
        if (m_options.ramp_time > 0.0)
        {
            Log::info("CavityGasModel: 升压时间 {:.6e} s（压力按 p0·t/t_ramp 线性上升）",
                m_options.ramp_time);
        }
    }

    void CavityGasModel::UpdateRingGeometry(const dynamic_array<DEMParticle>& particles)
    {
        const index_type n = static_cast<index_type>(m_ring.size());
        if (n < 3) return;

        const double cx = m_options.center.x();
        const double cy = m_options.center.y();

        // 1) 角位置与内切半径
        for (index_type i = 0; i < n; ++i)
        {
            const DEMParticle& p = particles[m_ring[i]];
            const double dx = p.pos.x() - cx;
            const double dy = p.pos.y() - cy;
            const double d = std::sqrt(dx * dx + dy * dy);
            m_theta[i] = std::atan2(dy, dx);
            m_rho[i] = std::max(0.0, d - p.radius);
            m_normal[3 * i + 0] = d > 1.0e-300 ? dx / d : 1.0;
            m_normal[3 * i + 1] = d > 1.0e-300 ? dy / d : 0.0;
            m_normal[3 * i + 2] = 0.0;
        }

        // 2) 沿环序单调展开（消除 atan2 的 ±π 跳变）。
        //    这样 θ_0..θ_{n−1} 单调递增、跨度略小于 2π；环的"接缝"由下面 Δθ 的
        //    cyclic 取法（θ_{−1} = θ_{n−1}−2π、θ_n = θ_0+2π）自然闭合，
        //    且 Σ_i Δθ_i ≡ 2π 恒成立。
        for (index_type i = 1; i < n; ++i)
        {
            m_theta[i] = m_theta[i - 1] + WrapToPi(m_theta[i] - m_theta[i - 1]);
        }

        // 3) 角向分片宽度 Δθ_i = ½(θ_{i+1} − θ_{i−1})，接缝处 cyclic
        m_ring_max_gap = 0.0;
        for (index_type i = 0; i < n; ++i)
        {
            const double theta_prev = (i == 0) ? m_theta[n - 1] - 2.0 * PI : m_theta[i - 1];
            const double theta_next = (i + 1 == n) ? m_theta[0] + 2.0 * PI : m_theta[i + 1];
            m_dtheta[i] = 0.5 * (theta_next - theta_prev);
            m_ring_max_gap = std::max(m_ring_max_gap, theta_next - m_theta[i]);
            m_chord[i] = m_rho[i] * m_dtheta[i];
        }

        // 4) 几何面积（与追踪面积交叉核对用）
        double twice_area = 0.0;
        for (index_type i = 0; i < n; ++i)
        {
            const index_type j = (i + 1) % n;
            const double xi = m_rho[i] * std::cos(m_theta[i]);
            const double yi = m_rho[i] * std::sin(m_theta[i]);
            const double xj = m_rho[j] * std::cos(m_theta[j]);
            const double yj = m_rho[j] * std::sin(m_theta[j]);
            twice_area += xi * yj - xj * yi;
        }
        m_area_geom = 0.5 * std::fabs(twice_area);
    }

    double CavityGasModel::CurrentPressure(double time) const
    {
        const double gamma = std::max(1.0 + 1.0e-6, m_options.polytropic_index);
        if (m_options.ramp_time > 0.0 && time < m_options.ramp_time)
        {
            return m_options.pressure_initial * std::max(0.0, time) / m_options.ramp_time;
        }
        const double volume = std::max(1.0e-300, m_volume);
        return std::max(0.0, (gamma - 1.0) * m_internal_energy / volume);
    }

    void CavityGasModel::ApplyPressureForce(dynamic_array<DEMParticle>& particles,
                                           double time)
    {
        m_force_sum = 0.0;
        if (!m_options.enabled || m_vented || m_ring.empty())
        {
            return;
        }
        UpdateRingGeometry(particles);

        double pressure = CurrentPressure(time);
        if (m_options.pressure_cap > 0.0)
        {
            pressure = std::min(pressure, m_options.pressure_cap);
        }
        m_pressure = pressure;
        if (!(pressure > 0.0))
        {
            return;
        }

        const index_type n = static_cast<index_type>(m_ring.size());
        for (index_type i = 0; i < n; ++i)
        {
            const double segment_force = pressure * m_chord[i] * m_thickness;
            Eigen::Vector3d f(segment_force * m_normal[3 * i + 0],
                              segment_force * m_normal[3 * i + 1],
                              segment_force * m_normal[3 * i + 2]);
            particles[m_ring[i]].force += f;
            m_force_sum += f.norm();
        }
    }

    void CavityGasModel::UpdateState(const dynamic_array<DEMParticle>& particles,
                                     double dt, double time)
    {
        if (!m_options.enabled || m_vented || m_ring.empty() || !(dt > 0.0))
        {
            return;
        }
        const index_type n = static_cast<index_type>(m_ring.size());
        const double gamma = std::max(1.0 + 1.0e-6, m_options.polytropic_index);

        // --- 体积变化：用与受力**完全相同**的 ℓ_i / n̂_i ---
        // 半隐式 Euler 下位移就是当前速度乘 dt，故 W = Σ F·Δu = p·ΔV 逐位成立。
        double delta_area = 0.0;
        for (index_type i = 0; i < n; ++i)
        {
            const Eigen::Vector3d du = particles[m_ring[i]].vel * dt;
            const double un = du.x() * m_normal[3 * i + 0] + du.y() * m_normal[3 * i + 1];
            delta_area += m_chord[i] * un;
        }
        const double delta_volume = delta_area * m_thickness;
        const double work = m_pressure * delta_volume;

        m_volume += delta_volume;
        m_internal_energy -= work;
        m_gas_work += work;
        if (m_volume < 1.0e-300)
        {
            m_volume = 1.0e-300;
        }

        // --- 升压段：压力按时间给定，注入的能量单独记账 ---
        const double next_time = time + dt;
        if (m_options.ramp_time > 0.0 && next_time <= m_options.ramp_time)
        {
            const double target_pressure =
                m_options.pressure_initial * next_time / m_options.ramp_time;
            const double target_energy = target_pressure * m_volume / (gamma - 1.0);
            if (target_energy > m_internal_energy)
            {
                m_external_work += target_energy - m_internal_energy;
                m_internal_energy = target_energy;
            }
        }

        if (m_internal_energy < 0.0)
        {
            m_internal_energy = 0.0;
        }
        m_pressure = std::max(0.0, (gamma - 1.0) * m_internal_energy / m_volume);
        if (m_options.pressure_cap > 0.0)
        {
            m_pressure = std::min(m_pressure, m_options.pressure_cap);
        }

        // --- 泄压：气腔已经扩到自由面（面积超过判据）⇒ 气体逃逸 ---
        if (m_options.vent_area_ratio > 0.0 && m_options.disk_radius > 0.0)
        {
            const double limit = m_options.vent_area_ratio * PI
                * m_options.disk_radius * m_options.disk_radius;
            if (m_volume / m_thickness >= limit)
            {
                m_vented = true;
                m_vent_time = next_time;
                m_pressure = 0.0;
                m_internal_energy = 0.0;
                Log::info("CavityGasModel: 气腔面积达到 {:.6e} m^2（判据 {:.6e} m^2），"
                    "t = {:.6e} s 起判定气体泄出，加载结束",
                    m_volume / m_thickness, limit, next_time);
            }
        }

        // 边界环裂缝预警（环不再闭合时几何离散退化）
        if (!m_warned_ring && m_ring_max_gap > 1.5 * 2.0 * PI / static_cast<double>(n))
        {
            Log::warn("CavityGasModel: 边界环最大角间隙 {:.4f} rad，"
                "已超过环平均间隙的 1.5 倍（{} 个粒子），气腔壁面开始离散化退化",
                m_ring_max_gap, n);
            m_warned_ring = true;
        }
    }
} // namespace zaran
