/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file MappingSimulation.h
 * \brief MappingSimulation: standalone case runner for particle -> grid mapping.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "BasicType.h"
#include "ParticleGridMapper.h"
namespace zaran
{
    /// @brief 粒子 → 网格映射的独立算例控制器（`task.simulation = "MAPPING"`）
    ///
    /// 它**不做时间推进**：把 `particles.csv`（与 DEM 同一个文件格式）里的粒子，
    /// 按一串逐级加密的均匀网格映射成固相体积分数 α，逐帧写成 Tecplot ASCII，
    /// 供在 Tecplot 里看 α 云图、并与粒子点云叠加对比。
    ///
    /// 存在的意义：阶段 2 的映射器此前只有独立的基准驱动（命令行参数喂参数），
    /// 不能像其它算例那样"给一个目录就能跑"。这个类补上那一环，顺带把
    /// 体积守恒残差、α 范围、光滑性指标随网格的变化打进日志与报告文件，
    /// 让 Δx/d_p 的方法选择准则可以当场量出来。
    class MappingSimulation
    {
    public:
        explicit MappingSimulation(std::string work_dir);
        ~MappingSimulation() = default;

        /// @brief 读参数 → 读粒子 → 逐级映射 → 输出 Tecplot + 报告
        void Run();

    protected:
        /// @brief 从 `[mapping]` 段读全部参数并对取值做校验
        void ReadParameters();
        /// @brief 读粒子文件（复用 `ReadDEMParticle`，与 DEM 输入格式完全一致）
        void ReadParticles();
        /// @brief 按 `[mapping]` 的 box/n 构造第 level 级的网格
        ParticleGridMapper::GridSpec MakeGrid(int level) const;
        /// @brief 输出一个帧（α 场 + 可选的粒子点云）
        void WriteFrame(int frame_index, int level, ParticleGridMapper::Method method,
            const ParticleGridMapper::GridSpec& grid,
            const dynamic_array<double>& alpha,
            const ParticleGridMapper::Stats& stats,
            bool append_header);
        /// @brief 写机器可读的汇总报告（`result/mapping_report.csv`）
        void WriteReport() const;

        struct FrameRecord
        {
            int         index = 0;      ///< 帧号（文件后缀）
            std::string method;         ///< "subcell" / "exact"
            int         n_cell = 0;     ///< 每方向单元数
            double      dx = 0.0;       ///< 网格间距
            double      dx_over_dp = 0.0; ///< Δx/d_p
            double      residual = 0.0; ///< |ΣαV - ΣV_p| / ΣV_p
            double      alpha_min = 0.0;
            double      alpha_max = 0.0;
            double      max_jump = 0.0; ///< 相邻单元 |Δα| 最大值
            double      smooth_shift = 0.0;
            double      l1_vs_exact = -1.0; ///< 与同级 Exact 结果的 L1（仅 subcell 帧有效）
            std::string file;           ///< 输出文件名（不含目录）
        };

        std::string m_work_dir;
        std::string m_result_folder;

        // --- [mapping] 参数 ---
        std::string m_particle_file = "particles.csv";
        double m_lo[3] = { -0.03, -0.03, -0.03 };   ///< 物理计算域（含半格）
        double m_hi[3] = { 0.03, 0.03, 0.03 };
        int    m_dim = 3;
        int    m_n0 = 4;                ///< 最粗一级每方向单元数
        int    m_refine = 2;            ///< 每级加密倍数
        int    m_levels = 4;            ///< 输出级数
        std::string m_method = "both";  ///< subcell | exact | both
        int    m_sub = 2;               ///< SubCell 每方向细分数
        int    m_smooth_passes = 0;
        double m_smooth_blend = 1.0;
        bool   m_write_particles = true;

        // --- 运行时状态 ---
        dynamic_array<ParticleGridMapper::Particle> m_particles;
        double m_particle_volume = 0.0;     ///< Σ V_p
        double m_dp = 0.0;                  ///< 代表性粒径（2 × 最小半径）
        int    m_outside_count = 0;         ///< 中心落在计算域外的粒子数
        dynamic_array<FrameRecord> m_records;
    };
} // namespace zaran
