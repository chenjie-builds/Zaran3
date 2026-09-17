/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file DEMSolver.h
 * \brief DEMSolver class: performs one DEM time step.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "FieldSolver.h"
#include "DEMFieldData.h"
#include "DEMSolverParam.h"
#include "ContactModel.h"
#include <algorithm>
#include <cstdint>
namespace zaran
{
    /// @brief DEM 求解器，继承 FieldSolver 基类
    ///
    /// 每调用一次 Solve() 推进一个 DEM 时间步：
    ///   ZeroForce → ContactDetection → CalcContactForce → CalcWallForce
    ///   → CalcGravity → Integrate
    class DEMSolver : public FieldSolver
    {
    public:
        DEMSolver(index_type index,
                  const std::string& name,
                  shared_ptr<DEMSolverParam> para,
                  shared_ptr<DEMFieldData>   dem_data);
        ~DEMSolver() override = default;

        // --- FieldSolver 接口 ---
        void Init()        override;
        void Solve()       override;
        void Preprocess()  override;
        void Postprocess() override;
        void InitSolver()  override;

        void InitField()                             override;
        void BoundaryCondition()                     override {}
        void UpdateField()                           override {}
        void BackupField(std::string& back_folder)   override;

        // --- 专用接口 ---
        /// @brief 备份粒子状态到文件（CSV）
        void BackupField(const std::string& back_folder) const;

        DEMFieldData*    GetDEMData()  const { return m_dem_data.get(); }
        DEMSolverParam*  GetDEMParam() const;

    protected:
        // 主流程各子步
        void ZeroForce();
        void ContactDetection();
        void CalcContactForce();
        void CalcWallForce();
        void CalcBondForce();
        void CalcThermalReaction();
        // --- 重力与积分 ---
        // 在启用反应功能时，该函数会根据粒子的 reaction_progress 和 gas_pressure 计算并施加气相压力力，将力以相反方向分别作用到成对颗粒上。
        // 若存在有效的气体 Voronoi 面，则按 Voronoi 面面积与局部开口比例计算；否则退化为按键合面的 conduction_area 和 rest_length 进行近似，并结合参考长度/当前距离修正作用力。
        void CalcGasPressureForce();
        double JwlPressure(const DEMParticle& particle, double gas_volume) const;
        void UpdateGasSolidState(DEMParticle& particle) const;
        void UpdateGasVoronoiMesh(bool force = false);
        void CalcGravity();
        void Integrate();

        // --- 接触重叠限制：不允许两个粒子重叠过近（参考 LSM linear_elastic_spring_mod.f90 /
        //     force_between_discrete_lattices.f90 的 rn_limit / rn_rebound）---
        /// @brief 深压缩刚性回弹。当 dist < ratio·(r_a+r_b) 且两端仍相向运动时，
        /// 基于动量守恒按（准）弹性碰撞直接改写两端的法向速度；
        /// 规定运动（kinematic / 非 dynamic）粒子按无限质量处理，自身速度不变。
        /// @return 本次是否发生了回弹（供诊断/统计）
        bool ApplyContactRebound(DEMParticle& pa, DEMParticle& pb, const DEMContact& contact) const;
        /// @brief 过深压缩时非物理放大法向排斥力：F_n *= (ratio·(r_a+r_b)/dist)^20。
        /// 须在 CalcNormalForce 之后、CalcTangentialForce 之前调用，
        /// 使切向库仑摩擦上限用的是放大后的法向力（与 LSM 一致）。
        void AmplifyDeepOverlapForce(DEMContact& contact, double sum_radius) const;

        // --- 均匀网格（链表式单元表）：接触检测与初始弹簧网络共用 ---
        /// @brief 建立单元表；cell 为期望单元边长（会放大以限制单元总数），返回实际边长
        double BuildUniformGrid(const dynamic_array<DEMParticle>& particles, double cell);
        /// @brief 查询点所属单元索引（须先调用 BuildUniformGrid）
        index_type GridCellOf(const Eigen::Vector3d& q) const;

        // --- 可选功能：初始弹簧连接网络（脆性材料断裂/破坏）---
        /// @brief 在仿真开始时按几何邻近关系自动建立弹簧连接，追加到键合列表。
        /// 与逐步接触检测相互独立：只在 t=0 执行一次，之后每步仅按 CalcBondForce 计算。
        void BuildSpringNetwork();

        // --- 可选功能：逐键强度异质性（脆性材料的离散强度）---
        /// @brief 为每条键分配 strength_scale（Weibull 抽样，均值归一为 1）。
        /// 在 InitField 中所有键（文件给定 + 自动建链）就绪后执行一次；
        /// dem.bond_weibull_modulus ≤ 0 时不做任何改动，保持均质（全部为 1）。
        void AssignBondStrengthScales();

    private:
        /// @brief 接触对唯一键：打包的 64 位整数（最高位区分"粒子-墙"）。
        /// 旧实现用 std::tuple 作为 std::map/std::set 的键，每个接触每步都要做
        /// 红黑树查找 + 节点分配，是万级粒子下的主要开销之一。
        using ContactKey = std::uint64_t;

        /// @brief 接触历史（切向弹簧位移）开放寻址哈希表。
        /// 采用"双缓冲 + 代标记"：每步把上一步的表作为只读历史，写进当前表；
        /// 步末交换，未出现的接触自然被淘汰。重置为 O(1)（只递增代标记），
        /// 全程无节点分配。
        struct ContactHistoryTable
        {
            struct Entry
            {
                ContactKey      key     = 0;
                Eigen::Vector3d delta_t = Eigen::Vector3d::Zero();
                std::uint32_t   gen     = 0; ///< 等于本表 gen 时该槽被占用
            };

            dynamic_array<Entry> slots;
            std::uint32_t        gen   = 0;
            index_type           count = 0;

            static std::size_t Hash(ContactKey x)
            {
                x ^= x >> 33; x *= 0xff51afd7ed558ccdULL;
                x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL;
                x ^= x >> 33;
                return static_cast<std::size_t>(x);
            }

            /// @brief O(1) 清空（代标记递增；回绕时做一次真正的清零）
            void Reset()
            {
                if (++gen == 0)
                {
                    for (auto& e : slots) e.gen = 0;
                    gen = 1;
                }
                count = 0;
            }

            void Reserve(index_type need)
            {
                index_type cap = 16;
                while (cap < need * 2) cap <<= 1;
                if (static_cast<index_type>(slots.size()) < cap)
                {
                    slots.assign(cap, Entry{});
                    gen = 0;
                    Reset();
                }
            }

            const Eigen::Vector3d* Find(ContactKey key) const
            {
                if (slots.empty()) return nullptr;
                const std::size_t mask = slots.size() - 1;
                std::size_t s = Hash(key) & mask;
                while (slots[s].gen == gen)
                {
                    if (slots[s].key == key) return &slots[s].delta_t;
                    s = (s + 1) & mask;
                }
                return nullptr;
            }

            void InsertRaw(ContactKey key, const Eigen::Vector3d& v)
            {
                const std::size_t mask = slots.size() - 1;
                std::size_t s = Hash(key) & mask;
                while (slots[s].gen == gen)
                {
                    if (slots[s].key == key) { slots[s].delta_t = v; return; }
                    s = (s + 1) & mask;
                }
                slots[s].key = key;
                slots[s].delta_t = v;
                slots[s].gen = gen;
                ++count;
            }

            void Insert(ContactKey key, const Eigen::Vector3d& v)
            {
                if (slots.empty() || (count + 1) * 2 > slots.size()) Grow();
                InsertRaw(key, v);
            }

            void Grow()
            {
                dynamic_array<Entry> old = std::move(slots);
                const std::uint32_t old_gen = gen;
                slots.assign(std::max<std::size_t>(16, old.size() * 2 + 8), Entry{});
                gen = 1;
                count = 0;
                for (const auto& e : old)
                    if (e.gen == old_gen) InsertRaw(e.key, e.delta_t);
            }
        };

        static ContactKey MakeContactKey(const DEMContact& contact);
        void RestoreContactHistory(DEMContact& contact) const;
        void SaveContactHistory(const DEMContact& contact);
        /// @brief 步末推进接触历史：交换双缓冲并清空当前表
        void AdvanceContactHistory();

        /// @brief 上一步的接触历史（只读，用于恢复 delta_t）
        ContactHistoryTable m_history_prev;
        /// @brief 当前步写入的接触历史
        ContactHistoryTable m_history_cur;

        // --- 均匀网格（cell list）接触检测缓冲，跨步复用，零每步分配 ---
        dynamic_array<index_type> m_grid_head; ///< 每个网格单元的首个粒子槽（-1 表示空）
        dynamic_array<index_type> m_grid_next; ///< 同单元内下一个粒子槽
        index_type         m_grid_nx  = 0;
        index_type         m_grid_ny  = 0;
        index_type         m_grid_nz  = 0;
        index_type         m_grid_nyz = 0;
        Eigen::Vector3d    m_grid_lo  = Eigen::Vector3d::Zero();
        double             m_grid_cell = 0.0;

        shared_ptr<DEMFieldData>   m_dem_data;
        unique_ptr<ContactModel>   m_contact_model;

        struct GasVoronoiCell
        {
            std::vector<Eigen::Vector2d> vertices;
            double area = 0.0;
        };
        struct GasVoronoiFace
        {
            index_type idx_a = 0;
            index_type idx_b = 0;
            double length = 0.0;
            double area = 0.0;
        };
        std::vector<GasVoronoiCell> m_gas_voronoi_cells;
        std::vector<GasVoronoiFace> m_gas_voronoi_faces;
        std::uint64_t m_dem_step = 0;
        bool m_voronoi_valid = false;
        double m_lattice_spacing = 0.0;

        // --- 热化学/反应每步复用的临时缓冲 ---
        // 原实现在每步构造 7 个大小为 N 的 std::vector，产生持续的
        // 分配/释放；改为成员缓冲 + assign，容量复用、不再分配。
        dynamic_array<double> m_tmp_energy_delta;
        dynamic_array<double> m_tmp_surface_delta;
        dynamic_array<double> m_tmp_core_delta;
        dynamic_array<double> m_tmp_length_ratio_sum;
        dynamic_array<int>    m_tmp_length_ratio_count;
        dynamic_array<double> m_tmp_gas_pressure_sum;
        dynamic_array<int>    m_tmp_gas_pressure_count;
        /// @brief 本步各粒子由【接触法向阻尼 + 接触切向摩擦 + 键合阻尼】不可逆耗散
        /// 的能量 (J)。在 ZeroForce 中清零/扩容，由 CalcBondForce / CalcContactForce /
        /// CalcWallForce 累加（每对按 0.5/0.5 分给两端），在 CalcThermalReaction 中
        /// 并入温度的能量账 energy_delta。
        dynamic_array<double> m_tmp_dissipation;

        /// @brief 备份输出的复用文本缓冲（先拼完整内容再一次性写盘，
        /// 避免逐字段 `ostream <<` 的格式化开销）
        mutable std::string m_out_buffer;
    };
} // namespace zaran
