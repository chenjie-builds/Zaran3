/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file DEMCFDSimulation.h
 * \brief One-way coupled DEM-CFD case runner (`task.solver = "EulerDEM"`).
 * \author Chen Jie.
 *
 *  阶段 3：DEM-CFD 单向耦合算例（见 docs/TWO_PHASE_STAGE3_DEM_CFD.md）
 *  ---------------------------------------------------------------
 *  流体侧**完全复用**两相 Euler 求解器（`EulerTwoPhaseStructUniform`，
 *  `task.solver` 只是多了个名字），DEM 侧完全复用粒子的几何与映射工具。本类只做三件事：
 *
 *    1. 把 particles.csv 里的粒子映射成气相体积分数 ε_g 并注入求解器
 *       （走 `VolumeFractionType::External`，粒子不动 ⇒ ε 只算一次）；
 *    2. 每个时间步之后，用同一套权重把网格上的 (ρ,u,p,μ) 与 ∇p 插值到粒子位置，
 *       经 `DragForceModel` 得到曳力 + 压力梯度力（**只算力，不反馈给流体**）；
 *    3. 把逐粒子的气态、Re、C_d、受力写成 `particle_force.csv` 与 Tecplot 点云。
 *
 *  时间推进不重写：`NSFieldSimulation` 的主循环 + 三个耦合钩子（PrepareCoupling /
 *  CouplingPreStep / CouplingPostStep）就是全部。**双向偶合**要做的正是在
 *  `CouplingPreStep` 里把上一时刻的相间反作用力散列进 `CalcSourceResidual`、
 *  并按子循环推进 DEM —— 本类的结构不需要变。
 *
 *  关于"为什么粒子不运动也要算时间步"
 *  ----------------------------------
 *  单向耦合下粒子的确不动，但**受力是时间的函数**：激波扫过粒子时曳力从 0 阶跃到
 *  某个量级，压力梯度力在激波/接触间断/稀疏波处才出现。所以本算例的输出是一张
 *  "力-时间历程"，它的正确性必须对照精确 Riemann 解来判，而不是只看末态。
 */

#pragma once
#include "BasicType.h"
#include "FieldSimulation.h"
#include "DEMCFDCoupler.h"
// 两相场（NSFieldStructTwoPhase）已经把求解器/参数/数据管理器三个头一起带进来了
#include "NSFieldStructTwoPhase.h"

#include <string>
#include <vector>

namespace zaran
{
	/// @brief DEM-CFD 单向耦合算例控制器。
	class DEMCFDSimulation : public NSFieldSimulation
	{
	public:
		explicit DEMCFDSimulation(const shared_ptr<FieldManager>& field_manager);
		virtual ~DEMCFDSimulation();

	protected:
		// ---------------- 三个耦合钩子 ----------------
		/// @brief 读参数/粒子 → 建耦合器 → 把 ε_g 注入两相数据管理器的 volume_fraction
		void PrepareCoupling() override;
		/// @brief 每个时间步之前：重建 ε（颗粒运动时）→ 评估受力 → 散列反作用力源项
		void CouplingPreStep(int iter) override;
		/// @brief 每个时间步之后：推进颗粒 → 在当前气态上重算受力 → 输出
		void CouplingPostStep(int iter) override;

	public:
		/// @brief 时间推进完全复用基类，只在收尾时补一份机器可读的耦合报告
		void SolveField() override;

	private:
		/// @brief 从场管理器取出两相求解器/参数/数据管理器（做严格类型校验）
		void BindField();
		/// @brief 读 `[coupling]` 段
		void ReadParameters();
		/// @brief 读 `particles.csv`（复用 DEM 的读取器，格式与 DEM 算例一致）
		void ReadParticles();
		/// @brief 把耦合器算出的 ε_g 写进求解器的体积分数数组（物理节点）
		void InjectVolumeFraction();
		/// @brief 由当前流场的 (ρ,p) 反算温度并刷新无量纲粘性系数数组（Sutherland，走求解器的 Gas）
		void UpdateViscosityField();
		/// @brief 评估所有粒子的受力（并缓存到 m_results）
		void EvaluateForces();
		/// @brief 颗粒运动后重建 ε 并注入（双向耦合每个时间步做一次）
		void RemapVolumeFraction();
		/// @brief 写 `result/particle_force.csv`（追加）
		void WriteForceHistory(int iter, double time);
		/// @brief 写 `result/particles_force_<iter>.dat`（Tecplot 点云，力矢量可视化用）
		void WriteForceFrame(int iter, double time);
		/// @brief 写 `result/coupling_report.csv`（机器可读的关键量）
		void WriteReport() const;

		// ---------------- 参数 ----------------
		std::string m_particle_file = "particles.csv";
		std::string m_drag_model_name = "gidaspow";
		std::string m_mapping_method = "subcell";
		int    m_mapping_sub = 2;
		int    m_smooth_passes = 0;
		double m_smooth_blend = 1.0;
		bool   m_volume_fraction_from_particles = true;
		bool   m_freeze_particles = true;
		double m_particle_density_fallback = 1000.0;
		double m_eps_min_clip = 1.0e-3;
		int    m_force_write_interval = 1;   ///< 0 表示跟随 iteration.write_interval
		bool   m_write_coupling_frames = true;

		// ---- 双向耦合（阶段 4）----
		std::string m_particle_motion_name = "fixed";
		std::string m_integrator_name = "exponential";
		bool   m_react_on_fluid = false;
		bool   m_react_work = true;             ///< 相间力是否对气相做功（S·u 进能量方程）
		std::string m_react_pg_mode = "auto";   ///< auto | on | off
		int    m_particle_substeps = 0;         ///< 0 = 自动
		double m_particle_max_step_param = 0.25;
		double m_body_acceleration[3] = { 0.0, 0.0, 0.0 };
		bool   m_remap_each_step = true;

		DEMCFDCoupler::Options m_options;
		DEMCFDCoupler::MotionOptions m_motion;

		// ---------------- 运行时状态 ----------------
		std::shared_ptr<NSFieldStructTwoPhase> m_field;
		std::shared_ptr<EulerTwoPhaseStructUniform> m_solver;
		std::shared_ptr<FlowSolverParamTwoPhase> m_para;
		std::shared_ptr<DataManagerNSTwoPhase> m_data_manager;

		DEMCFDCoupler m_coupler;
		std::vector<DEMCFDCoupler::ParticleResult> m_results;    ///< **步末**（输出用）
		std::vector<DEMCFDCoupler::ParticleResult> m_results_pre;///< **步初**（散射进流体的那一份）
		std::vector<double> m_eps_nodes;      ///< 注入求解器的 ε（含 ghost 的完整节点数组）
		std::vector<double> m_viscosity;      ///< 无量纲 μ（节点）
		std::vector<double> m_source_star;    ///< 相间反作用力源项（3*N，无量纲单位体积力）
		DEMCFDCoupler::ReactionSummary m_scatter_summary;

		double m_eps_min = 1.0, m_eps_max = 1.0;
		double m_alpha_max = 0.0;
		double m_partition_residual = 0.0;
		int    m_outside_count = 0;
		double m_dx_over_dp[3] = { 0.0, 0.0, 0.0 };
		bool   m_header_written = false;
		bool   m_two_way_active = false;      ///< 是否真的启用反作用/运动

		/// @brief 反应性/动量账本（全部无量纲）
		double m_source_impulse[3] = { 0.0, 0.0, 0.0 };   ///< Σ_cells src·V_cell·dt
		double m_particle_impulse[3] = { 0.0, 0.0, 0.0 }; ///< −Σ_p F*_p·dt
		double m_max_scatter_residual = 0.0;              ///< 两者相对偏差的最大值
		double m_weight_deficit = 0.0;
		int    m_last_substeps = 0;
		int    m_steps_advanced = 0;          ///< 实际推进颗粒的步数

		/// @brief 随时间累计的"最大 |F|"与"末态 ΣF_x"，用于报告与合理性判断
		double m_max_force = 0.0;
	double m_sum_abs_force = 0.0;   ///< 上一步 Σ_p |F|（L1 尺度，用于残差归一化）
		double m_max_particle_speed = 0.0;
		double m_last_sum_force[3] = { 0.0, 0.0, 0.0 };
		double m_last_time = 0.0;
	double m_last_time_for_dt = 0.0; ///< 上一个输出时刻（用来取 current_time 的增量作为真实 dt）
		int    m_last_iter = 0;
		int    m_react_pg_resolved = 0;       ///< 实际是否反作用了压力梯度力
	};
}
