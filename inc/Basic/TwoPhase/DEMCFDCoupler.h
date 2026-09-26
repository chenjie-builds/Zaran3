/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file DEMCFDCoupler.h
 * \brief One-way DEM-CFD coupling: particle -> grid voidage, grid -> particle gas state,
 *        and the inter-phase force on every particle.
 * \author Chen Jie.
 *
 *  阶段 3 的"相间数据接口"（见 docs/TWO_PHASE_STAGE3_DEM_CFD.md）
 *  -----------------------------------------------------------
 *  本类把两条链路显式地写在一起，是 CFD 侧与 DEM 侧**唯一**的接触面：
 *
 *      DEM → CFD      粒子几何 → ParticleGridMapper → 固相体积分数 α_s
 *                     → 气相体积分数 ε_g = 1 − α_s（CFD 物理节点）
 *
 *      CFD → DEM      网格上的 (ρ,u,v,w,p,μ) 与 ∇p → 与体积分配**同一套权重**插值
 *                     到粒子位置 → DragForceModel → 曳力 + 压力梯度力
 *
 *  单向耦合的三个要点
 *  ------------------
 *  1. **体积分数是"冻结"的几何量**：粒子不动 ⇒ ε_g 只算一次，之后每个时间步都用同一
 *     个数组。求解器侧由 `VolumeFractionType::External` 接收它（不是每步重算的解析场）。
 *  2. **力不反馈**：本类只输出作用在粒子上的力，唯一参与 CFD 的是 ε。因此
 *     `ε_g ≡ 1`（`volume_fraction_from_particles = false`）时整条 CFD 链路必须与
 *     单相 Euler 求解器**逐位相同** —— 这是"单向耦合没有扰动流体"的直接检验。
 *  3. **单位只在两个地方出现**：求解器内部是无量纲量，颗粒参数是 SI。本类用求解器
 *     自己的参考量（Dimensionless）做换算，换算关系写在 UnitScales 里，避免"另一处
 *     手写的参考速度"与求解器漂移。
 *
 *  单元约定（与阶段 2 一致，必须与 CFD 侧对齐）
 *  -------------------------------------------
 *  CFD 是**节点型有限体积**：物理量在节点上，节点 i 承担体积 dx·dy·dz 的控制体，
 *  即"单元中心 = 节点坐标"。映射器单元 c 与 CFD 物理节点一一对应：
 *
 *      CFD 物理节点 (is + c_i, js + c_j, ks + c_k)  ←→  映射器单元 (c_i, c_j, c_k)
 *
 *  其中 (is, js, ks) 是 CFD 不含 ghost 的计算范围起点（= ghost 层数）。
 *  CFD 的控制体总范围比节点跨度多出半个单元：
 *      [x0 − dx/2, x0 + (n_i − ½)·dx] × …
 *  粒子必须落在这个范围内，否则它的体积不会（完整地）落到网格上。
 */

#pragma once
#include "BasicType.h"
#include "DragForceModel.h"
#include "ParticleGridMapper.h"

#include <string>
#include <set>
#include <vector>

namespace zaran
{
	/// @brief 单向耦合的 DEM-CFD 相间数据接口。
	class DEMCFDCoupler
	{
	public:
		// ------------------------------------------------------------
		// 粒子（几何 + 动力学）
		// ------------------------------------------------------------
		/// @brief 参与耦合的粒子。单向耦合阶段速度恒为 0（保持静止），
		///        但接口保留 velocity，双向耦合时只需填上真实速度即可。
		struct Particle
		{
			index_type id = 0;                 ///< 全局编号（写进报告）
			int        group = 0;              ///< 分组（预留：按组给不同材料参数）
			double     x = 0.0, y = 0.0, z = 0.0; ///< 质心（SI, m）
			double     radius = 0.0;           ///< 半径（SI, m）
			double     density = 0.0;          ///< 密度（SI, kg/m³）
			double     velocity[3] = { 0.0, 0.0, 0.0 }; ///< 速度（SI, m/s）
		};

		// ------------------------------------------------------------
		// CFD 网格视图（只读）
		// ------------------------------------------------------------
		/// @brief CFD 的均匀结构网格视图。所有数组都是"节点型"线性索引 = Idx(i,j,k)。
		struct GridView
		{
			int    dim = 3;                    ///< 只支持 3（薄板算例仍按 3 维处理）
			int    ni = 0, nj = 0, nk = 0;      ///< **含 ghost** 的节点总数
			int    is = 0, ie = 0;               ///< 不含 ghost 的 i 范围
			int    js = 0, je = 0;
			int    ks = 0, ke = 0;
			double dx = 1.0, dy = 1.0, dz = 1.0;
			/// @brief 物理节点 (is, js, ks) 的坐标（= 第 0 个映射器单元的中心）
			double x0 = 0.0, y0 = 0.0, z0 = 0.0;

			/// @brief 节点线性索引（必须与 CFD 侧的 IdProxyStruct 一致）
			long long Idx(int i, int j, int k) const
			{
				return static_cast<long long>(i)
					+ static_cast<long long>(ni) * (static_cast<long long>(j)
						+ static_cast<long long>(nj) * static_cast<long long>(k));
			}
			int  NodeCount() const { return ni * nj * nk; }
			int  PhysicalNi() const { return ie - is + 1; }
			int  PhysicalNj() const { return je - js + 1; }
			int  PhysicalNk() const { return ke - ks + 1; }
			double CellVolume() const { return dx * dy * dz; }
			bool Valid() const
			{
				return dim == 3 && ni > 0 && nj > 0 && nk > 0
					&& ie >= is && je >= js && ke >= ks && dx > 0.0 && dy > 0.0 && dz > 0.0;
			}
		};

		// ------------------------------------------------------------
		// 单位换算（由求解器的 Dimensionless 给出）
		// ------------------------------------------------------------
		/// @brief 无量纲 → SI 的换算因子。求解器内部一律无量纲：
		///        ρ* = ρ/ρ_ref，u* = u/a_ref，p* = p/(ρ_ref a_ref²)，
		///        T* = T/T_ref，μ* = μ/(ρ_ref a_ref L_ref)，x* = x/L_ref。
		///        力：F* = F/(ρ_ref a_ref² L_ref²)。
		struct UnitScales
		{
			double ref_length = 1.0;       ///< L_ref (m)
			double ref_density = 1.0;      ///< ρ_ref (kg/m³)
			double ref_velocity = 1.0;     ///< a_ref = √(γ R T_ref) (m/s)
			double ref_temperature = 1.0;  ///< T_ref (K)
			double ref_gamma = 1.4;        ///< γ（用于把 p*/ρ* 还原成 T* = γp*/ρ*）

			double RefPressure() const { return ref_density * ref_velocity * ref_velocity; }
			double RefViscosity() const { return ref_density * ref_velocity * ref_length; }
			double RefForce() const
			{
				return ref_density * ref_velocity * ref_velocity * ref_length * ref_length;
			}
		};

		// ------------------------------------------------------------
		// 建模选项
		// ------------------------------------------------------------
		struct Options
		{
			DragForceModel::Model model = DragForceModel::Model::Gidaspow;
			DragForceModel::Options drag;

			ParticleGridMapper::Method mapping = ParticleGridMapper::Method::SubCell;
			int    mapping_sub = 2;            ///< SubCell 每轴细分数
			int    smooth_passes = 0;          ///< 保守光顺次数（0 = 不启用）
			double smooth_blend = 1.0;

			/// @brief ε_g 是否由粒子给出。false ⇒ 恒为 1（严格"无体积排斥"的单向耦合）
			bool   volume_fraction_from_particles = true;
			/// @brief ε_g 的下限裁剪（写进 CFD 前）
			double eps_min_clip = 1.0e-3;
		};

		// ------------------------------------------------------------
		// CFD → 粒子 插值出来的气态样本（全部 SI，逐列写进 particle_force.csv）
		// ------------------------------------------------------------
		struct GasSample
		{
			double density = 0.0;              ///< ρ_g (kg/m³)
			double velocity[3] = { 0.0, 0.0, 0.0 }; ///< u_g (m/s)
			double pressure = 0.0;             ///< p_g (Pa)
			double temperature = 0.0;          ///< T_g (K)
			double viscosity = 0.0;            ///< μ_g (Pa·s)
			double eps_gas = 1.0;              ///< 插值到粒子处的空隙率
			double grad_p[3] = { 0.0, 0.0, 0.0 };   ///< ∇p (Pa/m)
			/// @brief Σw（应为 1；< 1 说明粒子有一部分在计算域外）
			double weight_sum = 0.0;
			long long cells_touched = 0;       ///< 该粒子触及的单元数
		};

		/// @brief 逐粒子的耦合结果
		struct ParticleResult
		{
			GasSample            gas;          ///< 插值得到的当地气态（SI）
			DragForceModel::Result force;      ///< 受力（曳力 + 压力梯度力，单位：N）
		};

		/// @brief CFD 的原始变量视图（全部为**无量纲**节点数组，长度 = ni*nj*nk）
		struct PrimitiveView
		{
			const double* density = nullptr;                        ///< primitive_0
			const double* velocity[3] = { nullptr, nullptr, nullptr }; ///< primitive_1..3
			const double* pressure = nullptr;                       ///< primitive_4
			/// @brief 粘性系数（无量纲）。为 nullptr 时由本类用 γ、p、ρ 反算温度后
			///        交给调用方提供的 Sutherland 钩子；为了"与实际求解器完全一致"，
			///        正常路径都是直接传入求解器算好的数组。
			const double* viscosity = nullptr;

			bool Valid() const
			{
				return density != nullptr && pressure != nullptr
					&& velocity[0] != nullptr && velocity[1] != nullptr
					&& velocity[2] != nullptr;
			}
		};

		// ------------------------------------------------------------
		// 双向耦合：颗粒运动 + 反作用力
		// ------------------------------------------------------------
		/// @brief 颗粒运动模式
		enum class MotionMode
		{
			/// @brief 保持静止（阶段 3 的单向耦合行为）
			Fixed,
			/// @brief 按输入速度匀速运动（"颗粒在静止流场中运动"的算例：把颗粒运动
			///        从求解里拿出来，单独考察它对流场的反作用）
			Prescribed,
			/// @brief 按 m·du/dt = F_drag + F_∇p + m·a_body 积分（"冲击波推动颗粒"）
			Dynamic,
		};

		/// @brief 颗粒时间积分器
		enum class Integrator
		{
			/// @brief 把曳力线性化在步初状态上、对线性 ODE **解析积分**：
			///        `u ← u_∞ + (u − u_∞)·exp(−k·Δt)`（k = V_p·β/ε_s/m_p，u_∞ = u_g + a/k）。
			///        对 Stokes 曳力 + 冻结气态是**精确解**，且无条件稳定（推荐）
			Exponential,
			/// @brief 显式中点法（RK2）。对刚性曳力需要子步
			Rk2,
		};

		/// @brief 颗粒运动选项
		struct MotionOptions
		{
			MotionMode mode = MotionMode::Fixed;
			Integrator integrator = Integrator::Exponential;
			/// @brief 固定子步数（<= 0 表示按 `k·Δt ≤ max_step_parameter` 自动确定）
			int    substeps = 0;
			double max_step_parameter = 0.25;
			/// @brief 体积力/质量（SI, m/s²）—— 例如重力 (0,-9.8,0)。
			/// @details 这是"外部"施加在颗粒上的力，不参与流体的动量账本
			///          （反作用力只来自流体→颗粒的相间力）。
			double body_acceleration[3] = { 0.0, 0.0, 0.0 };
		};

		/// @brief 反作用力散列的统计（用于核对牛顿第三定律）
		struct ReactionSummary
		{
			/// @brief Σ_particle F*（无量纲，流体 → 颗粒的力）
			double particle_force_dimensionless[3] = { 0.0, 0.0, 0.0 };
			/// @brief Σ_cells src*·V*_cell（无量纲，散列到网格上的总力）
			double scattered_force_dimensionless[3] = { 0.0, 0.0, 0.0 };
			/// @brief Σ_p Σ_d |F*_d|（L1 尺度）——**散列残差的分母必须用它**：
			///        粒子布置对称时 ΣF* 本身可以恒等于 0，拿它当分母会得到 0/0 的
			///        虚假大残差（实测：对称布置下 ΣF* ≈ 1e-22，残差虚报 277%）。
			double particle_force_l1 = 0.0;
			/// @brief Σw 偏离 1 的总量（>0 说明有粒子跨出计算域 ⇒ 反作用被少算）
			double weight_deficit = 0.0;
		};

	public:
		DEMCFDCoupler() = default;
		~DEMCFDCoupler() = default;

		/// @brief 设定 CFD 网格视图、单位换算与建模选项；同时按 CFD 物理范围构造映射器
		void Init(const GridView& grid, const UnitScales& scales, const Options& options);

		/// @brief 设定粒子（会校验：半径、密度为正，且质心落在计算域控制体内）
		/// @return 质心落在计算域外的粒子数
		int SetParticles(const std::vector<Particle>& particles);

		/// @brief 计算 CFD 物理节点上的气相体积分数 ε_g
		/// @details 写入 `eps_nodes`（长度 = ni*nj*nk，含 ghost 的完整节点数组）的
		///          **物理节点**；ghost 部分保持原值，由求解器做零梯度外推。
		///          ε_g = 1 − α_s，α_s 由 ParticleGridMapper 给出（体积精确配分）。
		/// @param[out] eps_min / eps_max 实际写入的 ε 范围
		/// @param[in]  verbose 是否打日志（每个时间步重映射时传 false，避免刷屏）
		void ComputeGasVolumeFraction(std::vector<double>& eps_nodes,
			double& eps_min, double& eps_max, bool verbose = true) const;

		/// @brief 评估全部粒子的受力
		/// @param prim   无量纲原始变量视图
		/// @param outputs 输出，长度 = 粒子数；同时给出插值出的气态样本与受力
		void EvaluateParticleForces(const PrimitiveView& prim,
			std::vector<ParticleResult>& outputs) const;

		/// @brief 相间反作用力散列（**双向耦合的核心**）
		/// @details 用与体积分配/插值**完全相同**的权重，把颗粒受到的流体作用力取负后
		///          摊回网格节点。构造保证（Σw = 1 时逐位成立）：
		///              Σ_cells src*·V*_cell = −Σ_particle F*
		///          输出是**无量纲的单位体积力**（求解器内部单位），可直接塞给
		///          `EulerTwoPhaseStructUniform::SetInterphaseMomentumSource()`。
		/// @param include_pressure_gradient_force 是否连 −V_p∇p 一起反作用。
		///        ⚠ 只有气相动量方程取**完全守恒**的 ∇(εp) 形式
		///          （`two_phase.porosity_gradient_force = 0`）时才该为 true；
		///          取 ε∇p 形式（=1，默认）时为 false —— 那种情况下 −V_p∇p 已经是
		///          气相方程里 −ε_g∇p 的对偶项，再加一次就是双计。
		void ScatterReaction(const std::vector<ParticleResult>& results,
			bool include_pressure_gradient_force,
			std::vector<double>& src_dimensionless, ReactionSummary& summary) const;

		/// @brief 把颗粒推进一个时间步（`mode != Dynamic` 时什么都不做）
		/// @return 实际使用的子步数
		int AdvanceParticles(double dt, const MotionOptions& options,
			const std::vector<ParticleResult>& results);

		/// @brief 计算物理节点的压强梯度（中心差分，物理范围边界退化为单侧差分）
		/// @details `EvaluateParticleForces` 在需要压力梯度力时会自己调用它；
		///          公开出来是为了让校验脚本/单元测试能单独核对这一项。
		void ComputePressureGradient(const double* pressure_dimensionless,
			std::vector<double>& grad_p_dimensionless) const;

		// ------------------------------------------------------------
		// 访问
		// ------------------------------------------------------------
		const std::vector<Particle>& GetParticles() const { return m_particles; }
		/// @brief 可写粒子列表（双向耦合需要更新速度与位置）
		std::vector<Particle>& MutableParticles() { return m_particles; }
		const GridView& GetGrid() const { return m_grid; }
		const Options& GetOptions() const { return m_options; }
		const UnitScales& GetScales() const { return m_scales; }
		const ParticleGridMapper& GetMapper() const { return m_mapper; }
		/// @brief 颗粒质量 m_p = ρ_p·V_p (SI)
		double ParticleMass(size_t index) const;
		/// @brief 无量纲力的换算：F* = F_SI/F_ref
		double ForceToDimensionless(double force_si) const
		{
			return force_si / m_scales.RefForce();
		}
		/// @brief Σ_cells src·V*_cell（核对牛顿第三定律用）
		void SumSourceMomentum(const std::vector<double>& src_dimensionless,
			double out[3]) const;
		/// @brief 粒子总体积 ΣV_p（SI）
		double GetParticleVolume() const { return m_particle_volume; }
		/// @brief 代表性粒径 d_p（取最小半径 ⇒ 对应最苛刻的 Δx/d_p）
		double GetRepDiameter() const { return m_rep_diameter; }
		/// @brief 计算域控制体范围（比节点跨度多出半格）
		void GetControlVolumeBox(double lo[3], double hi[3]) const;
		/// @brief Δx/d_p、Δy/d_p、Δz/d_p（体积分数法的适用性自检）
		void GetResolutionRatios(double ratios[3]) const;

	private:
		/// @brief 把粒子质心映射到映射器单元坐标（用于落域判定）
		bool ParticleInside(const Particle& p) const;

		GridView m_grid;
		UnitScales m_scales;
		Options m_options;
		/// @brief mutable：ComputeSolidFraction 会写内部统计量，而 ComputeGasVolumeFraction
		///        是 const 方法（对外语义上"只是查询"）
		mutable ParticleGridMapper m_mapper;

		std::vector<Particle> m_particles;
		double m_particle_volume = 0.0;   ///< ΣV_p（SI）
		double m_rep_diameter = 1.0;      ///< d_p（SI）
		int m_outside_count = 0;

		// 复用缓冲
		mutable std::vector<long long> m_buf_cells;
		mutable std::vector<double> m_buf_weights;
		mutable std::vector<double> m_alpha;      ///< 映射器给出的固相体积分数（映射器网格）
		mutable std::vector<double> m_grad_p;     ///< 无量纲 ∇p，长度 3*N
		/// @brief 已经告警过"跑出计算域"的粒子 id（避免每步刷屏）
		mutable std::set<index_type> m_warned_outside;
	};
}
