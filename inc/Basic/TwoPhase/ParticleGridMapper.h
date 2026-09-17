/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file ParticleGridMapper.h
 * \brief Particle-to-grid mapping on a uniform Cartesian grid: solid volume fraction,
 *        momentum-conserving force scattering, and its dual interpolation.
 * \author Chen Jie.
 *
 *  阶段 2 的核心基础设施（见 docs/TWO_PHASE_STAGE1_MAPPING.md）
 *  ---------------------------------------------------------
 *  体积分数法（ε-weighted Euler）要求把 DEM 粒子的体积**精确地、尽量光滑地**
 *  分配到流体的均匀笛卡尔网格上。本类只做几何映射，不含任何物理，因此可以被
 *  单元测试到机器精度：
 *
 *    ComputeSolidFraction : 固相体积分数 α_s = V_solid/V_cell（ε_g = 1 − α_s）
 *    ScatterParticleVector: 把粒子量（如相间力）散列到网格，**Σ_cell = Σ_particle**
 *    GatherToParticle     : 把网格量插值到粒子位置，与 Scatter 用同一套权重 ⇒ 严格对偶
 *
 *  两种分配方法（都由构造保证**体积精确守恒**）：
 *    SubCell(n) 把粒子切成 n×n×n 个小立方体，每个小立方体把 V_p/n³ 整体记到包含其
 *               中心的那个单元。n=2 即 CFD-DEM 经典的"八分体/粒子中心法"，
 *               每个粒子最多触及 8 个单元，最便宜。n 越大越接近精确。
 *    Exact      球–单元交集体积，按 x 方向断点分段做高斯积分、内层用闭式
 *               "矩形∩圆盘"面积。每个粒子最多触及 27 个单元，精度接近机器零，
 *               代价约两个量级，适合做参考解与中小规模算例。
 *
 *  关于"单元"的约定（必须与 CFD 侧一致）
 *  --------------------------------------
 *  CFD 求解器把物理量存在**网格节点**上，残差是 dU_i/dt = -(F_{i+1/2}-F_{i-1/2})/Δx，
 *  即每个节点承担一个体积为 dx·dy·dz 的**控制体**，单元中心就是节点坐标。
 *  因此本映射器的单元中心 = CFD 节点坐标；整个计算域是
 *  [x0-dx/2, x0+(ni-1/2)dx] × ...，比节点跨度多出半个单元。DEM 侧的计算盒
 *  应当覆盖这个范围（阶段 3 的耦合驱动负责校验）。
 *
 *  守恒性为什么是精确的
 *  --------------------
 *  - SubCell：V_p 被切成 n³ 份、每份整体落到唯一一个单元，既不重也不漏；
 *  - Exact  ：球被所有相交单元切成互不重叠的完备分划，Σ 交集体积 ≡ V_p；
 *  - 光顺   ：通量型 Laplacian（面通量成对加减同一个数），Σα 逐位不变。
 */

#pragma once
#include "BasicType.h"

#include <cmath>

namespace zaran
{
	/// @brief 粒子 → 均匀笛卡尔网格的几何映射。
	class ParticleGridMapper
	{
	public:
		// ------------------------------------------------------------
		// 分配方法
		// ------------------------------------------------------------
		enum class Method
		{
			/// @brief n×n×n 子立方体，各占 V_p/n³，落到包含其中心的单元（n=2 即经典八分体法）
			SubCell,
			/// @brief 球–单元交集体积（按断点分段高斯积分 + 闭式矩形∩圆盘面积）
			Exact,
		};

		/// @brief 粒子（几何量）
		struct Particle
		{
			double x = 0.0, y = 0.0, z = 0.0;  ///< 质心坐标
			double radius = 0.0;               ///< 半径
		};

		/// @brief 网格描述（单元中心 = CFD 节点坐标）
		struct GridSpec
		{
			int ni = 1, nj = 1, nk = 1;
			/// @brief 空间维数（2 时把粒子当圆盘处理，单元面积 dx·dy）
			int dim = 3;
			double dx = 1.0, dy = 1.0, dz = 1.0;
			/// @brief 单元 (0,0,0) 的中心坐标
			double x0 = 0.0, y0 = 0.0, z0 = 0.0;

			/// @brief 单元体积（2 维时厚度取 1）
			double CellVolume() const
			{
				return (dim >= 3) ? (dx * dy * dz) : (dx * dy);
			}
			/// @brief 物理单元总数
			long long CellCount() const
			{
				return static_cast<long long>(ni) * nj * nk;
			}
		};

		/// @brief 一次映射的统计量（全部用于自检）
		struct Stats
		{
			double assigned_volume = 0.0;   ///< Σ α·V_cell（实际落到网格上的固体体积）
			double particle_volume = 0.0;   ///< Σ 4πr³/3（粒子总体积）
			double alpha_min = 0.0;
			double alpha_max = 0.0;
			double max_jump = 0.0;          ///< 相邻单元 |Δα| 的最大值（光滑性指标）
			long long cells_touched = 0;    ///< 被写过的单元数（去重前为写入次数）
			long long particle_cell_pairs = 0; ///< Σ 每个粒子触及的单元数
			double smoothing_shift = 0.0;   ///< 光顺引起的 Σα 相对变化（应 ~1e-16）
		};

	public:
		ParticleGridMapper() = default;
		explicit ParticleGridMapper(const GridSpec& grid) { Init(grid); }
		~ParticleGridMapper() = default;

		/// @brief 设定网格（会校验步长与维数）
		void Init(const GridSpec& grid);

		/// @brief 选择分配方法（默认 SubCell）
		void SetMethod(Method m) { m_method = m; }
		Method GetMethod() const { return m_method; }
		/// @brief SubCell 的每轴细分数（默认 2 = 经典八分体法）
		void SetSubdivisions(int n);
		int GetSubdivisions() const { return m_sub; }

		/// @brief 保守光顺：passes 次通量型 Laplacian 迭代，blend ∈ [0,1] 为混合权重
		/// @details 先算 Δα = Σ_邻居 (α_j - α_i)/(2·dim)，再 α += λ·Δα，面通量成对加减
		///          同一个数 ⇒ Σα 精确不变。最终结果取 (1-blend)·α_raw + blend·α_smooth
		///          —— 两者总和相同，因此混合仍然守恒，且不会引入负值（若原始场非负）。
		void SetSmoothing(int passes, double blend = 1.0);
		int GetSmoothingPasses() const { return m_smooth_passes; }
		double GetSmoothingBlend() const { return m_smooth_blend; }

		/// @brief 计算固相体积分数场
		/// @param particles  粒子列表
		/// @param alpha      输出，长度 = ni*nj*nk，索引 i + ni*(j + nj*k)
		/// @return 实际分配到的固体体积 Σα·V_cell（应等于 ΣV_p）
		double ComputeSolidFraction(const dynamic_array<Particle>& particles,
			dynamic_array<double>& alpha);

		/// @brief 把粒子上的矢量（如相间力）散列到网格
		/// @details 权重与 ComputeSolidFraction 完全一致 ⇒ Σ_cell F_cell = Σ_particle F_p
		/// @param particle_values 长度 = 3*particles.size()，按 (fx,fy,fz) 排列
		/// @param cell_values     输出，长度 = 3*ni*nj*nk，索引 (i + ni*(j + nj*k))*3 + c
		void ScatterParticleVector(const dynamic_array<Particle>& particles,
			const dynamic_array<double>& particle_values,
			dynamic_array<double>& cell_values) const;

		/// @brief 把网格上的矢量插值到粒子位置（Scatter 的对偶算子）
		/// @details 常数场插值后仍为该常数（Σ w = 1），这是权重归一化的直接检验
		/// @param cell_values     长度 = 3*ni*nj*nk
		/// @param particle_values 输出，长度 = 3*particles.size()
		void GatherToParticle(const dynamic_array<Particle>& particles,
			const dynamic_array<double>& cell_values,
			dynamic_array<double>& particle_values) const;

		/// @brief 上一次 ComputeSolidFraction 的统计
		const Stats& GetStats() const { return m_stats; }

		/// @brief 单元线性索引
		long long Idx(int i, int j, int k) const
		{
			return static_cast<long long>(i) + m_grid.ni * (static_cast<long long>(j) + m_grid.nj * static_cast<long long>(k));
		}

		// ------------------------------------------------------------
		// 供验证脚本调用的几何工具（静态，纯函数）
		// ------------------------------------------------------------
		/// @brief 球（半径 R，球心在原点）与轴对齐盒 [x0,x1]×[y0,y1]×[z0,z1] 的交集体积
		static double SphereBoxVolume(double R,
			double x0, double x1, double y0, double y1, double z0, double z1);

		/// @brief 圆盘（半径 ρ，圆心在原点）与轴对齐矩形 [y0,y1]×[z0,z1] 的相交面积
		static double DiskRectArea(double rho,
			double y0, double y1, double z0, double z1);

		/// @brief 正下界角区面积 { y ≥ u, z ≥ v, y²+z² ≤ ρ² }（u,v ≥ 0）
		static double CornerArea(double u, double v, double rho);

	private:
		/// @brief 按当前方法把一个粒子的单位体积分配到单元上（Σ w = 1）
		void BuildWeights(const Particle& p,
			dynamic_array<long long>& cells, dynamic_array<double>& weights) const;
		/// @brief SubCell 权重
		void BuildWeightsSubCell(const Particle& p,
			dynamic_array<long long>& cells, dynamic_array<double>& weights) const;
		/// @brief Exact 权重
		void BuildWeightsExact(const Particle& p,
			dynamic_array<long long>& cells, dynamic_array<double>& weights) const;
		/// @brief 通量型 Laplacian 光顺（原地，Σ 不变）
		void SmoothInPlace(dynamic_array<double>& alpha) const;

		bool InDomain(double x, double y, double z) const
		{
			return (x >= m_lo[0] && x <= m_hi[0]
				&& y >= m_lo[1] && y <= m_hi[1]
				&& z >= m_lo[2] && z <= m_hi[2]);
		}

		GridSpec m_grid;
		Method m_method = Method::SubCell;
		int m_sub = 2;
		int m_smooth_passes = 0;
		double m_smooth_blend = 0.0;

		/// @brief 单元的物理边界（含半厚度）
		double m_lo[3] = { 0.0, 0.0, 0.0 };
		double m_hi[3] = { 0.0, 0.0, 0.0 };

		Stats m_stats;

		// 复用缓冲，避免每次分配
		mutable dynamic_array<long long> m_buf_cells;
		mutable dynamic_array<double> m_buf_weights;
		mutable dynamic_array<double> m_buf_field;
	};
}
