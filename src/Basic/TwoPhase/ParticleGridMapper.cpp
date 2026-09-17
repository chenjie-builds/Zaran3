/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file ParticleGridMapper.cpp
 * \brief Particle-to-grid mapping on a uniform Cartesian grid.
 * \author Chen Jie.
 */

#include "ParticleGridMapper.h"
#include "Log.h"
#include "ZaranError.h"

#include <algorithm>

namespace zaran
{
	namespace
	{
		constexpr double kPi = 3.14159265358979323846;

		/// @brief Gauss-Legendre 节点与权重（Newton 迭代求 Legendre 多项式零点，只算一次）
		struct GaussLegendreTable
		{
			int order = 0;
			dynamic_array<double> nodes;
			dynamic_array<double> weights;

			void Build(int K)
			{
				if (order == K)
				{
					return;
				}
				order = K;
				nodes.assign(K, 0.0);
				weights.assign(K, 0.0);
				const int m = (K + 1) / 2;
				for (int i = 0; i < m; ++i)
				{
					// 初始近似（Abramowitz & Stegun 22.16.5）
					double t = std::cos(kPi * (i + 0.75) / (K + 0.5));
					double p_k = 1.0, p_km1 = 1.0;
					for (int iter = 0; iter < 200; ++iter)
					{
						// P_K(t) 与 P_{K-1}(t)
						p_km1 = 1.0;
						p_k = t;
						for (int n = 2; n <= K; ++n)
						{
							const double p_n = ((2.0 * n - 1.0) * t * p_k - (n - 1.0) * p_km1) / n;
							p_km1 = p_k;
							p_k = p_n;
						}
						const double dp = K * (t * p_k - p_km1) / (t * t - 1.0);
						const double dt = p_k / dp;
						t -= dt;
						if (std::fabs(dt) < 1.0e-16)
						{
							break;
						}
					}
					// 最后再算一次导数用于权重
					p_km1 = 1.0;
					p_k = t;
					for (int n = 2; n <= K; ++n)
					{
						const double p_n = ((2.0 * n - 1.0) * t * p_k - (n - 1.0) * p_km1) / n;
						p_km1 = p_k;
						p_k = p_n;
					}
					const double dp = K * (t * p_k - p_km1) / (t * t - 1.0);
					const double w = 2.0 / ((1.0 - t * t) * dp * dp);
					nodes[i] = -t;
					nodes[K - 1 - i] = t;
					weights[i] = w;
					weights[K - 1 - i] = w;
				}
			}
		};

		const GaussLegendreTable& GetGaussTable(int order)
		{
			static GaussLegendreTable table;
			table.Build(order);
			return table;
		}
	}

	// ==================================================================
	// 几何工具：矩形∩圆盘、球∩长方体
	// ==================================================================
	/// @brief 正下界角区面积 { y ≥ u, z ≥ v, y² + z² ≤ ρ² }（要求 u,v ≥ 0）
	/// @details ∫_u^Y (√(ρ²-y²) - v) dy，Y = √(ρ²-v²)，有闭式解。
	double ParticleGridMapper::CornerArea(double u, double v, double rho)
	{
		if (!(rho > 0.0))
		{
			return 0.0;
		}
		if (u < 0.0) u = 0.0;
		if (v < 0.0) v = 0.0;
		if (u >= rho || v >= rho)
		{
			return 0.0;
		}
		if (u * u + v * v >= rho * rho)
		{
			return 0.0;
		}
		const double Y = std::sqrt(std::max(0.0, rho * rho - v * v));
		const double uu = std::min(u, Y);
		const double su = std::sqrt(std::max(0.0, rho * rho - uu * uu));
		const double value = 0.5 * (Y * v) + 0.5 * rho * rho * std::acos(v / rho)
			- 0.5 * (uu * su) - 0.5 * rho * rho * std::asin(uu / rho)
			- v * (Y - uu);
		return std::max(0.0, value);
	}

	/// @brief 圆盘（半径 ρ，圆心在原点）与轴对齐矩形的相交面积
	/// @details 先按坐标轴把矩形切成 ≤4 个象限子矩形，再对每个子矩形用
	///          "1_{a≤x≤b} = 1_{x≥a} - 1_{x≥b}" 的四项容斥，落到 CornerArea 上。
	double ParticleGridMapper::DiskRectArea(double rho,
		double y0, double y1, double z0, double z1)
	{
		if (!(rho > 0.0))
		{
			return 0.0;
		}
		// 区间按 0 切分成若干"同号"子区间，并归一成 [a,b] ⊂ [0,∞)
		// 注意：最多 2 个子区间 × 2 个端点 = 4 个槽位，数组必须开 4（开 3 会栈越界）
		double ys[4], zs[4];
		int ny = 0, nz = 0;
		{
			const double cuts[3] = { y0, 0.0, y1 };
			// 只保留落在 [y0,y1] 内的切点，并排序去重
			double tmp[3];
			int n = 0;
			for (int i = 0; i < 3; ++i)
			{
				const double c = cuts[i];
				if (c >= y0 && c <= y1)
				{
					bool dup = false;
					for (int j = 0; j < n; ++j)
					{
						if (tmp[j] == c) { dup = true; break; }
					}
					if (!dup) tmp[n++] = c;
				}
			}
			std::sort(tmp, tmp + n);
			for (int i = 0; i + 1 < n; ++i)
			{
				ys[2 * ny] = std::min(std::fabs(tmp[i]), std::fabs(tmp[i + 1]));
				ys[2 * ny + 1] = std::max(std::fabs(tmp[i]), std::fabs(tmp[i + 1]));
				++ny;
			}
		}
		{
			const double cuts[3] = { z0, 0.0, z1 };
			double tmp[3];
			int n = 0;
			for (int i = 0; i < 3; ++i)
			{
				const double c = cuts[i];
				if (c >= z0 && c <= z1)
				{
					bool dup = false;
					for (int j = 0; j < n; ++j)
					{
						if (tmp[j] == c) { dup = true; break; }
					}
					if (!dup) tmp[n++] = c;
				}
			}
			std::sort(tmp, tmp + n);
			for (int i = 0; i + 1 < n; ++i)
			{
				zs[2 * nz] = std::min(std::fabs(tmp[i]), std::fabs(tmp[i + 1]));
				zs[2 * nz + 1] = std::max(std::fabs(tmp[i]), std::fabs(tmp[i + 1]));
				++nz;
			}
		}
		double total = 0.0;
		for (int a = 0; a < ny; ++a)
		{
			const double ya = ys[2 * a], yb = ys[2 * a + 1];
			for (int b = 0; b < nz; ++b)
			{
				const double za = zs[2 * b], zb = zs[2 * b + 1];
				total += CornerArea(ya, za, rho) - CornerArea(yb, za, rho)
					- CornerArea(ya, zb, rho) + CornerArea(yb, zb, rho);
			}
		}
		return std::max(0.0, total);
	}

	/// @brief 球（半径 R，球心在原点）与轴对齐长方体的交集体积
	/// @details 沿 x 做 1 维积分：V = ∫ A_rect∩disk(ρ(x)) dx，ρ(x) = √(R²-x²)。
	///          为消除端点处 √ 造成的半幂奇性，换成 x = R·sinθ，于是
	///          V = ∫ R·cosθ · A(R·cosθ) dθ。
	///
	///          **断点必须取全**，否则被积函数在段内仍有折点、高斯收敛很慢：
	///            ρ = |y_i|、|z_j|          —— 圆盘扫过矩形的四条边
	///            ρ = √(y_i² + z_j²)        —— 圆盘扫过矩形的四个角点
	///          对应的 θ = ±acos(ρ/R)。实测（纯配分恒等式 Σ V = 4πR³/3）：
	///            只取边界 ρ 且用 asin 断点：Δx/R=0.5 时残差 2.4e-8
	///            取全（边+角）且用 acos 断点、24 点高斯：残差 ~1e-13
	///          另外单元"整块在球内/球外"时直接返回解析值，既精确又省时间。
	double ParticleGridMapper::SphereBoxVolume(double R,
		double x0, double x1, double y0, double y1, double z0, double z1)
	{
		if (!(R > 0.0))
		{
			return 0.0;
		}
		// 整块在球内 / 球外的精确早退
		const double fx = std::max(std::fabs(x0), std::fabs(x1));
		const double fy = std::max(std::fabs(y0), std::fabs(y1));
		const double fz = std::max(std::fabs(z0), std::fabs(z1));
		if (fx * fx + fy * fy + fz * fz <= R * R)
		{
			return (x1 - x0) * (y1 - y0) * (z1 - z0);   // 整块在球内
		}
		auto nearest_abs = [](double a, double b)
		{
			return (a * b <= 0.0) ? 0.0 : std::min(std::fabs(a), std::fabs(b));
		};
		const double nx = nearest_abs(x0, x1);
		const double ny = nearest_abs(y0, y1);
		const double nz = nearest_abs(z0, z1);
		if (nx * nx + ny * ny + nz * nz >= R * R)
		{
			return 0.0;                                 // 整块在球外
		}

		// x 方向裁剪到球的范围，并转成 θ 区间
		const double xa = std::max(x0, -R);
		const double xb = std::min(x1, R);
		if (!(xb > xa))
		{
			return 0.0;
		}
		const double th_a = std::asin(std::max(-1.0, std::min(1.0, xa / R)));
		const double th_b = std::asin(std::max(-1.0, std::min(1.0, xb / R)));
		if (!(th_b > th_a))
		{
			return 0.0;
		}

		// θ 断点：ρ(θ) = R·cosθ 扫过四条边与四个角点
		const double ys[2] = { y0, y1 };
		const double zs[2] = { z0, z1 };
		double bp[32];
		int nbp = 0;
		bp[nbp++] = th_a;
		bp[nbp++] = th_b;
		auto add_rho_break = [&](double rho)
		{
			if (rho > 0.0 && rho < R)
			{
				const double t = std::acos(rho / R);
				if (t > th_a && t < th_b && nbp < 30) { bp[nbp++] = t; }
				if (-t > th_a && -t < th_b && nbp < 31) { bp[nbp++] = -t; }
			}
		};
		for (int i = 0; i < 2; ++i)
		{
			add_rho_break(std::fabs(ys[i]));
			add_rho_break(std::fabs(zs[i]));
			for (int j = 0; j < 2; ++j)
			{
				add_rho_break(std::sqrt(ys[i] * ys[i] + zs[j] * zs[j]));
			}
		}
		std::sort(bp, bp + nbp);

		constexpr int kOrder = 24;
		const GaussLegendreTable& gl = GetGaussTable(kOrder);
		double total = 0.0;
		for (int s = 0; s + 1 < nbp; ++s)
		{
			const double a = bp[s], b = bp[s + 1];
			if (!(b > a))
			{
				continue;
			}
			const double half = 0.5 * (b - a);
			const double mid = 0.5 * (a + b);
			double acc = 0.0;
			for (int q = 0; q < kOrder; ++q)
			{
				const double theta = mid + half * gl.nodes[q];
				const double ct = std::cos(theta);
				const double rho = R * ct;      // ≥ 0（θ ∈ [-π/2, π/2]）
				acc += gl.weights[q] * ct * DiskRectArea(rho, y0, y1, z0, z1);
			}
			total += R * half * acc;
		}
		return std::max(0.0, total);
	}

	// ==================================================================
	// 初始化
	// ==================================================================
	void ParticleGridMapper::Init(const GridSpec& grid)
	{
		m_grid = grid;
		if (m_grid.dim != 2 && m_grid.dim != 3)
		{
			throw ZaranError("ParticleGridMapper: only 2-D and 3-D grids are supported");
		}
		if (!(m_grid.dx > 0.0) || !(m_grid.dy > 0.0) || (m_grid.dim >= 3 && !(m_grid.dz > 0.0)))
		{
			throw ZaranError("ParticleGridMapper: non-positive grid spacing");
		}
		if (m_grid.ni < 1 || m_grid.nj < 1 || m_grid.nk < 1)
		{
			throw ZaranError("ParticleGridMapper: invalid grid size");
		}
		// 单元物理边界（节点型控制体：中心 ± 半间距）
		m_lo[0] = m_grid.x0 - 0.5 * m_grid.dx;
		m_hi[0] = m_grid.x0 + (m_grid.ni - 0.5) * m_grid.dx;
		m_lo[1] = m_grid.y0 - 0.5 * m_grid.dy;
		m_hi[1] = m_grid.y0 + (m_grid.nj - 0.5) * m_grid.dy;
		if (m_grid.dim >= 3)
		{
			m_lo[2] = m_grid.z0 - 0.5 * m_grid.dz;
			m_hi[2] = m_grid.z0 + (m_grid.nk - 0.5) * m_grid.dz;
		}
		else
		{
			m_lo[2] = -0.5;
			m_hi[2] = 0.5;
		}
		Log::info("ParticleGridMapper: dim={}, cells=({},{},{}), spacing=({:E},{:E},{:E}), "
			"cell_center0=({:E},{:E},{:E}), cell volume={:E} m^3",
			m_grid.dim, m_grid.ni, m_grid.nj, m_grid.nk,
			m_grid.dx, m_grid.dy, m_grid.dz, m_grid.x0, m_grid.y0, m_grid.z0,
			m_grid.CellVolume());
	}

	void ParticleGridMapper::SetSubdivisions(int n)
	{
		m_sub = std::max(1, n);
	}

	void ParticleGridMapper::SetSmoothing(int passes, double blend)
	{
		m_smooth_passes = std::max(0, passes);
		m_smooth_blend = std::max(0.0, std::min(1.0, blend));
		if (m_smooth_passes == 0)
		{
			m_smooth_blend = 0.0;
		}
	}

	// ==================================================================
	// 权重构造
	// ==================================================================
	void ParticleGridMapper::BuildWeights(const Particle& p,
		dynamic_array<long long>& cells, dynamic_array<double>& weights) const
	{
		cells.clear();
		weights.clear();
		if (m_method == Method::Exact)
		{
			BuildWeightsExact(p, cells, weights);
		}
		else
		{
			BuildWeightsSubCell(p, cells, weights);
		}
	}

	/// @brief 把粒子切成 n×n×n 个小立方体，每个小立方体把 V_p/n³ 整体记到包含其中心的单元。
	/// @details 这就是 CFD-DEM 里的 "particle-centre method"（n=2 时即经典八分体法）。
	///          每个小立方体只落一个单元 ⇒ 精确配分（既不重也不漏）。
	void ParticleGridMapper::BuildWeightsSubCell(const Particle& p,
		dynamic_array<long long>& cells, dynamic_array<double>& weights) const
	{
		const int n = m_sub;
		const int nz = (m_grid.dim >= 3) ? n : 1;
		const double count = static_cast<double>(n) * n * (m_grid.dim >= 3 ? n : 1);
		const double w = 1.0 / count;
		for (int a = 0; a < n; ++a)
		{
			const double ox = p.radius * (2.0 * (a + 0.5) / n - 1.0);
			const double x = p.x + ox;
			const int i = static_cast<int>(std::floor((x - m_grid.x0) / m_grid.dx + 0.5));
			if (i < 0 || i >= m_grid.ni)
			{
				continue;
			}
			for (int b = 0; b < n; ++b)
			{
				const double oy = p.radius * (2.0 * (b + 0.5) / n - 1.0);
				const double y = p.y + oy;
				const int j = static_cast<int>(std::floor((y - m_grid.y0) / m_grid.dy + 0.5));
				if (j < 0 || j >= m_grid.nj)
				{
					continue;
				}
				for (int c = 0; c < nz; ++c)
				{
					int k = 0;
					if (m_grid.dim >= 3)
					{
						const double oz = p.radius * (2.0 * (c + 0.5) / n - 1.0);
						const double z = p.z + oz;
						k = static_cast<int>(std::floor((z - m_grid.z0) / m_grid.dz + 0.5));
						if (k < 0 || k >= m_grid.nk)
						{
							continue;
						}
					}
					cells.push_back(Idx(i, j, k));
					weights.push_back(w);
				}
			}
		}
	}

	/// @brief 球–单元交集体积，按体积份额转成权重（Σ = 1 当粒子完全在域内）
	void ParticleGridMapper::BuildWeightsExact(const Particle& p,
		dynamic_array<long long>& cells, dynamic_array<double>& weights) const
	{
		const double r = p.radius;
		if (!(r > 0.0))
		{
			return;
		}
		const double particle_volume = (m_grid.dim >= 3)
			? (4.0 / 3.0) * kPi * r * r * r : kPi * r * r;
		if (!(particle_volume > 0.0))
		{
			return;
		}
		// 覆盖范围内的单元索引（注意：不能写成 std::max(int, double)，模板推导会失败）
		auto cell_of = [](double center, double h, double base, int n) -> int
		{
			const double idx = std::floor((center - base) / h + 0.5);
			const double v = std::min(static_cast<double>(n - 1),
				std::max(0.0, idx));
			return static_cast<int>(v);
		};
		const int i0 = cell_of(p.x - r, m_grid.dx, m_grid.x0, m_grid.ni);
		const int i1 = cell_of(p.x + r, m_grid.dx, m_grid.x0, m_grid.ni);
		const int j0 = cell_of(p.y - r, m_grid.dy, m_grid.y0, m_grid.nj);
		const int j1 = cell_of(p.y + r, m_grid.dy, m_grid.y0, m_grid.nj);
		int k0 = 0, k1 = 0;
		if (m_grid.dim >= 3)
		{
			k0 = cell_of(p.z - r, m_grid.dz, m_grid.z0, m_grid.nk);
			k1 = cell_of(p.z + r, m_grid.dz, m_grid.z0, m_grid.nk);
		}
		const double hx = 0.5 * m_grid.dx, hy = 0.5 * m_grid.dy, hz = 0.5 * m_grid.dz;
		for (int i = i0; i <= i1; ++i)
		{
			const double xc = m_grid.x0 + i * m_grid.dx;
			for (int j = j0; j <= j1; ++j)
			{
				const double yc = m_grid.y0 + j * m_grid.dy;
				for (int k = k0; k <= k1; ++k)
				{
					double v = 0.0;
					if (m_grid.dim >= 3)
					{
						const double zc = m_grid.z0 + k * m_grid.dz;
						v = SphereBoxVolume(r,
							xc - p.x - hx, xc - p.x + hx,
							yc - p.y - hy, yc - p.y + hy,
							zc - p.z - hz, zc - p.z + hz);
					}
					else
					{
						v = DiskRectArea(r,
							xc - p.x - hx, xc - p.x + hx,
							yc - p.y - hy, yc - p.y + hy);
					}
					if (v > 0.0)
					{
						cells.push_back(Idx(i, j, k));
						weights.push_back(v / particle_volume);
					}
				}
			}
		}
	}

	// ==================================================================
	// 固相体积分数
	// ==================================================================
	double ParticleGridMapper::ComputeSolidFraction(const dynamic_array<Particle>& particles,
		dynamic_array<double>& alpha)
	{
		const long long ncells = m_grid.CellCount();
		alpha.assign(static_cast<size_t>(ncells), 0.0);
		m_stats = Stats();

		const double inv_cell_volume = 1.0 / m_grid.CellVolume();
		double assigned = 0.0;
		double particle_volume = 0.0;
		long long pairs = 0;

		for (size_t ip = 0; ip < particles.size(); ++ip)
		{
			const Particle& p = particles[ip];
			const double vp = (m_grid.dim >= 3)
				? (4.0 / 3.0) * kPi * p.radius * p.radius * p.radius
				: kPi * p.radius * p.radius;
			if (!(vp > 0.0))
			{
				continue;
			}
			particle_volume += vp;
			BuildWeights(p, m_buf_cells, m_buf_weights);
			pairs += static_cast<long long>(m_buf_cells.size());
			for (size_t q = 0; q < m_buf_cells.size(); ++q)
			{
				const double dv = m_buf_weights[q] * vp;
				alpha[static_cast<size_t>(m_buf_cells[q])] += dv * inv_cell_volume;
				assigned += dv;
			}
		}

		// 保守光顺：先留一份原始场，再通量型 Laplacian，最后按 blend 混合
		if (m_smooth_passes > 0 && m_smooth_blend > 0.0)
		{
			double sum_before = 0.0;
			for (double a : alpha)
			{
				sum_before += a;
			}
			m_buf_field = alpha;
			SmoothInPlace(alpha);
			if (m_smooth_blend < 1.0)
			{
				for (size_t c = 0; c < alpha.size(); ++c)
				{
					alpha[c] = (1.0 - m_smooth_blend) * m_buf_field[c] + m_smooth_blend * alpha[c];
				}
			}
			double sum_after = 0.0;
			for (double a : alpha)
			{
				sum_after += a;
			}
			m_stats.smoothing_shift = (sum_before != 0.0)
				? std::fabs(sum_after / sum_before - 1.0) : 0.0;
		}

		// 统计
		double amin = 1.0e300, amax = -1.0e300;
		long long touched = 0;
		for (double a : alpha)
		{
			amin = std::min(amin, a);
			amax = std::max(amax, a);
			if (a > 0.0)
			{
				++touched;
			}
		}
		double max_jump = 0.0;
		for (int k = 0; k < m_grid.nk; ++k)
		{
			for (int j = 0; j < m_grid.nj; ++j)
			{
				for (int i = 0; i < m_grid.ni; ++i)
				{
					const double a0 = alpha[static_cast<size_t>(Idx(i, j, k))];
					if (i + 1 < m_grid.ni)
					{
						max_jump = std::max(max_jump,
							std::fabs(a0 - alpha[static_cast<size_t>(Idx(i + 1, j, k))]));
					}
					if (j + 1 < m_grid.nj)
					{
						max_jump = std::max(max_jump,
							std::fabs(a0 - alpha[static_cast<size_t>(Idx(i, j + 1, k))]));
					}
					if (k + 1 < m_grid.nk)
					{
						max_jump = std::max(max_jump,
							std::fabs(a0 - alpha[static_cast<size_t>(Idx(i, j, k + 1))]));
					}
				}
			}
		}

		m_stats.assigned_volume = assigned;
		m_stats.particle_volume = particle_volume;
		m_stats.alpha_min = (ncells > 0) ? amin : 0.0;
		m_stats.alpha_max = (ncells > 0) ? amax : 0.0;
		m_stats.max_jump = max_jump;
		m_stats.cells_touched = touched;
		m_stats.particle_cell_pairs = pairs;
		return assigned;
	}

	/// @brief 通量型 Laplacian 光顺：每个面通量成对加减同一个数 ⇒ Σα 逐位不变。
	void ParticleGridMapper::SmoothInPlace(dynamic_array<double>& alpha) const
	{
		const double k = 1.0 / (2.0 * static_cast<double>(m_grid.dim));
		const size_t ncells = alpha.size();
		dynamic_array<double> delta(ncells, 0.0);

		auto add_face = [&](long long c0, long long c1)
		{
			const double f = (alpha[static_cast<size_t>(c0)] - alpha[static_cast<size_t>(c1)]) * k;
			delta[static_cast<size_t>(c0)] -= f;
			delta[static_cast<size_t>(c1)] += f;
		};

		for (int pass = 0; pass < m_smooth_passes; ++pass)
		{
			std::fill(delta.begin(), delta.end(), 0.0);
			for (int kk = 0; kk < m_grid.nk; ++kk)
			{
				for (int jj = 0; jj < m_grid.nj; ++jj)
				{
					for (int ii = 0; ii < m_grid.ni; ++ii)
					{
						const long long c0 = Idx(ii, jj, kk);
						if (ii + 1 < m_grid.ni) add_face(c0, Idx(ii + 1, jj, kk));
						if (jj + 1 < m_grid.nj) add_face(c0, Idx(ii, jj + 1, kk));
						if (kk + 1 < m_grid.nk) add_face(c0, Idx(ii, jj, kk + 1));
					}
				}
			}
			for (size_t c = 0; c < ncells; ++c)
			{
				alpha[c] += delta[c];
			}
		}
	}

	// ==================================================================
	// 力散列与其对偶插值
	// ==================================================================
	void ParticleGridMapper::ScatterParticleVector(const dynamic_array<Particle>& particles,
		const dynamic_array<double>& particle_values,
		dynamic_array<double>& cell_values) const
	{
		if (particle_values.size() < particles.size() * 3)
		{
			throw ZaranError("ParticleGridMapper::ScatterParticleVector: particle_values too short");
		}
		cell_values.assign(static_cast<size_t>(m_grid.CellCount()) * 3, 0.0);
		dynamic_array<long long> cells;
		dynamic_array<double> weights;
		for (size_t ip = 0; ip < particles.size(); ++ip)
		{
			BuildWeights(particles[ip], cells, weights);
			if (cells.empty())
			{
				continue;
			}
			for (size_t q = 0; q < cells.size(); ++q)
			{
				const double w = weights[q];
				double* dst = &cell_values[static_cast<size_t>(cells[q]) * 3];
				dst[0] += w * particle_values[ip * 3 + 0];
				dst[1] += w * particle_values[ip * 3 + 1];
				dst[2] += w * particle_values[ip * 3 + 2];
			}
		}
	}

	void ParticleGridMapper::GatherToParticle(const dynamic_array<Particle>& particles,
		const dynamic_array<double>& cell_values,
		dynamic_array<double>& particle_values) const
	{
		if (cell_values.size() < static_cast<size_t>(m_grid.CellCount()) * 3)
		{
			throw ZaranError("ParticleGridMapper::GatherToParticle: cell_values too short");
		}
		particle_values.assign(particles.size() * 3, 0.0);
		dynamic_array<long long> cells;
		dynamic_array<double> weights;
		for (size_t ip = 0; ip < particles.size(); ++ip)
		{
			BuildWeights(particles[ip], cells, weights);
			double acc[3] = { 0.0, 0.0, 0.0 };
			for (size_t q = 0; q < cells.size(); ++q)
			{
				const double* src = &cell_values[static_cast<size_t>(cells[q]) * 3];
				const double w = weights[q];
				acc[0] += w * src[0];
				acc[1] += w * src[1];
				acc[2] += w * src[2];
			}
			particle_values[ip * 3 + 0] = acc[0];
			particle_values[ip * 3 + 1] = acc[1];
			particle_values[ip * 3 + 2] = acc[2];
		}
	}
}
