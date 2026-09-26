/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file DEMCFDCoupler.cpp
 * \brief One-way DEM-CFD coupling: particle -> grid voidage, grid -> particle gas state,
 *        and the inter-phase force on every particle.
 */

#include "DEMCFDCoupler.h"
#include "Log.h"
#include "ZaranError.h"

#include <algorithm>
#include <cmath>

namespace zaran
{
	namespace
	{
		int ClampInt(int v, int lo, int hi)
		{
			return (v < lo) ? lo : ((v > hi) ? hi : v);
		}
	}

	// ==================================================================
	// 初始化
	// ==================================================================
	void DEMCFDCoupler::Init(const GridView& grid, const UnitScales& scales, const Options& options)
	{
		if (!grid.Valid())
		{
			throw ZaranError("DEMCFDCoupler: invalid CFD grid view "
				"(only 3-D uniform structured grids are supported)");
		}
		if (!(scales.ref_length > 0.0) || !(scales.ref_density > 0.0)
			|| !(scales.ref_velocity > 0.0) || !(scales.ref_temperature > 0.0))
		{
			throw ZaranError("DEMCFDCoupler: non-positive reference value in UnitScales");
		}
		m_grid = grid;
		m_scales = scales;
		m_options = options;

		// 映射器网格 = CFD 的**物理节点**范围：
		//   CFD 物理节点 (is + ci, js + cj, ks + ck) ←→ 映射器单元 (ci, cj, ck)
		ParticleGridMapper::GridSpec spec;
		spec.dim = 3;
		spec.ni = grid.PhysicalNi();
		spec.nj = grid.PhysicalNj();
		spec.nk = grid.PhysicalNk();
		spec.dx = grid.dx;
		spec.dy = grid.dy;
		spec.dz = grid.dz;
		spec.x0 = grid.x0;
		spec.y0 = grid.y0;
		spec.z0 = grid.z0;
		m_mapper.Init(spec);
		m_mapper.SetMethod(options.mapping);
		m_mapper.SetSubdivisions(options.mapping_sub);
		m_mapper.SetSmoothing(options.smooth_passes, options.smooth_blend);

		double lo[3], hi[3];
		GetControlVolumeBox(lo, hi);
		Log::info("DEMCFDCoupler: CFD nodes (with ghost) = ({},{},{}), physical range "
			"= i[{}..{}] j[{}..{}] k[{}..{}], spacing = ({:E},{:E},{:E})",
			grid.ni, grid.nj, grid.nk, grid.is, grid.ie, grid.js, grid.je, grid.ks, grid.ke,
			grid.dx, grid.dy, grid.dz);
		Log::info("DEMCFDCoupler: 控制体范围 = [{:E},{:E}] x [{:E},{:E}] x [{:E},{:E}] m "
			"(比节点跨度多出半格)", lo[0], hi[0], lo[1], hi[1], lo[2], hi[2]);
		Log::info("DEMCFDCoupler: 曳力模型 = {}（稠密阈值 ε_g={}, 稀疏支 Re {}空隙率, "
			"压力梯度力 {}）, 体积分数 {}",
			DragForceModel::ModelName(options.model), options.drag.dense_threshold,
			options.drag.dilute_re_uses_voidage ? "含" : "不含",
			options.drag.include_pressure_gradient_force ? "开" : "关",
			options.volume_fraction_from_particles ? "由粒子给出" : "恒为 1（无体积排斥）");
	}

	int DEMCFDCoupler::SetParticles(const std::vector<Particle>& particles)
	{
		m_particles.clear();
		m_particles.reserve(particles.size());
		m_particle_volume = 0.0;
		m_rep_diameter = 1.0e30;
		m_outside_count = 0;

		for (const auto& p : particles)
		{
			if (!(p.radius > 0.0) || !(p.density > 0.0))
			{
				Log::warn("DEMCFDCoupler: 跳过非法粒子 id={}（radius={:E}, density={:E}）",
					p.id, p.radius, p.density);
				continue;
			}
			m_particles.push_back(p);
			m_particle_volume += ParticleGridMapper::ParticleVolume(3, p.radius);
			m_rep_diameter = std::min(m_rep_diameter, 2.0 * p.radius);
			if (!ParticleInside(p))
			{
				++m_outside_count;
			}
		}
		if (m_particles.empty())
		{
			throw ZaranError("DEMCFDCoupler: no valid particle");
		}
		if (m_rep_diameter <= 0.0)
		{
			m_rep_diameter = 1.0;
		}

		double lo[3], hi[3];
		GetControlVolumeBox(lo, hi);
		double ratios[3];
		GetResolutionRatios(ratios);
		Log::info("DEMCFDCoupler: {} 个粒子，ΣV_p = {:E} m^3，d_p(min) = {:E} m，"
			"质心在计算域外 = {}",
			m_particles.size(), m_particle_volume, m_rep_diameter, m_outside_count);
		Log::info("DEMCFDCoupler: Δx/d_p = {:.4f}, Δy/d_p = {:.4f}, Δz/d_p = {:.4f} "
			"（体积分数法要求 ≳ 3）", ratios[0], ratios[1], ratios[2]);
		if (m_outside_count > 0)
		{
			Log::warn("DEMCFDCoupler: {} 个粒子的质心落在控制体范围外，这部分体积不会落到网格上，"
				"体积分数会偏小、该粒子的受力插值权重也会小于 1", m_outside_count);
		}
		return m_outside_count;
	}

	bool DEMCFDCoupler::ParticleInside(const Particle& p) const
	{
		double lo[3], hi[3];
		GetControlVolumeBox(lo, hi);
		const double c[3] = { p.x, p.y, p.z };
		for (int d = 0; d < 3; ++d)
		{
			if (c[d] < lo[d] || c[d] > hi[d])
			{
				return false;
			}
		}
		return true;
	}

	void DEMCFDCoupler::GetControlVolumeBox(double lo[3], double hi[3]) const
	{
		const double h[3] = { m_grid.dx, m_grid.dy, m_grid.dz };
		const double c0[3] = { m_grid.x0, m_grid.y0, m_grid.z0 };
		const int n[3] = { m_grid.PhysicalNi(), m_grid.PhysicalNj(), m_grid.PhysicalNk() };
		for (int d = 0; d < 3; ++d)
		{
			lo[d] = c0[d] - 0.5 * h[d];
			hi[d] = c0[d] + (n[d] - 0.5) * h[d];
		}
	}

	void DEMCFDCoupler::GetResolutionRatios(double ratios[3]) const
	{
		ratios[0] = m_grid.dx / m_rep_diameter;
		ratios[1] = m_grid.dy / m_rep_diameter;
		ratios[2] = m_grid.dz / m_rep_diameter;
	}

	// ==================================================================
	// DEM → CFD：固相体积分数 → 气相体积分数
	// ==================================================================
	void DEMCFDCoupler::ComputeGasVolumeFraction(std::vector<double>& eps_nodes,
		double& eps_min, double& eps_max, bool verbose) const
	{
		const long long node_num = m_grid.NodeCount();
		eps_nodes.assign(static_cast<size_t>(node_num), 1.0);

		// ε ≡ 1 的退化模式：整个数组恒为 1，CFD 侧数值上等价于单相求解器
		if (!m_options.volume_fraction_from_particles)
		{
			eps_min = 1.0;
			eps_max = 1.0;
			if (verbose)
			{
				Log::info("DEMCFDCoupler: volume_fraction_from_particles = false ⇒ "
					"ε_g ≡ 1（严格单向耦合，CFD 与单相 Euler 逐位退化）");
			}
			return;
		}

		// 映射器：粒子 → 固相体积分数 α_s（体积精确配分：Σ α_s V_cell = Σ V_p）
		std::vector<ParticleGridMapper::Particle> pts;
		pts.reserve(m_particles.size());
		for (const auto& p : m_particles)
		{
			ParticleGridMapper::Particle q;
			q.x = p.x;
			q.y = p.y;
			q.z = p.z;
			q.radius = p.radius;
			pts.push_back(q);
		}
		m_mapper.ComputeSolidFraction(pts, m_alpha);

		const ParticleGridMapper::Stats& stats = m_mapper.GetStats();
		eps_min = 1.0;
		eps_max = 0.0;
		long long written = 0;
		const double eps_floor = m_options.eps_min_clip;
		for (int ck = 0; ck < m_grid.PhysicalNk(); ++ck)
		{
			for (int cj = 0; cj < m_grid.PhysicalNj(); ++cj)
			{
				for (int ci = 0; ci < m_grid.PhysicalNi(); ++ci)
				{
					const long long cell = m_mapper.Idx(ci, cj, ck);
					double eps = 1.0 - m_alpha[static_cast<size_t>(cell)];
					if (eps < eps_floor)
					{
						eps = eps_floor;
					}
					if (eps > 1.0)
					{
						eps = 1.0;
					}
					eps_nodes[static_cast<size_t>(m_grid.Idx(m_grid.is + ci, m_grid.js + cj,
						m_grid.ks + ck))] = eps;
					eps_min = std::min(eps_min, eps);
					eps_max = std::max(eps_max, eps);
					++written;
				}
			}
		}

		if (verbose)
		{
			const double resid = (m_particle_volume > 0.0)
				? std::fabs(stats.assigned_volume - m_particle_volume) / m_particle_volume
				: 0.0;
			Log::info("DEMCFDCoupler: 气相体积分数写毕，{} 个物理节点，ε_g ∈ [{:.6f}, {:.6f}]，"
				"固相体积分数 α_max = {:.6f}", written, eps_min, eps_max, stats.alpha_max);
			Log::info("DEMCFDCoupler: 体积配分残差 |ΣαV − ΣV_p|/ΣV_p = {:.3e}", resid);
		}
	}

	// ==================================================================
	// 双向耦合：反作用力散列 + 颗粒推进
	// ==================================================================
	double DEMCFDCoupler::ParticleMass(size_t index) const
	{
		if (index >= m_particles.size())
		{
			return 0.0;
		}
		const Particle& p = m_particles[index];
		return p.density * ParticleGridMapper::ParticleVolume(3, p.radius);
	}

	void DEMCFDCoupler::SumSourceMomentum(const std::vector<double>& src_dimensionless,
		double out[3]) const
	{
		out[0] = out[1] = out[2] = 0.0;
		if (src_dimensionless.size() != static_cast<size_t>(3 * m_grid.NodeCount()))
		{
			return;
		}
		const double v_cell_star = m_grid.CellVolume()
			/ (m_scales.ref_length * m_scales.ref_length * m_scales.ref_length);
		for (int ck = 0; ck < m_grid.PhysicalNk(); ++ck)
		{
			for (int cj = 0; cj < m_grid.PhysicalNj(); ++cj)
			{
				for (int ci = 0; ci < m_grid.PhysicalNi(); ++ci)
				{
					const index_type idx = m_grid.Idx(m_grid.is + ci, m_grid.js + cj,
						m_grid.ks + ck);
					for (int d = 0; d < 3; ++d)
					{
						out[d] += src_dimensionless[static_cast<size_t>(3 * idx + d)]
							* v_cell_star;
					}
				}
			}
		}
	}

	void DEMCFDCoupler::ScatterReaction(const std::vector<ParticleResult>& results,
		bool include_pressure_gradient_force,
		std::vector<double>& src_dimensionless, ReactionSummary& summary) const
	{
		summary = ReactionSummary();
		src_dimensionless.assign(static_cast<size_t>(3 * m_grid.NodeCount()), 0.0);

		const double l3 = m_scales.ref_length * m_scales.ref_length * m_scales.ref_length;
		const double v_cell_si = m_grid.CellVolume();
		if (!(v_cell_si > 0.0))
		{
			throw ZaranError("DEMCFDCoupler::ScatterReaction: zero cell volume");
		}
		// S* = F* / V*_cell，其中 F* = F_SI / F_ref，V*_cell = V_cell / L_ref³
		//    ⇒ S* = F_SI · L_ref³ / (F_ref · V_cell)
		// ⚠ 这里只能对 F_ref **除一次**：写成 `ForceToDimensionless(f) * (L³/(F_ref·V))`
		//   会变成 1/F_ref²，反作用力被缩小 5 个量级（实测：源冲量 4.8e-13 vs
		//   颗粒冲量 5.2e-08）。判据就是"散列进去的动量 = 颗粒受力的负值"。
		const double factor = l3 / (m_scales.RefForce() * v_cell_si);

		const bool has_grad = (m_grad_p.size() == static_cast<size_t>(3 * m_grid.NodeCount()));

		for (size_t ip = 0; ip < m_particles.size() && ip < results.size(); ++ip)
		{
			const Particle& p = m_particles[ip];
			const DragForceModel::Result& r = results[ip].force;

			// 流体施加在颗粒上的合力（反作用要取负）
			double f_si[3];
			for (int d = 0; d < 3; ++d)
			{
				f_si[d] = r.drag[d];
			}
			if (include_pressure_gradient_force)
			{
				// 只有气相方程取完全守恒的 ∇(εp) 形式时才需要把浮力也反作用回去
				for (int d = 0; d < 3; ++d)
				{
					f_si[d] += r.pressure_gradient[d];
				}
			}
			(void)has_grad;
			for (int d = 0; d < 3; ++d)
			{
				summary.particle_force_dimensionless[d] += ForceToDimensionless(f_si[d]);
				summary.particle_force_l1 += std::fabs(ForceToDimensionless(f_si[d]));
			}

			ParticleGridMapper::Particle q;
			q.x = p.x;
			q.y = p.y;
			q.z = p.z;
			q.radius = p.radius;
			m_mapper.ComputeWeights(q, m_buf_cells, m_buf_weights);

			double w_sum = 0.0;
			for (size_t t = 0; t < m_buf_cells.size(); ++t)
			{
				const double w = m_buf_weights[t];
				w_sum += w;
				const long long cell = m_buf_cells[t];
				const int ci = static_cast<int>(cell % m_grid.PhysicalNi());
				const int cj = static_cast<int>((cell / m_grid.PhysicalNi())
					% m_grid.PhysicalNj());
				const int ck = static_cast<int>(cell / (static_cast<long long>(
					m_grid.PhysicalNi()) * m_grid.PhysicalNj()));
				const index_type idx = m_grid.Idx(m_grid.is + ci, m_grid.js + cj,
					m_grid.ks + ck);
				for (int d = 0; d < 3; ++d)
				{
					src_dimensionless[static_cast<size_t>(3 * idx + d)]
						-= w * f_si[d] * factor;
				}
			}
			summary.weight_deficit += std::max(0.0, 1.0 - w_sum);
		}

		SumSourceMomentum(src_dimensionless, summary.scattered_force_dimensionless);
	}

	int DEMCFDCoupler::AdvanceParticles(double dt, const MotionOptions& options,
		const std::vector<ParticleResult>& results)
	{
		if (options.mode != MotionMode::Dynamic || !(dt > 0.0))
		{
			return 0;
		}
		// ⚠ 求解器内部的 dt 是**无量纲**的（t* = t·a_ref/L_ref），而颗粒状态是 SI。
		//   直接拿无量纲 dt 去乘 SI 速度，位移会被放大 a_ref/L_ref ≈ 3.3e2 倍
		//   （实测：颗粒在 0.15 的无量纲时间里"跑"了 0.34 m，直接飞出计算域）。
		const double dt_si = dt * (m_scales.ref_length / m_scales.ref_velocity);
		double max_sub = 1.0;
		for (size_t i = 0; i < m_particles.size(); ++i)
		{
			const Particle& p = m_particles[i];
			const double m = p.density * ParticleGridMapper::ParticleVolume(3, p.radius);
			if (!(m > 0.0) || i >= results.size())
			{
				continue;
			}
			const double vp = ParticleGridMapper::ParticleVolume(3, p.radius);
			const double k = vp * results[i].force.beta_over_es / m;   // [1/s]
			max_sub = std::max(max_sub, k * dt_si / std::max(1.0e-12,
				options.max_step_parameter));
		}
		int n_sub = options.substeps > 0 ? options.substeps
			: static_cast<int>(std::ceil(max_sub));
		n_sub = std::min(std::max(n_sub, 1), 100000);
		const double h = dt_si / static_cast<double>(n_sub);

		for (size_t i = 0; i < m_particles.size(); ++i)
		{
			Particle& p = m_particles[i];
			if (i >= results.size())
			{
				continue;
			}
			const GasSample& g = results[i].gas;
			const DragForceModel::Result& r = results[i].force;
			const double vp = ParticleGridMapper::ParticleVolume(3, p.radius);
			const double m = p.density * vp;
			if (!(m > 0.0))
			{
				continue;
			}
			// 曳力线性化系数 k = V_p·(β/ε_s)/m 与"外力"加速度 a_ext
			const double k = vp * r.beta_over_es / m;
			double a_ext[3];
			for (int d = 0; d < 3; ++d)
			{
				a_ext[d] = r.pressure_gradient[d] / m + options.body_acceleration[d];
			}
			for (int s = 0; s < n_sub; ++s)
			{
				const double u_old[3] = { p.velocity[0], p.velocity[1], p.velocity[2] };
				for (int d = 0; d < 3; ++d)
				{
					if (options.integrator == Integrator::Rk2)
					{
						const double k1 = k * (g.velocity[d] - p.velocity[d]) + a_ext[d];
						const double u_mid = p.velocity[d] + 0.5 * h * k1;
						const double k2 = k * (g.velocity[d] - u_mid) + a_ext[d];
						p.velocity[d] += h * k2;
					}
					else
					{
						// 对线性 ODE du/dt = k(u_g − u) + a_ext 的解析解
						const double kh = k * h;
						if (kh > 1.0e-8)
						{
							const double u_inf = g.velocity[d] + a_ext[d] / k;
							p.velocity[d] = u_inf
								+ (p.velocity[d] - u_inf) * std::exp(-kh);
						}
						else
						{
							p.velocity[d] += (k * (g.velocity[d] - p.velocity[d])
								+ a_ext[d]) * h;
						}
					}
				}
				// 位置用梯形法（与速度的精度同阶）
				p.x += h * 0.5 * (u_old[0] + p.velocity[0]);
				p.y += h * 0.5 * (u_old[1] + p.velocity[1]);
				p.z += h * 0.5 * (u_old[2] + p.velocity[2]);
			}
			// 跑出计算域必须点名：体积分数与反作用都会跟着少算
			if (!ParticleInside(p))
			{
				if (m_warned_outside.insert(p.id).second)
				{
					Log::warn("DEMCFDCoupler: 颗粒 id={} 已运动到计算域外（x={:.6f}, y={:.6f}, z={:.6f}）——"
						"它的体积分数与反作用力都会被少算，请减小 t_end 或调整初值",
						p.id, p.x, p.y, p.z);
				}
			}
		}
		return n_sub;
	}


	// ==================================================================
	// 压强梯度
	// ==================================================================
	void DEMCFDCoupler::ComputePressureGradient(const double* pressure_dimensionless,
		std::vector<double>& grad_p_dimensionless) const
	{
		if (pressure_dimensionless == nullptr)
		{
			throw ZaranError("DEMCFDCoupler::ComputePressureGradient: null pressure");
		}
		const long long node_num = m_grid.NodeCount();
		grad_p_dimensionless.assign(static_cast<size_t>(3 * node_num), 0.0);

		const double h[3] = { m_grid.dx, m_grid.dy, m_grid.dz };
		const int lo[3] = { m_grid.is, m_grid.js, m_grid.ks };
		const int hi[3] = { m_grid.ie, m_grid.je, m_grid.ke };

		for (int k = m_grid.ks; k <= m_grid.ke; ++k)
		{
			for (int j = m_grid.js; j <= m_grid.je; ++j)
			{
				for (int i = m_grid.is; i <= m_grid.ie; ++i)
				{
					const long long idx = m_grid.Idx(i, j, k);
					int c[3] = { i, j, k };
					for (int dir = 0; dir < 3; ++dir)
					{
						int ap[3] = { c[0], c[1], c[2] };
						int am[3] = { c[0], c[1], c[2] };
						ap[dir] = ClampInt(c[dir] + 1, lo[dir], hi[dir]);
						am[dir] = ClampInt(c[dir] - 1, lo[dir], hi[dir]);
						// 物理范围边界处退化为单侧差分（不依赖 ghost 的边界处理）
						const double denom = static_cast<double>(ap[dir] - am[dir]) * h[dir];
						const double dp = pressure_dimensionless[m_grid.Idx(ap[0], ap[1], ap[2])]
							- pressure_dimensionless[m_grid.Idx(am[0], am[1], am[2])];
						grad_p_dimensionless[static_cast<size_t>(3 * idx + dir)] = dp / denom;
					}
				}
			}
		}
	}

	// ==================================================================
	// CFD → DEM：插值 + 相间力
	// ==================================================================
	void DEMCFDCoupler::EvaluateParticleForces(const PrimitiveView& prim,
		std::vector<ParticleResult>& outputs) const
	{
		if (!prim.Valid() || prim.viscosity == nullptr)
		{
			throw ZaranError("DEMCFDCoupler::EvaluateParticleForces: incomplete primitive view");
		}
		outputs.assign(m_particles.size(), ParticleResult());

		// ∇p：先算到 CFD 节点上，再按同一套权重插值到粒子位置。
		// 只有需要压力梯度力时才做（单向耦合下它与 ε 一样是"每步的流场后处理"）。
		if (m_options.drag.include_pressure_gradient_force)
		{
			ComputePressureGradient(prim.pressure, m_grad_p);
		}
		const bool have_grad_p = (m_grad_p.size() == static_cast<size_t>(3 * m_grid.NodeCount()));
		const bool use_particles = m_options.volume_fraction_from_particles;
		if (use_particles && m_alpha.empty())
		{
			throw ZaranError("DEMCFDCoupler::EvaluateParticleForces: "
				"ComputeGasVolumeFraction() must be called first");
		}

		const double fac_rho = m_scales.ref_density;
		const double fac_vel = m_scales.ref_velocity;
		const double fac_p = m_scales.RefPressure();
		const double fac_mu = m_scales.RefViscosity();
		const double fac_grad_p = fac_p / m_scales.ref_length;

		for (size_t ip = 0; ip < m_particles.size(); ++ip)
		{
			const Particle& p = m_particles[ip];

			ParticleGridMapper::Particle q;
			q.x = p.x;
			q.y = p.y;
			q.z = p.z;
			q.radius = p.radius;
			m_mapper.ComputeWeights(q, m_buf_cells, m_buf_weights);

			GasSample s;
			double rho = 0.0, p_dim = 0.0, mu_dim = 0.0, eps = 0.0;
			double vel[3] = { 0.0, 0.0, 0.0 };
			double gp[3] = { 0.0, 0.0, 0.0 };

			for (size_t t = 0; t < m_buf_cells.size(); ++t)
			{
				const double w = m_buf_weights[t];
				const long long cell = m_buf_cells[t];
				// 映射器单元 → CFD 物理节点
				const int ci = static_cast<int>(cell % m_grid.PhysicalNi());
				const int cj = static_cast<int>((cell / m_grid.PhysicalNi()) % m_grid.PhysicalNj());
				const int ck = static_cast<int>(cell / (static_cast<long long>(m_grid.PhysicalNi())
					* m_grid.PhysicalNj()));
				const long long idx = m_grid.Idx(m_grid.is + ci, m_grid.js + cj, m_grid.ks + ck);

				rho += w * prim.density[idx];
				for (int d = 0; d < 3; ++d)
				{
					vel[d] += w * prim.velocity[d][idx];
				}
				p_dim += w * prim.pressure[idx];
				mu_dim += w * prim.viscosity[idx];
				eps += use_particles
					? w * (1.0 - m_alpha[static_cast<size_t>(cell)])
					: w;
				if (have_grad_p)
				{
					for (int d = 0; d < 3; ++d)
					{
						gp[d] += w * m_grad_p[static_cast<size_t>(3 * idx + d)];
					}
				}
				s.weight_sum += w;
			}
			s.cells_touched = static_cast<long long>(m_buf_cells.size());

			// ---- 无量纲 → SI ----
			s.density = rho * fac_rho;
			for (int d = 0; d < 3; ++d)
			{
				s.velocity[d] = vel[d] * fac_vel;
			}
			s.pressure = p_dim * fac_p;
			// T* = γ p*/ρ*（与 Gas::CalcTemperature 同一口径），再乘 T_ref
			s.temperature = (rho > 0.0)
				? (m_scales.ref_gamma * p_dim / rho) * m_scales.ref_temperature
				: 0.0;
			s.viscosity = mu_dim * fac_mu;
			s.eps_gas = (eps > 0.0) ? std::min(1.0, eps) : 1.0;
			for (int d = 0; d < 3; ++d)
			{
				s.grad_p[d] = gp[d] * fac_grad_p;
			}

			// ---- 相间力（SI 下评估 ⇒ 输出直接是 N） ----
			DragForceModel::GasState gas;
			gas.density = s.density;
			gas.pressure = s.pressure;
			gas.temperature = s.temperature;
			gas.viscosity = s.viscosity;
			for (int d = 0; d < 3; ++d)
			{
				gas.velocity[d] = s.velocity[d];
			}
			DragForceModel::ParticleState part;
			part.radius = p.radius;
			part.density = p.density;
			for (int d = 0; d < 3; ++d)
			{
				part.velocity[d] = p.velocity[d];
			}

			outputs[ip].gas = s;
			outputs[ip].force = DragForceModel::Evaluate(m_options.model, m_options.drag,
				gas, part, s.eps_gas, s.grad_p);
		}
	}
}
