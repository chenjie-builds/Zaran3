/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file EulerTwoPhaseStructUniform.cpp
 * \brief Two-phase (volume-fraction weighted) Euler solver on a uniform Cartesian grid.
 */

#include "EulerTwoPhaseStructUniform.h"
#include "Log.h"
#include "ZaranError.h"
#include "MathBasic.h"

#include <cmath>

namespace zaran
{
	namespace
	{
		constexpr int kEqNum = 5;
		/// @brief 动量方程在守恒量里的起始下标（rho*u, rho*v, rho*w）
		constexpr int kMomentumOffset = 1;
		/// @brief 能量方程在守恒量里的下标
		constexpr int kEnergyOffset = 4;

		/// @brief 不含 ghost 的计算范围，统一转成 int 供 OpenMP 使用
		///        （MSVC 对循环条件里直接出现 static_cast 会报 C1001）
		struct LoopBounds
		{
			int is = 0, ie = 0, js = 0, je = 0, ks = 0, ke = 0;

			explicit LoopBounds(GridStruct& grid)
			{
				index_type a, b, c, d, e, f;
				grid.GetRange(a, b, c, d, e, f);
				is = static_cast<int>(a);
				ie = static_cast<int>(b);
				js = static_cast<int>(c);
				je = static_cast<int>(d);
				ks = static_cast<int>(e);
				ke = static_cast<int>(f);
			}
		};

		/// @brief 全节点范围（含 ghost）
		struct NodeCounts
		{
			int ni = 0, nj = 0, nk = 0;

			explicit NodeCounts(GridStruct& grid)
			{
				ni = static_cast<int>(grid.GetNi());
				nj = static_cast<int>(grid.GetNj());
				nk = static_cast<int>(grid.GetNk());
			}
		};
	}

	// ==================================================================
	// 构造与访问
	// ==================================================================
	EulerTwoPhaseStructUniform::EulerTwoPhaseStructUniform(index_type index, string name,
		shared_ptr<FlowSolverParam> para, shared_ptr<GridBase> grid,
		shared_ptr<DataManagerNS> data_manager)
		: EulerSolverStructUniform(index, name, para, grid, data_manager)
	{
		m_tp_data_manager = std::dynamic_pointer_cast<DataManagerNSTwoPhase>(data_manager);
		if (!m_tp_data_manager)
		{
			throw ZaranError("EulerTwoPhaseStructUniform requires DataManagerNSTwoPhase");
		}
		m_tp_para = std::dynamic_pointer_cast<FlowSolverParamTwoPhase>(para);
		if (!m_tp_para)
		{
			throw ZaranError("EulerTwoPhaseStructUniform requires FlowSolverParamTwoPhase");
		}
	}

	EulerTwoPhaseStructUniform::~EulerTwoPhaseStructUniform()
	{
	}

	DataManagerNSTwoPhase* EulerTwoPhaseStructUniform::GetTwoPhaseDataManager()
	{
		return m_tp_data_manager.get();
	}

	FlowSolverParamTwoPhase* EulerTwoPhaseStructUniform::GetTwoPhasePara()
	{
		return m_tp_para.get();
	}

	const double* EulerTwoPhaseStructUniform::GetVolumeFraction() const
	{
		return m_tp_data_manager->GetVolumeFraction();
	}

	const double* EulerTwoPhaseStructUniform::GetVolumeFractionField() const
	{
		// 基类的通量加权与变量转换都通过这个钩子取 ε
		return m_tp_data_manager->GetVolumeFraction();
	}

	// ==================================================================
	// 体积分数初值
	// ==================================================================
	double EulerTwoPhaseStructUniform::EvaluateVolumeFraction(double x, double y, double z) const
	{
		(void)y;
		(void)z;
		const auto type = m_tp_para->GetVolumeFractionType();
		if (type == VolumeFractionType::Step)
		{
			return (x <= m_tp_para->GetVolumeFractionXStep())
				? m_tp_para->GetVolumeFractionLeft()
				: m_tp_para->GetVolumeFractionRight();
		}
		if (type == VolumeFractionType::Sine)
		{
			const double k = 2.0 * PI / m_tp_para->GetVolumeFractionWavelength();
			return m_tp_para->GetVolumeFractionMean()
				+ m_tp_para->GetVolumeFractionAmplitude()
				* std::sin(k * x + m_tp_para->GetVolumeFractionPhase());
		}
		return m_tp_para->GetVolumeFractionValue();
	}

	void EulerTwoPhaseStructUniform::InitVolumeFraction()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto data_manager = GetTwoPhaseDataManager();
		IdProxyStruct& idx_proxy = GetIdxProxy();
		const NodeCounts counts(*GetGrid());
		const int ni = counts.ni, nj = counts.nj, nk = counts.nk;
		double* volume_fraction = data_manager->GetVolumeFraction();

		double eps_min = LARGE_NUMBER, eps_max = -LARGE_NUMBER;
		double x_first = 0.0, x_last = 0.0;
		for (int k = 0; k < nk; ++k)
		{
			for (int j = 0; j < nj; ++j)
			{
				for (int i = 0; i < ni; ++i)
				{
					const auto coord = node->GetCoord(i, j, k);
					const double eps = EvaluateVolumeFraction(coord[0], coord[1], coord[2]);
					volume_fraction[idx_proxy(i, j, k)] = eps;
					eps_min = Min(eps_min, eps);
					eps_max = Max(eps_max, eps);
					if (k == 0 && j == 0)
					{
						if (i == 0) { x_first = eps; }
						if (i == ni - 1) { x_last = eps; }
					}
				}
			}
		}

		Log::info("Two-phase volume fraction initialized: type={}, eps in [{:E}, {:E}], "
			"eps(x_min)={:E}, eps(x_max)={:E}",
			static_cast<int>(m_tp_para->GetVolumeFractionType()),
			eps_min, eps_max, x_first, x_last);
	}

	void EulerTwoPhaseStructUniform::FillVolumeFractionGhost()
	{
		auto grid = GetGrid();
		auto data_manager = GetTwoPhaseDataManager();
		IdProxyStruct& idx_proxy = GetIdxProxy();
		auto bound_map = grid->GetBoundMap();
		const size_t ghost_level = grid->GetGhostLevel();
		double* volume_fraction = data_manager->GetVolumeFraction();

		for (auto& boundary : bound_map->GetBoundMap())
		{
			const auto& bound_name = boundary.first;
			if (bound_name == "hole")
			{
				continue;
			}
			auto& bound = boundary.second;
			for (size_t iBound = 0; iBound < bound.size(); ++iBound)
			{
				index_type i_bound, j_bound, k_bound;
				bound[iBound].GetIdx(i_bound, j_bound, k_bound);
				const index_type idx_bound = idx_proxy(i_bound, j_bound, k_bound);
				const double eps_bound = volume_fraction[idx_bound];
				const auto direction = bound[iBound].GetDirectionSrc();
				for (size_t iGhost = 1; iGhost <= ghost_level; ++iGhost)
				{
					const index_type i_ghost = i_bound + iGhost * direction[0];
					const index_type j_ghost = j_bound + iGhost * direction[1];
					const index_type k_ghost = k_bound + iGhost * direction[2];
					volume_fraction[idx_proxy(i_ghost, j_ghost, k_ghost)] = eps_bound;
				}
			}
		}
	}

	void EulerTwoPhaseStructUniform::InitField()
	{
		// ε 必须先于基类初值铺好：基类 InitField() 末尾会调用虚函数 Prim2Cons()，
		// 而本类的 Prim2Cons 需要 ε 场（cons = ε · cons_intrinsic）。
		if (m_tp_para->GetVolumeFractionType() == VolumeFractionType::External)
		{
			// 阶段 3：ε 由 DEM-CFD 耦合器在 solver->Init() 之前写进数据管理器的
			// volume_fraction（物理节点）。这里只负责把 ghost 层补齐——
			// 不能直接调用基类 InitField()，否则它会用解析场把耦合器写好的值覆盖掉。
			FillVolumeFractionGhost();
			double eps_min = LARGE_NUMBER, eps_max = -LARGE_NUMBER;
			auto grid = GetGrid();
			IdProxyStruct& idx_proxy = GetIdxProxy();
			const NodeCounts counts(*GetGrid());
			index_type is, ie, js, je, ks, ke;
			grid->GetRange(is, ie, js, je, ks, ke);
			const double* volume_fraction = GetVolumeFraction();
			for (int k = static_cast<int>(ks); k <= static_cast<int>(ke); ++k)
			{
				for (int j = static_cast<int>(js); j <= static_cast<int>(je); ++j)
				{
					for (int i = static_cast<int>(is); i <= static_cast<int>(ie); ++i)
					{
						const double eps = volume_fraction[idx_proxy(i, j, k)];
						eps_min = Min(eps_min, eps);
						eps_max = Max(eps_max, eps);
					}
				}
			}
			Log::info("Two-phase volume fraction: type=external（由耦合器注入）, "
				"物理节点 ε_g ∈ [{:E}, {:E}]", eps_min, eps_max);
			if (eps_max > 1.0 + 1.0e-12)
			{
				Log::warn("EulerTwoPhaseStructUniform: 注入的 ε_g 上界 {:E} > 1，"
					"耦合侧应把 ε 限制在 [eps_min, 1]", eps_max);
			}
		}
		else
		{
			InitVolumeFraction();
		}
		EulerSolverStructUniform::InitField();
	}

	// ==================================================================
	// 边界与残差
	// ==================================================================
	void EulerTwoPhaseStructUniform::BoundaryCondition()
	{
		// ε 的虚拟节点必须先铺好：基类填充流场虚拟节点时要按 ε 加权存守恒量
		FillVolumeFractionGhost();
		EulerSolverStructUniform::BoundaryCondition();
	}

	/// @brief 孔隙率梯度动量源 +p∇ε。
	/// @details 动量方程取 ε∇p 形式（DEM-CFD 惯用），而左侧写的是守恒的
	///          ∇·(ε ρ u u + ε p I)，两者相差一项 p∇ε：
	///              ∂(ε ρ u)/∂t = -∇·(ε ρ u u + ε p I) + p∇ε
	///          ε 的中心差分与界面的算术平均恰好抵消（见文档），于是
	///          "压强均匀 + 静止" 是离散精确的定常解。
	void EulerTwoPhaseStructUniform::AddPorosityGradientSource()
	{
		if (m_tp_para->GetPorosityGradientForce() == 0)
		{
			return;
		}
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		IdProxyStruct& idx_proxy = GetIdxProxy();
		const LoopBounds bounds(*GetGrid());
		const int is = bounds.is, ie = bounds.ie;
		const int js = bounds.js, je = bounds.je;
		const int ks = bounds.ks, ke = bounds.ke;
		const dimension_type dim = grid->GetDim();
		const double* volume_fraction = GetVolumeFraction();
		const double* pressure = data_manager->GetPrim(ID_PRESSURE);
		double inv_dx[3] = { 1.0 / GetDx(), 1.0 / GetDy(), 1.0 / GetDz() };

		for (int dir = 0; dir < static_cast<int>(dim); ++dir)
		{
			double* res_mom = data_manager->GetResidual(kMomentumOffset + dir);
#ifdef USE_OMP
#pragma omp parallel for collapse(3)
#endif
			for (int k = ks; k <= ke; ++k)
			{
				for (int j = js; j <= je; ++j)
				{
					for (int i = is; i <= ie; ++i)
					{
						const index_type idx = idx_proxy(i, j, k);
						int c_p[3] = { i, j, k };
						int c_m[3] = { i, j, k };
						c_p[dir] += 1;
						c_m[dir] -= 1;
						const double eps_p = volume_fraction[idx_proxy(c_p[0], c_p[1], c_p[2])];
						const double eps_m = volume_fraction[idx_proxy(c_m[0], c_m[1], c_m[2])];
						const double dedx = (eps_p - eps_m) * 0.5 * inv_dx[dir];
						// 注意：eps_p == eps_m 时上式精确等于 0.0（减出 0，乘 0 仍是 0）
						res_mom[idx] += pressure[idx] * dedx;
					}
				}
			}
		}
	}

	void EulerTwoPhaseStructUniform::CalcConvectionResidual()
	{
		// 基类的 AddFluxResidual 已经通过 GetVolumeFractionField() 做了 ε 加权
		EulerSolverStructUniform::CalcConvectionResidual();
		AddPorosityGradientSource();
	}

	/// @brief 相间反作用力源项（DEM-CFD **双向耦合**）
	/// @details 离散的 ε 加权动量方程是
	///              ∂(ε ρ u)/∂t + ∇·(ε ρ u u + ε p I) = p∇ε + S
	///          这里加的 S 就是"单位**混合物**体积上、流体受到的来自颗粒的力"。
	///          它与连续方程里 −ε_g∇p 的分工（重要，写错就是双计）：
	///              - 颗粒受的**压力梯度力** −V_p∇p 是 −ε_g∇p 的对偶项，
	///                它已经隐含在气相的 −ε_g∇p 里了（两者之和 = −∇p，正是混合物的压强项）；
	///              - 颗粒受的**曳力**在气相方程里没有任何对应项，必须由 S 补上。
	///          所以 S = −Σ_p F_drag,p / V_cell，**不包含** −V_p∇p。
	///          （若把 porosity_gradient_force 关掉、改用完全守恒的 ∇(εp) 形式，
	///            那就必须把压力梯度力也反作用回去 —— 耦合器会按该开关自动决定并告警。）
	void EulerTwoPhaseStructUniform::CalcSourceResidual()
	{
		// 基类目前是空实现；保留调用以便以后加真正的体积源（如重力）
		EulerSolverStructUniform::CalcSourceResidual();

		const double* src = m_interphase_source;
		if (src == nullptr)
		{
			return;
		}
		auto data_manager = GetDataManager();
		IdProxyStruct& idx_proxy = GetIdxProxy();
		const LoopBounds bounds(*GetGrid());
		const int is = bounds.is, ie = bounds.ie;
		const int js = bounds.js, je = bounds.je;
		const int ks = bounds.ks, ke = bounds.ke;
		// 相间力对气相做功能：∂(εE)/∂t + ... = S·u（可以关，见 SetInterphaseWork 的说明）
		double* res_energy = m_interphase_work ? data_manager->GetResidual(kEnergyOffset) : nullptr;
		// 做功项用的速度取与原变量同一时刻（S 是滞后量，u 用**步初**的节点速度），
		// 因此整个 RK 阶段这一项是常量，Σb_i = 1 ⇒ 一步的总功 = Σ_cells S·u·V·Δt，可精确核算。
		const double* vel[3] = {
			data_manager->GetPrim(ID_VELOCITY_X),
			data_manager->GetPrim(ID_VELOCITY_Y),
			data_manager->GetPrim(ID_VELOCITY_Z) };
		for (int dir = 0; dir < 3; ++dir)
		{
			double* res_mom = data_manager->GetResidual(kMomentumOffset + dir);
#ifdef USE_OMP
#pragma omp parallel for collapse(3)
#endif
			for (int k = ks; k <= ke; ++k)
			{
				for (int j = js; j <= je; ++j)
				{
					for (int i = is; i <= ie; ++i)
					{
						const index_type idx = idx_proxy(i, j, k);
						const double s = src[3 * idx + dir];
						res_mom[idx] += s;
						if (res_energy != nullptr)
						{
							res_energy[idx] += s * vel[dir][idx];
						}
					}
				}
			}
		}
	}

	// ==================================================================
	// 物理性检查
	// ==================================================================
	void EulerTwoPhaseStructUniform::CheckPrimtive()
	{
		EulerSolverStructUniform::CheckPrimtive();

		auto grid = GetGrid();
		IdProxyStruct& idx_proxy = GetIdxProxy();
		const LoopBounds bounds(*GetGrid());
		const int is = bounds.is, ie = bounds.ie;
		const int js = bounds.js, je = bounds.je;
		const int ks = bounds.ks, ke = bounds.ke;
		const double* volume_fraction = GetVolumeFraction();
		const double eps_min_allowed = m_tp_para->GetVolumeFractionMin();

		int bad_count = 0;
		double eps_min = LARGE_NUMBER, eps_max = -LARGE_NUMBER;
#ifdef USE_OMP
#pragma omp parallel for collapse(3) \
	reduction(+ : bad_count) reduction(min : eps_min) reduction(max : eps_max)
#endif
		for (int k = ks; k <= ke; ++k)
		{
			for (int j = js; j <= je; ++j)
			{
				for (int i = is; i <= ie; ++i)
				{
					const double eps = volume_fraction[idx_proxy(i, j, k)];
					eps_min = Min(eps_min, eps);
					eps_max = Max(eps_max, eps);
					if (!(eps >= eps_min_allowed) || eps > 1.0)
					{
						++bad_count;
					}
				}
			}
		}
		if (bad_count > 0)
		{
			Log::warn("EulerTwoPhaseStructUniform: {} nodes have volume fraction outside "
				"[{:E}, 1.0]. eps in [{:E}, {:E}]",
				bad_count, eps_min_allowed, eps_min, eps_max);
		}
	}
}
