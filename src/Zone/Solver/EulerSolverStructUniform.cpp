/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file EulerSolverStructUniform.cpp
 * \brief Euler equations solver on a uniform Cartesian structured grid.
 * \author Chen Jie.
 */

#include "EulerSolverStructUniform.h"
#include "Log.h"
#include "ZaranError.h"
#include "Limiter.h"
#include "MathBasic.h"
#include "GlobalData.h"

#include <cmath>
#include <vector>

namespace zaran
{
	// 方程个数（rho, rho*u, rho*v, rho*w, E）
	namespace
	{
		constexpr int kEqNum = 5;

		/// @brief 不含 ghost 的计算范围，统一转换成 int。
		/// @details OpenMP 的循环边界必须是简单变量：MSVC 对循环条件里直接出现
		///          static_cast<int>(...) 会报内部编译器错误 C1001。因此先把范围
		///          取成普通 int 局部变量，再用于 #pragma omp parallel for。
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

		/// @brief 全节点范围（含 ghost），同样转换成 int 供 OpenMP 使用
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
	EulerSolverStructUniform::EulerSolverStructUniform(index_type index, string name,
		shared_ptr<FlowSolverParam> para, shared_ptr<GridBase> grid,
		shared_ptr<DataManagerNS> data_manager)
		: NSSolver(index, name, para, grid, data_manager)
	{
		m_grid = std::dynamic_pointer_cast<GridStruct>(grid);
		if (!m_grid)
		{
			throw ZaranError("EulerSolverStructUniform requires a structured grid (GridStruct)");
		}
		m_data_manager = std::dynamic_pointer_cast<DataManagerNSStruct>(data_manager);
		if (!m_data_manager)
		{
			throw ZaranError("EulerSolverStructUniform requires DataManagerNSStruct");
		}
		m_para = std::dynamic_pointer_cast<FlowSolverParamUniform>(para);
		if (!m_para)
		{
			throw ZaranError("EulerSolverStructUniform requires FlowSolverParamUniform");
		}
		m_idx_proxy = make_shared<IdProxyStruct>(m_grid->GetNi(), m_grid->GetNj(), m_grid->GetNk());
	}

	EulerSolverStructUniform::~EulerSolverStructUniform()
	{
	}

	DataManagerNSStruct* EulerSolverStructUniform::GetDataManager()
	{
		return m_data_manager.get();
	}

	IdProxyStruct& EulerSolverStructUniform::GetIdxProxy()
	{
		return *m_idx_proxy;
	}

	FlowSolverParamUniform* EulerSolverStructUniform::GetPara()
	{
		return m_para.get();
	}

	GridStruct* EulerSolverStructUniform::GetGrid()
	{
		return m_grid.get();
	}

	// ==================================================================
	// 初值
	// ==================================================================
	void EulerSolverStructUniform::InitField()
	{
		auto init_type = GetPara()->GetInitFieldType();
		if (init_type == InitFieldType::Vortex)
		{
			InitFieldVortex();
			Prim2Cons();
			Log::info("Flow Field Initialize Finished! (Isentropic Vortex, rest frame)");
		}
		else if (init_type == InitFieldType::Riemann1D)
		{
			InitFieldRiemann1D();
			Prim2Cons();
			Log::info("Flow Field Initialize Finished! (1-D Riemann / shock tube)");
		}
		else
		{
			// 其余初值类型（远场 / 无速度远场 / 备份 / 爆炸）由基类统一处理
			NSSolver::InitField();
		}
	}

	void EulerSolverStructUniform::InitFieldFarfield()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto data_manager = GetDataManager();
		IdProxyStruct& idx_proxy = GetIdxProxy();
		auto ni = grid->GetNi();
		auto nj = grid->GetNj();
		auto nk = grid->GetNk();
		auto para = GetPara();

		double prim_far[kEqNum];
		prim_far[0] = para->GetInflowDensity();
		prim_far[1] = para->GetInflowVelocityX();
		prim_far[2] = para->GetInflowVelocityY();
		prim_far[3] = para->GetInflowVelocityZ();
		prim_far[4] = para->GetInflowPressure();

		for (index_type k = 0; k < nk; ++k)
		{
			for (index_type j = 0; j < nj; ++j)
			{
				for (index_type i = 0; i < ni; ++i)
				{
					// 远场初值允许沿 x 方向给出线性分布（便于做精度验证）
					double prim_local[kEqNum];
					for (int eq = 0; eq < kEqNum; ++eq)
					{
						prim_local[eq] = prim_far[eq];
					}
					data_manager->SetPrim(idx_proxy(i, j, k), prim_local);
				}
			}
		}
	}

	void EulerSolverStructUniform::InitFieldFarFieldZeroVel()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto data_manager = GetDataManager();
		IdProxyStruct& idx_proxy = GetIdxProxy();
		auto ni = grid->GetNi();
		auto nj = grid->GetNj();
		auto nk = grid->GetNk();
		auto para = GetPara();

		double prim[kEqNum];
		prim[0] = para->GetInflowDensity();
		prim[1] = 0.0;
		prim[2] = 0.0;
		prim[3] = 0.0;
		prim[4] = para->GetInflowPressure();

		for (index_type k = 0; k < nk; ++k)
		{
			for (index_type j = 0; j < nj; ++j)
			{
				for (index_type i = 0; i < ni; ++i)
				{
					data_manager->SetPrim(idx_proxy(i, j, k), prim);
				}
			}
		}
	}

	void EulerSolverStructUniform::InitFieldBackup()
	{
		// 与 NSSolverStruct 保持一致：续算读盘留待后续实现
		Log::warn("EulerSolverStructUniform: Backup initial field is not implemented yet, "
			"falling back to far-field.");
		InitFieldFarfield();
	}

	void EulerSolverStructUniform::InitFieldExplosion()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto data_manager = GetDataManager();
		IdProxyStruct& idx_proxy = GetIdxProxy();
		auto ni = grid->GetNi();
		auto nj = grid->GetNj();
		auto nk = grid->GetNk();
		auto para = GetPara();

		double prim_far[kEqNum];
		prim_far[0] = 1.0;
		prim_far[1] = 0.0;
		prim_far[2] = 0.0;
		prim_far[3] = 0.0;
		prim_far[4] = para->GetInflowPressure();

		const double pressure_ratio = GlobalData::GetDouble("init.explosion.pressure");
		double prim_inner[kEqNum];
		prim_inner[0] = 1.0;
		prim_inner[1] = 0.0;
		prim_inner[2] = 0.0;
		prim_inner[3] = 0.0;
		prim_inner[4] = pressure_ratio * prim_far[4];

		const double center[3] = {
			GlobalData::GetDouble("init.explosion.center_x"),
			GlobalData::GetDouble("init.explosion.center_y"),
			GlobalData::GetDouble("init.explosion.center_z")
		};
		const double radius = GlobalData::GetDouble("init.explosion.radius");
		Log::info("Explosion pressure ratio: {}, radius: {}", pressure_ratio, radius);

		for (index_type k = 0; k < nk; ++k)
		{
			for (index_type j = 0; j < nj; ++j)
			{
				for (index_type i = 0; i < ni; ++i)
				{
					const auto coord = node->GetCoord(i, j, k);
					const double dist = std::sqrt(
						(coord[0] - center[0]) * (coord[0] - center[0]) +
						(coord[1] - center[1]) * (coord[1] - center[1]) +
						(coord[2] - center[2]) * (coord[2] - center[2]));
					data_manager->SetPrim(idx_proxy(i, j, k),
						(dist <= radius) ? prim_inner : prim_far);
				}
			}
		}
	}

	/// @brief 二维等熵涡（Hu & Shu 1999），涡心 (5,5)、强度 beta=5、远场静止。
	/// @details 与 NSSolverStruct::InitFieldVortex 完全一致，便于两种求解器互相印证。
	///          注意取远场静止（u=0）而不是标准版的远场 u=1：周期边界下静止涡不会
	///          靠近边界，因而参考解里不需要叠加周期像（详见 docs 中的说明）。
	void EulerSolverStructUniform::InitFieldVortex()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto data_manager = GetDataManager();
		IdProxyStruct& idx_proxy = GetIdxProxy();
		auto ni = grid->GetNi();
		auto nj = grid->GetNj();
		auto nk = grid->GetNk();

		const double beta = 5.0;
		const double gamma = GetGas()->GetGamma();
		const double xc = 5.0;
		const double yc = 5.0;
		double prim[kEqNum];

		for (index_type k = 0; k < nk; ++k)
		{
			for (index_type j = 0; j < nj; ++j)
			{
				for (index_type i = 0; i < ni; ++i)
				{
					const auto coord = node->GetCoord(i, j, k);
					const double x = coord[0];
					const double y = coord[1];
					const double r2 = (x - xc) * (x - xc) + (y - yc) * (y - yc);
					prim[0] = std::pow(1.0 - (gamma - 1.0) * beta * beta * std::exp(1.0 - r2)
						/ (8.0 * gamma * PI * PI), 1.0 / (gamma - 1.0));
					prim[4] = std::pow(prim[0], gamma);
					prim[1] = beta * std::exp(0.5 * (1.0 - r2)) / (2.0 * PI) * (-y + yc);
					prim[2] = beta * std::exp(0.5 * (1.0 - r2)) / (2.0 * PI) * (x - xc);
					prim[3] = 0.0;
					data_manager->SetPrim(idx_proxy(i, j, k), prim);
				}
			}
		}
	}

	/// @brief 一维 Riemann（激波管）初值：按 x 方向膜片位置分两段常值。
	/// @details 相关参数在控制文件的 [init.riemann] 节给出。
	void EulerSolverStructUniform::InitFieldRiemann1D()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto data_manager = GetDataManager();
		IdProxyStruct& idx_proxy = GetIdxProxy();
		auto ni = grid->GetNi();
		auto nj = grid->GetNj();
		auto nk = grid->GetNk();

		const double x_diaphragm = GlobalData::GetDouble("init.riemann.x_diaphragm");
		double prim_left[kEqNum], prim_right[kEqNum];
		prim_left[0] = GlobalData::GetDouble("init.riemann.left_density");
		prim_left[1] = GlobalData::GetDouble("init.riemann.left_velocity_x");
		prim_left[2] = GlobalData::GetDouble("init.riemann.left_velocity_y");
		prim_left[3] = GlobalData::GetDouble("init.riemann.left_velocity_z");
		prim_left[4] = GlobalData::GetDouble("init.riemann.left_pressure");
		prim_right[0] = GlobalData::GetDouble("init.riemann.right_density");
		prim_right[1] = GlobalData::GetDouble("init.riemann.right_velocity_x");
		prim_right[2] = GlobalData::GetDouble("init.riemann.right_velocity_y");
		prim_right[3] = GlobalData::GetDouble("init.riemann.right_velocity_z");
		prim_right[4] = GlobalData::GetDouble("init.riemann.right_pressure");

		Log::info("Riemann 1-D initial field: diaphragm x = {}, left = (rho {}, u {}, p {}), right = (rho {}, u {}, p {})",
			x_diaphragm, prim_left[0], prim_left[1], prim_left[4],
			prim_right[0], prim_right[1], prim_right[4]);

		for (index_type k = 0; k < nk; ++k)
		{
			for (index_type j = 0; j < nj; ++j)
			{
				for (index_type i = 0; i < ni; ++i)
				{
					const auto coord = node->GetCoord(i, j, k);
					const bool is_left = (coord[0] <= x_diaphragm);
					data_manager->SetPrim(idx_proxy(i, j, k), is_left ? prim_left : prim_right);
				}
			}
		}
	}

	// ==================================================================
	// 网格：均匀性校验与步长记录
	// ==================================================================
	void EulerSolverStructUniform::CalcCoordTransCoef()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		const LoopBounds bounds(*m_grid);
		const int is = bounds.is, ie = bounds.ie;
		const int js = bounds.js, je = bounds.je;
		const int ks = bounds.ks, ke = bounds.ke;
		const dimension_type dim = grid->GetDim();

		const auto origin = node->GetCoord(is, js, ks);
		const auto next_x = node->GetCoord(is + 1, js, ks);
		m_dx = next_x[0] - origin[0];
		m_dy = 1.0;
		m_dz = 1.0;
		if (dim >= TWO_DIM)
		{
			const auto next_y = node->GetCoord(is, js + 1, ks);
			m_dy = next_y[1] - origin[1];
		}
		if (dim >= THREE_DIM)
		{
			const auto next_z = node->GetCoord(is, js, ks + 1);
			m_dz = next_z[2] - origin[2];
		}

		if (!(m_dx > 0.0) || !(m_dy > 0.0) || !(m_dz > 0.0))
		{
			Log::error("Euler uniform solver: non-positive grid spacing dx={}, dy={}, dz={}",
				m_dx, m_dy, m_dz);
			throw ZaranError("EulerSolverStructUniform: invalid grid spacing");
		}

		// 校验网格是否严格均匀（本求解器的前提条件）
		double max_dev = 0.0;
		for (index_type k = ks; k <= ke; ++k)
		{
			for (index_type j = js; j <= je; ++j)
			{
				for (index_type i = is; i < ie; ++i)
				{
					const auto a = node->GetCoord(i, j, k);
					const auto b = node->GetCoord(i + 1, j, k);
					max_dev = Max(max_dev, std::fabs((b[0] - a[0]) - m_dx) / m_dx);
				}
			}
		}
		if (dim >= TWO_DIM)
		{
			for (index_type k = ks; k <= ke; ++k)
			{
				for (index_type j = js; j < je; ++j)
				{
					for (index_type i = is; i <= ie; ++i)
					{
						const auto a = node->GetCoord(i, j, k);
						const auto b = node->GetCoord(i, j + 1, k);
						max_dev = Max(max_dev, std::fabs((b[1] - a[1]) - m_dy) / m_dy);
					}
				}
			}
		}
		if (dim >= THREE_DIM)
		{
			for (index_type k = ks; k < ke; ++k)
			{
				for (index_type j = js; j <= je; ++j)
				{
					for (index_type i = is; i <= ie; ++i)
					{
						const auto a = node->GetCoord(i, j, k);
						const auto b = node->GetCoord(i, j, k + 1);
						max_dev = Max(max_dev, std::fabs((b[2] - a[2]) - m_dz) / m_dz);
					}
				}
			}
		}

		Log::info("Euler uniform grid: dim={}, nodes=({},{},{}), dx={:E}, dy={:E}, dz={:E}, max spacing deviation={:E}",
			static_cast<int>(dim), grid->GetNi(), grid->GetNj(), grid->GetNk(),
			m_dx, m_dy, m_dz, max_dev);

		if (m_para->GetRequireUniformGrid() && max_dev > 1.0e-6)
		{
			Log::error("Euler uniform solver requires a strictly uniform grid, "
				"but the maximum relative spacing deviation is {:E}.", max_dev);
			throw ZaranError("EulerSolverStructUniform: grid is not uniform");
		}
	}

	// ==================================================================
	// 变量转换与残差清零
	// ==================================================================
	void EulerSolverStructUniform::Prim2Cons()
	{
		auto gas = GetGas();
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		IdProxyStruct& idx_proxy = GetIdxProxy();
		const NodeCounts counts(*m_grid);
		const int ni = counts.ni, nj = counts.nj, nk = counts.nk;
		// OpenMP 要求循环边界是简单变量
		const int ni_i = ni;
		const int nj_i = nj;
		const int nk_i = nk;
		// 用字面量尺寸 + private 子句，与 NSSolverStruct::Prim2Cons 写法一致
		double prim[5], cons[5];

#ifdef USE_OMP
#pragma omp parallel for collapse(3) private(prim, cons)
#endif
		for (int k = 0; k < nk_i; ++k)
		{
			for (int j = 0; j < nj_i; ++j)
			{
				for (int i = 0; i < ni_i; ++i)
				{
					const index_type idx = idx_proxy(i, j, k);
					for (int eq = 0; eq < kEqNum; ++eq)
					{
						prim[eq] = data_manager->GetPrim(eq, idx);
					}
					gas->Prim2Cons(prim, cons);
					data_manager->SetCons(idx, cons);
					// 注意：DataManagerNS 只为 SetPrim / SetCons / SetResidual 提供了指针重载，
					// SetConsOld 需要逐分量写入
					for (int eq = 0; eq < kEqNum; ++eq)
					{
						data_manager->SetConsOld(eq, idx, cons[eq]);
					}
				}
			}
		}
	}

	void EulerSolverStructUniform::Cons2Prim()
	{
		auto gas = GetGas();
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		IdProxyStruct& idx_proxy = GetIdxProxy();
		const NodeCounts counts(*m_grid);
		const int ni = counts.ni, nj = counts.nj, nk = counts.nk;
		const int ni_i = ni;
		const int nj_i = nj;
		const int nk_i = nk;
		double prim[5], cons[5];

#ifdef USE_OMP
#pragma omp parallel for collapse(3) private(prim, cons)
#endif
		for (int k = 0; k < nk_i; ++k)
		{
			for (int j = 0; j < nj_i; ++j)
			{
				for (int i = 0; i < ni_i; ++i)
				{
					const index_type idx = idx_proxy(i, j, k);
					for (int eq = 0; eq < kEqNum; ++eq)
					{
						cons[eq] = data_manager->GetCons(eq, idx);
					}
					gas->Cons2Prim(cons, prim);
					data_manager->SetPrim(idx, prim);
				}
			}
		}
	}

	void EulerSolverStructUniform::ZeroResidual()
	{
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		IdProxyStruct& idx_proxy = GetIdxProxy();
		const NodeCounts counts(*m_grid);
		const int ni = counts.ni, nj = counts.nj, nk = counts.nk;
		const double res_zero[kEqNum] = { 0.0, 0.0, 0.0, 0.0, 0.0 };

#ifdef USE_OMP
#pragma omp parallel for collapse(3)
#endif
		for (int k = 0; k < nk; ++k)
		{
			for (int j = 0; j < nj; ++j)
			{
				for (int i = 0; i < ni; ++i)
				{
					data_manager->SetResidual(idx_proxy(i, j, k), res_zero);
				}
			}
		}
	}

	// ==================================================================
	// 空间离散：界面通量与对流残差
	// ==================================================================
	double EulerSolverStructUniform::ApplyLimiter(double slope_left, double slope_right)
	{
		switch (GetPara()->GetLimiterType())
		{
		case LimiterType::none:
			// 不做限制：中心差分斜率（会产生振荡，仅用于对照）
			return 0.5 * (slope_left + slope_right);
		case LimiterType::first_order:
			return 0.0;
		case LimiterType::barth:
			// 框架中的 barth 是面向非结构网格的几何限制器；
			// 在均匀网格的 MUSCL 重构中对应更保守的 minmod。
			return LimiterMinMod(slope_left, slope_right);
		case LimiterType::vk:
			// Venkatakrishnan 限制器追求"尽量接近二阶"，在 MUSCL 的斜率限制
			// 语境下用同样平滑的 VanLeer 限制器对应（这也是默认配置）。
			return LimiterVanLeer(slope_left, slope_right);
		default:
			break;
		}
		return LimiterVanLeer(slope_left, slope_right);
	}

	void EulerSolverStructUniform::ReconstructFace(const double* const prim[kEqNum],
		index_type idx_m, index_type idx_0, index_type idx_p, index_type idx_pp,
		bool first_order, double* prim_left, double* prim_right)
	{
		if (first_order)
		{
			for (int eq = 0; eq < kEqNum; ++eq)
			{
				prim_left[eq] = prim[eq][idx_0];
				prim_right[eq] = prim[eq][idx_p];
			}
			return;
		}
		// MUSCL：slope 表示"斜率 x dx"，即相邻节点值之差经限制器作用后的结果
		for (int eq = 0; eq < kEqNum; ++eq)
		{
			const double slope_0 = ApplyLimiter(prim[eq][idx_0] - prim[eq][idx_m],
				prim[eq][idx_p] - prim[eq][idx_0]);
			const double slope_p = ApplyLimiter(prim[eq][idx_p] - prim[eq][idx_0],
				prim[eq][idx_pp] - prim[eq][idx_p]);
			prim_left[eq] = prim[eq][idx_0] + 0.5 * slope_0;
			prim_right[eq] = prim[eq][idx_p] - 0.5 * slope_p;
		}
	}

	void EulerSolverStructUniform::AddFluxResidual(int dir, bool first_order)
	{
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		IdProxyStruct& idx_proxy = GetIdxProxy();
		const LoopBounds bounds(*m_grid);
		const int is = bounds.is, ie = bounds.ie;
		const int js = bounds.js, je = bounds.je;
		const int ks = bounds.ks, ke = bounds.ke;

		const double gamma = GetGas()->GetGamma();
		const double dx = (dir == 0) ? m_dx : ((dir == 1) ? m_dy : m_dz);
		const double inv_dx = 1.0 / dx;

		const double* prim[kEqNum];
		for (int eq = 0; eq < kEqNum; ++eq)
		{
			prim[eq] = data_manager->GetPrim(eq);
		}
		double* residual[kEqNum];
		for (int eq = 0; eq < kEqNum; ++eq)
		{
			residual[eq] = data_manager->GetResidual(eq);
		}

		Eigen::Vector3d norm = Eigen::Vector3d::Zero();
		norm(dir) = 1.0;

		// 物理节点范围（不含 ghost）
		index_type lo[3] = { is, js, ks };
		index_type hi[3] = { ie, je, ke };

		// 沿 dir 的面 a+1/2（a 从 is-1 到 ie，使区间 [is,ie] 的左右面都被覆盖）
		index_type face_lo[3] = { is, js, ks };
		face_lo[dir] = lo[dir] - 1;

		// ---- 阶段 1：逐面重构 + Riemann 通量，写入半点通量数组 ----
		for (index_type k = face_lo[2]; k <= hi[2]; ++k)
		{
			for (index_type j = face_lo[1]; j <= hi[1]; ++j)
			{
				for (index_type i = face_lo[0]; i <= hi[0]; ++i)
				{
					index_type c_l[3] = { i, j, k };
					index_type c_r[3] = { i, j, k };
					c_r[dir] += 1;
					const index_type idx_l = idx_proxy(c_l[0], c_l[1], c_l[2]);
					const index_type idx_r = idx_proxy(c_r[0], c_r[1], c_r[2]);

					double prim_left[kEqNum], prim_right[kEqNum];
					if (first_order)
					{
						ReconstructFace(prim, idx_l, idx_l, idx_r, idx_r, true,
							prim_left, prim_right);
					}
					else
					{
						index_type c_m[3] = { i, j, k };
						c_m[dir] -= 1;
						index_type c_pp[3] = { i, j, k };
						c_pp[dir] += 2;
						ReconstructFace(prim,
							idx_proxy(c_m[0], c_m[1], c_m[2]), idx_l, idx_r,
							idx_proxy(c_pp[0], c_pp[1], c_pp[2]), false,
							prim_left, prim_right);
					}

					RiemannSolverPara riemann_para;
					for (int eq = 0; eq < kEqNum; ++eq)
					{
						riemann_para.prim_left(eq) = prim_left[eq];
						riemann_para.prim_right(eq) = prim_right[eq];
					}
					riemann_para.gamma_left = gamma;
					riemann_para.gamma_right = gamma;
					// 注意：HLLC 等求解器会原地修改 norm/nt，因此每次都必须重新设置
					riemann_para.norm = norm;
					riemann_para.nt = 0.0;
					m_riemann_solver->Solver(riemann_para);

					for (int eq = 0; eq < kEqNum; ++eq)
					{
						data_manager->SetMidNodeFlux(eq, dir, idx_l, riemann_para.flux(eq));
					}
				}
			}
		}

		// ---- 阶段 2：通量差分，累加到节点残差 ----
		// dU_i/dt = -( F_{i+1/2} - F_{i-1/2} ) / dx ，其中半点通量以左节点为索引
#ifdef USE_OMP
#pragma omp parallel for collapse(3)
#endif
		for (int k = ks; k <= ke; ++k)
		{
			for (int j = js; j <= je; ++j)
			{
				for (int i = is; i <= ie; ++i)
				{
					int c[3] = { i, j, k };
					int c_m[3] = { i, j, k };
					c_m[dir] -= 1;
					const index_type idx = idx_proxy(c[0], c[1], c[2]);
					const index_type idx_m = idx_proxy(c_m[0], c_m[1], c_m[2]);
					for (int eq = 0; eq < kEqNum; ++eq)
					{
						const double flux_right = data_manager->GetMidNodeFlux(eq, dir, idx);
						const double flux_left = data_manager->GetMidNodeFlux(eq, dir, idx_m);
						residual[eq][idx] -= (flux_right - flux_left) * inv_dx;
					}
				}
			}
		}
	}

	void EulerSolverStructUniform::CalcConvectionResidual()
	{
		const bool first_order = UseFirstOrder();
		const dimension_type dim = m_grid->GetDim();
		AddFluxResidual(0, first_order);
		if (dim >= TWO_DIM)
		{
			AddFluxResidual(1, first_order);
		}
		if (dim >= THREE_DIM)
		{
			AddFluxResidual(2, first_order);
		}
	}

	// ==================================================================
	// 时间推进
	// ==================================================================
	bool EulerSolverStructUniform::UseFirstOrder()
	{
		auto para = GetPara();
		if (para->GetReconstructionOrder() <= 1)
		{
			return true;
		}
		const int first_order_steps = para->GetFirstOrderSteps();
		if (first_order_steps > 0)
		{
			const int current_iter = GlobalData::IsExist("iteration.current_iter")
				? GlobalData::GetInt("iteration.current_iter") : 0;
			return current_iter < first_order_steps;
		}
		return false;
	}

	void EulerSolverStructUniform::CalcTimeStepLocal()
	{
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		IdProxyStruct& idx_proxy = GetIdxProxy();
		const LoopBounds bounds(*m_grid);
		const int is = bounds.is, ie = bounds.ie;
		const int js = bounds.js, je = bounds.je;
		const int ks = bounds.ks, ke = bounds.ke;
		const dimension_type dim = grid->GetDim();
		const double cfl = GetPara()->GetCflNumber();
		const double gamma = GetGas()->GetGamma();
		const double inv_dx = 1.0 / m_dx;
		const double inv_dy = 1.0 / m_dy;
		const double inv_dz = 1.0 / m_dz;

		const double* density = data_manager->GetPrim(ID_DENSITY);
		const double* pressure = data_manager->GetPrim(ID_PRESSURE);
		const double* velocity_x = data_manager->GetPrim(ID_VELOCITY_X);
		const double* velocity_y = data_manager->GetPrim(ID_VELOCITY_Y);
		const double* velocity_z = data_manager->GetPrim(ID_VELOCITY_Z);

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
					const double c = std::sqrt(gamma * pressure[idx] / density[idx]);
					double lamda = (std::fabs(velocity_x[idx]) + c) * inv_dx;
					if (dim >= TWO_DIM)
					{
						lamda += (std::fabs(velocity_y[idx]) + c) * inv_dy;
					}
					if (dim >= THREE_DIM)
					{
						lamda += (std::fabs(velocity_z[idx]) + c) * inv_dz;
					}
					data_manager->SetTimeStep(idx, cfl / lamda);
				}
			}
		}
	}

	void EulerSolverStructUniform::CalcMinTimeStep(double& min_dt)
	{
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		IdProxyStruct& idx_proxy = GetIdxProxy();
		const LoopBounds bounds(*m_grid);
		const int is = bounds.is, ie = bounds.ie;
		const int js = bounds.js, je = bounds.je;
		const int ks = bounds.ks, ke = bounds.ke;
		double local_min = LARGE_NUMBER;

#ifdef USE_OMP
#pragma omp parallel for collapse(3) reduction(min : local_min)
#endif
		for (int k = ks; k <= ke; ++k)
		{
			for (int j = js; j <= je; ++j)
			{
				for (int i = is; i <= ie; ++i)
				{
					const index_type idx = idx_proxy(i, j, k);
					const double dt = data_manager->GetTimeStep(idx);
					if (dt < local_min)
					{
						local_min = dt;
					}
				}
			}
		}
		min_dt = local_min;
	}

	void EulerSolverStructUniform::ReduceTimeStep(double& dt)
	{
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		IdProxyStruct& idx_proxy = GetIdxProxy();
		const NodeCounts counts(*m_grid);
		const int ni = counts.ni, nj = counts.nj, nk = counts.nk;

#ifdef USE_OMP
#pragma omp parallel for collapse(3)
#endif
		for (int k = 0; k < nk; ++k)
		{
			for (int j = 0; j < nj; ++j)
			{
				for (int i = 0; i < ni; ++i)
				{
					data_manager->SetTimeStep(idx_proxy(i, j, k), dt);
				}
			}
		}
	}

	/// @brief SSP-RK 时间推进，与 NSSolverStruct::RungeKutta 使用同一套 RK 系数约定：
	///        cons = rk_coef[0] * cons_old + rk_coef[1] * cons + rk_coef[2] * dt * residual
	///        均匀网格的坐标雅可比恒为 1，故不出现 jacobi 因子。
	void EulerSolverStructUniform::RungeKutta()
	{
		auto grid = GetGrid();
		auto para = GetPara();
		auto data_manager = GetDataManager();
		IdProxyStruct& idx_proxy = GetIdxProxy();
		const LoopBounds bounds(*m_grid);
		const int is = bounds.is, ie = bounds.ie;
		const int js = bounds.js, je = bounds.je;
		const int ks = bounds.ks, ke = bounds.ke;
		const int rk_step = para->GetRkStep();

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
					for (int eq = 0; eq < kEqNum; ++eq)
					{
						data_manager->SetConsOld(eq, idx, data_manager->GetCons(eq, idx));
					}
				}
			}
		}

		for (int iStep = 0; iStep < rk_step; ++iStep)
		{
			BoundaryCondition();
			CalcResidual();
			const auto& rk_coef = para->GetRkCoef(iStep);

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
						const double dt = data_manager->GetTimeStep(idx);
						for (int eq = 0; eq < kEqNum; ++eq)
						{
							data_manager->SetCons(eq, idx,
								rk_coef[0] * data_manager->GetConsOld(eq, idx) +
								rk_coef[1] * data_manager->GetCons(eq, idx) +
								rk_coef[2] * dt * data_manager->GetResidual(eq, idx));
						}
					}
				}
			}
			// 更新原始变量，供下一个 RK 阶段（以及结果输出）使用
			UpdateField();
		}
	}

	// ==================================================================
	// 边界条件
	// ==================================================================
	/// @brief 按网格的边界映射填充虚拟节点。
	/// @details 支持的类型与 mesh.inp 中使用的名字一致：
	///          - "outlet"   零梯度外推（也是激波管的入流/出流边界）
	///          - "inlet"    给定来流常值
	///          - "wall"     固壁（法向动量反号）
	///          - "hole"     内部挖洞，跳过
	///          未识别的名字会直接报错，避免静默地使用错误边界。
	void EulerSolverStructUniform::BoundaryCondition()
	{
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		IdProxyStruct& idx_proxy = GetIdxProxy();
		auto bound_map = grid->GetBoundMap();
		const size_t ghost_level = grid->GetGhostLevel();
		const dimension_type dim = grid->GetDim();

		double prim_inflow[kEqNum] = {
			GetPara()->GetInflowDensity(),
			GetPara()->GetInflowVelocityX(),
			GetPara()->GetInflowVelocityY(),
			GetPara()->GetInflowVelocityZ(),
			GetPara()->GetInflowPressure()
		};
		double cons_inflow[kEqNum];
		GetGas()->Prim2Cons(prim_inflow, cons_inflow);

		for (auto& boundary : bound_map->GetBoundMap())
		{
			const auto& bound_name = boundary.first;
			auto& bound = boundary.second;
			if (bound_name == "hole")
			{
				continue;
			}
			if (bound_name != "outlet" && bound_name != "inlet" && bound_name != "wall")
			{
				Log::error("EulerSolverStructUniform: unsupported boundary type '{}'. "
					"Supported: outlet | inlet | wall.", bound_name);
				throw ZaranError("EulerSolverStructUniform: unsupported boundary type: " + bound_name);
			}

			for (size_t iBound = 0; iBound < bound.size(); ++iBound)
			{
				index_type i_bound, j_bound, k_bound;
				bound[iBound].GetIdx(i_bound, j_bound, k_bound);
				const index_type idx_bound = idx_proxy(i_bound, j_bound, k_bound);
				const auto direction = bound[iBound].GetDirectionSrc();

				// 边界节点自身的原始变量（outlet / wall 需要用它外推或镜像）
				double prim_bound[kEqNum], cons_bound[kEqNum];
				for (int eq = 0; eq < kEqNum; ++eq)
				{
					prim_bound[eq] = data_manager->GetPrim(eq, idx_bound);
				}
				GetGas()->Prim2Cons(prim_bound, cons_bound);

				for (size_t iGhost = 1; iGhost <= ghost_level; ++iGhost)
				{
					const index_type i_ghost = i_bound + iGhost * direction[0];
					const index_type j_ghost = j_bound + iGhost * direction[1];
					const index_type k_ghost = k_bound + iGhost * direction[2];
					const index_type idx_ghost = idx_proxy(i_ghost, j_ghost, k_ghost);

					if (bound_name == "inlet")
					{
						data_manager->SetPrim(idx_ghost, prim_inflow);
						data_manager->SetCons(idx_ghost, cons_inflow);
					}
					else if (bound_name == "outlet")
					{
						data_manager->SetPrim(idx_ghost, prim_bound);
						data_manager->SetCons(idx_ghost, cons_bound);
					}
					else // wall
					{
						// 固壁：密度、压力、切向动量取镜像，法向动量反号
						// 边界方向指向计算域外，故法向分量为 -direction
						dimension_type normal_dir = 0;
						double normal_sign = 1.0;
						for (dimension_type iDim = 0; iDim < dim; ++iDim)
						{
							if (direction[iDim] != 0)
							{
								normal_dir = iDim;
								normal_sign = -static_cast<double>(direction[iDim]);
								break;
							}
						}
						double cons_ghost[kEqNum];
						for (int eq = 0; eq < kEqNum; ++eq)
						{
							cons_ghost[eq] = cons_bound[eq];
						}
						// 指向域外 => ghost 的法向速度与边界节点相反
						cons_ghost[1 + normal_dir] = -cons_bound[1 + normal_dir];
						double prim_ghost[kEqNum];
						GetGas()->Cons2Prim(cons_ghost, prim_ghost);
						(void)normal_sign;
						data_manager->SetPrim(idx_ghost, prim_ghost);
						data_manager->SetCons(idx_ghost, cons_ghost);
					}
				}
			}
		}
	}

	// ==================================================================
	// 物理性检查与备份
	// ==================================================================
	void EulerSolverStructUniform::CheckPrimtive()
	{
		auto grid = GetGrid();
		auto data_manager = GetDataManager();
		IdProxyStruct& idx_proxy = GetIdxProxy();
		const LoopBounds bounds(*m_grid);
		const int is = bounds.is, ie = bounds.ie;
		const int js = bounds.js, je = bounds.je;
		const int ks = bounds.ks, ke = bounds.ke;

		int bad_count = 0;
		double rho_min = LARGE_NUMBER, p_min = LARGE_NUMBER;
		double rho_max = -LARGE_NUMBER, p_max = -LARGE_NUMBER;

#ifdef USE_OMP
#pragma omp parallel for collapse(3) \
	reduction(+ : bad_count) reduction(min : rho_min, p_min) reduction(max : rho_max, p_max)
#endif
		for (int k = ks; k <= ke; ++k)
		{
			for (int j = js; j <= je; ++j)
			{
				for (int i = is; i <= ie; ++i)
				{
					const index_type idx = idx_proxy(i, j, k);
					const double rho = data_manager->GetPrim(ID_DENSITY, idx);
					const double p = data_manager->GetPrim(ID_PRESSURE, idx);
					rho_min = Min(rho_min, rho);
					rho_max = Max(rho_max, rho);
					p_min = Min(p_min, p);
					p_max = Max(p_max, p);
					if (!(rho > 0.0) || !(p > 0.0))
					{
						++bad_count;
					}
				}
			}
		}

		if (bad_count > 0)
		{
			Log::warn("EulerSolverStructUniform: {} nodes have non-physical primitive values "
				"(rho or p <= 0). rho in [{:E}, {:E}], p in [{:E}, {:E}]",
				bad_count, rho_min, rho_max, p_min, p_max);
		}
	}

	void EulerSolverStructUniform::CheckResidual()
	{
		// 残差检查由 FieldNS::CalcResidual 统一完成（输出残差文件），此处无需重复
	}

	void EulerSolverStructUniform::FixPrimtive()
	{
		// 与 NSSolverStruct 保持一致：暂不做强制修正。
		// 若后续算例出现强膨胀导致负压强，可在此处按下限截断或重置为来流值。
	}

	void EulerSolverStructUniform::BackupField(std::string& back_folder)
	{
		// 与 NSSolverStruct 保持一致：流场落盘备份留待后续实现。
		// 输出通道（Tecplot / 残差文件）已由 FieldSimulation 统一处理。
	}
}
