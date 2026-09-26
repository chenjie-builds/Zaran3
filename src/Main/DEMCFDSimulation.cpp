/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file DEMCFDSimulation.cpp
 * \brief One-way coupled DEM-CFD case runner (`task.solver = "EulerDEM"`).
 */

#include "DEMCFDSimulation.h"
#include "NSFieldStructTwoPhase.h"
#include "EulerTwoPhaseStructUniform.h"
#include "FlowSolverParamTwoPhase.h"
#include "DataManagerNSTwoPhase.h"
#include "GlobalData.h"
#include "Log.h"
#include "ZaranError.h"
#include "Visual.h"
#include "File.h"
#include "CommonPara.h"
#include "DEMFieldData.h"
#include "ReadDEMParticle.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

namespace zaran
{
	namespace
	{
		double ReadDoubleOr(const char* key, double fallback)
		{
			return GlobalData::IsExist(key) ? GlobalData::GetDouble(key) : fallback;
		}

		int ReadIntOr(const char* key, int fallback)
		{
			return GlobalData::IsExist(key) ? GlobalData::GetInt(key) : fallback;
		}

		bool ReadBoolOr(const char* key, bool fallback)
		{
			return GlobalData::IsExist(key) ? GlobalData::GetBool(key) : fallback;
		}

		std::string ReadStringOr(const char* key, const std::string& fallback)
		{
			return GlobalData::IsExist(key) ? GlobalData::GetString(key) : fallback;
		}

		/// @brief 整块落盘（与 MappingSimulation 同样的写法）
		bool WriteText(const std::string& path, const std::string& content)
		{
			std::ofstream out(path);
			if (!out.is_open())
			{
				Log::warn("DEMCFDSimulation: cannot open '{}' for writing", path);
				return false;
			}
			out.write(content.data(), static_cast<std::streamsize>(content.size()));
			return out.good();
		}

		/// @brief 粒子力历史 CSV 的表头（列顺序与 WriteForceHistory 严格一致）
		/// @details 前 36 列是阶段 3 的列（列名不变，便于旧脚本沿用）；
		///          后面 11 列是双向耦合新增的颗粒运动/反作用量。
		///          ⚠ 改这里就必须同步改 WriteForceHistory 的 printf 说明符个数 ——
		///            少/多一个不会报错，只会让后面的列整体错位（verify.py 的 V1 会抓）。
		const char* kForceHeader =
			"iter,time,id,group,x,y,z,radius,density,eps_gas,weight_sum,cells,"
			"rho_g,u_g,v_g,w_g,p_g,T_g,mu_g,gradpx,gradpy,gradpz,Re,Cd,branch,beta_over_es,"
			"Fx_drag,Fy_drag,Fz_drag,Fx_pg,Fy_pg,Fz_pg,Fx,Fy,Fz,F_mag,"
			"mass,u_p,v_p,w_p,a_px,a_py,a_pz,n_sub,Freact_x,Freact_y,Freact_z\n";
	}

	// ==================================================================
	// 时间推进：完全复用基类，只在收尾补一份报告
	// ==================================================================
	void DEMCFDSimulation::SolveField()
	{
		// 基类负责 Initialize()（前会调用 PrepareCoupling()）、时间循环与逐帧输出；
		// 本类只补一件事：循环结束后把关键量写成 coupling_report.csv。
		NSFieldSimulation::SolveField();
		WriteReport();
	}

	DEMCFDSimulation::DEMCFDSimulation(const shared_ptr<FieldManager>& field_manager)
		: NSFieldSimulation(field_manager)
	{
	}

	DEMCFDSimulation::~DEMCFDSimulation()
	{
	}

	// ==================================================================
	// 参数与输入
	// ==================================================================
	void DEMCFDSimulation::ReadParameters()
	{
		m_particle_file = ReadStringOr("coupling.particle_file", "particles.csv");
		m_drag_model_name = ReadStringOr("coupling.drag_model", "gidaspow");
		m_mapping_method = ReadStringOr("coupling.mapping_method", "subcell");
		m_mapping_sub = ReadIntOr("coupling.mapping_sub", 2);
		m_smooth_passes = ReadIntOr("coupling.smooth_passes", 0);
		m_smooth_blend = ReadDoubleOr("coupling.smooth_blend", 1.0);
		m_volume_fraction_from_particles =
			ReadBoolOr("coupling.volume_fraction_from_particles", true);
		m_freeze_particles = ReadBoolOr("coupling.freeze_particles", true);
		m_particle_density_fallback = ReadDoubleOr("coupling.particle_density", 1000.0);
		m_eps_min_clip = ReadDoubleOr("coupling.eps_min_clip", 1.0e-3);
		m_force_write_interval = ReadIntOr("coupling.force_write_interval", 1);
		m_write_coupling_frames = ReadBoolOr("coupling.write_coupling_frames", true);

		m_options.drag.dense_threshold = ReadDoubleOr("coupling.dense_threshold", 0.8);
		m_options.drag.dilute_re_uses_voidage =
			ReadBoolOr("coupling.dilute_re_uses_voidage", true);
		m_options.drag.include_pressure_gradient_force =
			ReadBoolOr("coupling.include_pressure_gradient_force", true);
		m_options.drag.voidage_exponent = ReadDoubleOr("coupling.voidage_exponent", 2.65);
		m_options.volume_fraction_from_particles = m_volume_fraction_from_particles;
		m_options.eps_min_clip = m_eps_min_clip;
		m_options.smooth_passes = m_smooth_passes;
		m_options.smooth_blend = m_smooth_blend;
		m_options.mapping_sub = m_mapping_sub;

		// ---- 双向耦合（阶段 4）----
		m_particle_motion_name = ReadStringOr("coupling.particle_motion", "fixed");
		m_integrator_name = ReadStringOr("coupling.particle_integrator", "exponential");
		m_react_on_fluid = ReadBoolOr("coupling.react_on_fluid", false);
		m_react_work = ReadBoolOr("coupling.react_interphase_work", true);
		m_react_pg_mode = ReadStringOr("coupling.react_pressure_gradient_force", "auto");
		m_particle_substeps = ReadIntOr("coupling.particle_substeps", 0);
		m_particle_max_step_param = ReadDoubleOr("coupling.particle_max_step_parameter", 0.25);
		m_body_acceleration[0] = ReadDoubleOr("coupling.body_acceleration_x", 0.0);
		m_body_acceleration[1] = ReadDoubleOr("coupling.body_acceleration_y", 0.0);
		m_body_acceleration[2] = ReadDoubleOr("coupling.body_acceleration_z", 0.0);
		m_remap_each_step = ReadBoolOr("coupling.remap_volume_fraction", true);

		if (m_particle_motion_name == "fixed")
		{
			m_motion.mode = DEMCFDCoupler::MotionMode::Fixed;
		}
		else if (m_particle_motion_name == "prescribed")
		{
			m_motion.mode = DEMCFDCoupler::MotionMode::Prescribed;
		}
		else if (m_particle_motion_name == "dynamic")
		{
			m_motion.mode = DEMCFDCoupler::MotionMode::Dynamic;
		}
		else
		{
			throw ZaranError("DEMCFDSimulation: [coupling] particle_motion must be "
				"fixed | prescribed | dynamic");
		}
		if (m_integrator_name == "exponential")
		{
			m_motion.integrator = DEMCFDCoupler::Integrator::Exponential;
		}
		else if (m_integrator_name == "rk2")
		{
			m_motion.integrator = DEMCFDCoupler::Integrator::Rk2;
		}
		else
		{
			throw ZaranError("DEMCFDSimulation: [coupling] particle_integrator must be "
				"exponential | rk2");
		}
		m_motion.substeps = m_particle_substeps;
		m_motion.max_step_parameter = m_particle_max_step_param;
		for (int d = 0; d < 3; ++d)
		{
			m_motion.body_acceleration[d] = m_body_acceleration[d];
		}
		if (m_motion.mode == DEMCFDCoupler::MotionMode::Dynamic && m_freeze_particles)
		{
			Log::warn("DEMCFDSimulation: particle_motion = dynamic 但 freeze_particles = true，"
				"输入速度仍会被强制为 0（动态积分会把速度算出来，这没有问题，"
				"但如果你想让颗粒带初速，请把 freeze_particles 设为 false）");
		}
		// 反作用"是否包含压力梯度力"：与气相动量方程的写法**必须配套**
		bool pg_on = (m_react_pg_mode == "on");
		if (m_react_pg_mode == "auto")
		{
			pg_on = (m_para->GetPorosityGradientForce() == 0);
		}
		else if (m_react_pg_mode != "off" && m_react_pg_mode != "on")
		{
			throw ZaranError("DEMCFDSimulation: [coupling] react_pressure_gradient_force "
				"must be auto | on | off");
		}
		m_react_pg_resolved = pg_on ? 1 : 0;
		if (m_react_on_fluid && m_react_pg_mode != "auto"
			&& pg_on != (m_para->GetPorosityGradientForce() == 0))
		{
			Log::warn("DEMCFDSimulation: react_pressure_gradient_force 与 "
				"two_phase.porosity_gradient_force 不配套（前者 {}，后者 {}）——"
				" 要么双计、要么漏掉浮力，请确认这是有意为之",
				pg_on ? "on" : "off", m_para->GetPorosityGradientForce());
		}
		m_two_way_active = m_react_on_fluid
			|| m_motion.mode == DEMCFDCoupler::MotionMode::Dynamic;

		if (!DragForceModel::ParseModel(m_drag_model_name, m_options.model))
		{
			throw ZaranError("DEMCFDSimulation: [coupling] drag_model '" + m_drag_model_name
				+ "' is unknown (stokes | schiller_naumann | wen_yu | ergun | gidaspow)");
		}
		if (m_mapping_method == "subcell")
		{
			m_options.mapping = ParticleGridMapper::Method::SubCell;
		}
		else if (m_mapping_method == "exact")
		{
			m_options.mapping = ParticleGridMapper::Method::Exact;
		}
		else
		{
			throw ZaranError("DEMCFDSimulation: [coupling] mapping_method must be subcell | exact");
		}
		if (m_mapping_sub < 1)
		{
			throw ZaranError("DEMCFDSimulation: [coupling] mapping_sub must be >= 1");
		}
		if (m_force_write_interval < 0)
		{
			m_force_write_interval = 0;
		}
		if (!(m_particle_density_fallback > 0.0))
		{
			throw ZaranError("DEMCFDSimulation: [coupling] particle_density must be positive");
		}

		Log::info("DEMCFDSimulation: particle_file={}, drag_model={}, mapping={}(n_sub={}), "
			"smoothing={}, ε from particles={}, freeze_particles={}, "
			"pressure_gradient_force={}, force_write_interval={}",
			m_particle_file, DragForceModel::ModelName(m_options.model), m_mapping_method,
			m_mapping_sub, m_smooth_passes, m_volume_fraction_from_particles, m_freeze_particles,
			m_options.drag.include_pressure_gradient_force ? "on" : "off",
			m_force_write_interval);
		Log::info("DEMCFDSimulation[双向耦合]: particle_motion={}, integrator={}, substeps={}, "
			"react_on_fluid={}, react_pressure_gradient_force={}（与 porosity_gradient_force={} "
			"配套 ⇒ 实际{}）, body_accel=({:g},{:g},{:g}) m/s^2, remap_ε={}",
			m_particle_motion_name, m_integrator_name, m_particle_substeps,
			m_react_on_fluid ? "true" : "false",
			m_react_pg_mode, m_para->GetPorosityGradientForce(),
			m_react_pg_resolved ? "反作用曳力 + 浮力" : "只反作用曳力",
			m_body_acceleration[0], m_body_acceleration[1], m_body_acceleration[2],
			m_remap_each_step);
	}

	void DEMCFDSimulation::BindField()
	{
		if (GetFieldManager()->GetFieldNum() < 1)
		{
			throw ZaranError("DEMCFDSimulation: empty field manager");
		}
		auto field = std::dynamic_pointer_cast<NSFieldStructTwoPhase>(
			GetFieldManager()->GetField(0));
		if (!field)
		{
			throw ZaranError("DEMCFDSimulation: 需要两相场 NSFieldStructTwoPhase"
				"（task.solver = \"EulerDEM\"）");
		}
		m_field = field;
		m_solver = field->GetTwoPhaseSolver();
		m_para = field->GetTwoPhaseSolverPara();
		m_data_manager = field->GetTwoPhaseDataManager();
		if (!m_solver || !m_para || !m_data_manager)
		{
			throw ZaranError("DEMCFDSimulation: inconsistent two-phase field");
		}

		// 外部注入 ε 是耦合的前提：ε 由耦合器写、求解器只补 ghost。
		// 若控制文件写成解析场，耦合器写的值会被 InitField() 覆盖掉，
		// 表现为"粒子完全不起作用"——这里直接拦下来。
		if (m_para->GetVolumeFractionType() != VolumeFractionType::External)
		{
			throw ZaranError("DEMCFDSimulation: 需要 init.volume_fraction.type = \"external\""
				"（ε_g 由 DEM 侧注入），当前是解析场");
		}
	}

	void DEMCFDSimulation::ReadParticles()
	{
		const std::string work_dir = GlobalData::GetString("work_dir");
		const std::string path = work_dir + "/" + m_particle_file;
		if (!IsFileExist(path))
		{
			throw ZaranError("DEMCFDSimulation: particle file not found: " + path);
		}
		ReadDEMParticle reader;
		std::vector<DEMParticle> raw;
		reader.ReadCSV(path, raw);
		if (raw.empty())
		{
			throw ZaranError("DEMCFDSimulation: particle file has no particle: " + path);
		}

		std::vector<DEMCFDCoupler::Particle> particles;
		particles.reserve(raw.size());
		int zero_velocity_disabled = 0;
		for (const auto& p : raw)
		{
			if (!(p.radius > 0.0))
			{
				continue;
			}
			DEMCFDCoupler::Particle q;
			q.id = p.id;
			q.group = p.group;
			q.x = p.pos.x();
			q.y = p.pos.y();
			q.z = p.pos.z();
			q.radius = p.radius;

			const double vp = ParticleGridMapper::ParticleVolume(3, p.radius);
			// 密度优先由 mass/V_p 给出（与 DEM 输入口径一致）；mass 为 0 时退回参数值
			q.density = (p.mass > 0.0 && vp > 0.0)
				? (p.mass / vp) : m_particle_density_fallback;

			for (int d = 0; d < 3; ++d)
			{
				q.velocity[d] = m_freeze_particles ? 0.0 : p.vel[d];
			}
			if (m_freeze_particles
				&& (std::fabs(p.vel.x()) + std::fabs(p.vel.y()) + std::fabs(p.vel.z()) > 0.0))
			{
				++zero_velocity_disabled;
			}
			particles.push_back(q);
		}
		if (particles.empty())
		{
			throw ZaranError("DEMCFDSimulation: no particle with a positive radius");
		}
		if (zero_velocity_disabled > 0)
		{
			Log::warn("DEMCFDSimulation: freeze_particles = true，{} 个粒子的输入速度被强制为 0"
				"（本阶段粒子保持静止，不推进位移）", zero_velocity_disabled);
		}

		m_outside_count = m_coupler.SetParticles(particles);
	}

	// ==================================================================
	// 耦合并入点
	// ==================================================================
	void DEMCFDSimulation::PrepareCoupling()
	{
		// 先绑场（ReadParameters 需要读 porosity_gradient_force 来判断反作用口径）
		BindField();
		ReadParameters();

		// ---- 1) 建立 CFD 网格视图（注意：在 solver->Init() 之前，m_dx 还没算出来，
		//         所以间距直接从节点坐标取，保证这一阶段也能用） ----
		auto grid = m_field->GetGrid();
		index_type is, ie, js, je, ks, ke;
		grid->GetRange(is, ie, js, je, ks, ke);
		auto node = grid->GetNode();

		DEMCFDCoupler::GridView view;
		view.dim = 3;
		view.ni = static_cast<int>(grid->GetNi());
		view.nj = static_cast<int>(grid->GetNj());
		view.nk = static_cast<int>(grid->GetNk());
		view.is = static_cast<int>(is);
		view.ie = static_cast<int>(ie);
		view.js = static_cast<int>(js);
		view.je = static_cast<int>(je);
		view.ks = static_cast<int>(ks);
		view.ke = static_cast<int>(ke);
		if (view.PhysicalNi() < 2 || view.PhysicalNj() < 2 || view.PhysicalNk() < 2)
		{
			throw ZaranError("DEMCFDSimulation: 每个方向的物理节点数至少要有 2 个");
		}
		const auto c000 = node->GetCoord(is, js, ks);
		const auto c100 = node->GetCoord(is + 1, js, ks);
		const auto c010 = node->GetCoord(is, js + 1, ks);
		const auto c001 = node->GetCoord(is, js, ks + 1);
		view.x0 = c000[0];
		view.y0 = c000[1];
		view.z0 = c000[2];
		view.dx = c100[0] - c000[0];
		view.dy = c010[1] - c000[1];
		view.dz = c001[2] - c000[2];

		// ---- 2) 单位换算：用求解器自己的参考量，避免"另一处手写的 a_ref" ----
		const Dimensionless& ref = m_para->GetDimensionless();
		DEMCFDCoupler::UnitScales scales;
		scales.ref_length = ref.GetRefLength();
		scales.ref_density = ref.GetRefDensity();
		scales.ref_velocity = ref.GetRefVelocity();
		scales.ref_temperature = ref.GetRefTemp();
		scales.ref_gamma = ref.GetRefGamma();
		Log::info("DEMCFDSimulation: 参考量 L_ref={:E} m, ρ_ref={:E} kg/m^3, "
			"a_ref={:E} m/s, T_ref={:E} K ⇒ p_ref={:E} Pa, F_ref={:E} N",
			scales.ref_length, scales.ref_density, scales.ref_velocity,
			scales.ref_temperature, scales.RefPressure(), scales.RefForce());

		// ---- 3) 建耦合器、灌粒子 ----
		m_coupler.Init(view, scales, m_options);
		ReadParticles();
		m_coupler.GetResolutionRatios(m_dx_over_dp);

		// ---- 4) 粒子 → ε_g，并注入求解器的体积分数数组 ----
		m_coupler.ComputeGasVolumeFraction(m_eps_nodes, m_eps_min, m_eps_max);
		InjectVolumeFraction();

		Log::info("DEMCFDSimulation: 耦合已就绪（流体 → 粒子：曳力 + 压力梯度力；"
			"粒子 → 流体：{}）",
			m_react_on_fluid
			? (m_react_work ? "反作用力（动量）+ 相间功（能量）" : "反作用力（动量，未加相间功）")
			: "仅体积排斥，无动量反馈");
	}

	void DEMCFDSimulation::InjectVolumeFraction()
	{
		double* volume_fraction = m_data_manager->GetVolumeFraction();
		if (volume_fraction == nullptr)
		{
			throw ZaranError("DEMCFDSimulation: DataManagerNSTwoPhase::GetVolumeFraction() "
				"返回空指针");
		}
		for (size_t i = 0; i < m_eps_nodes.size(); ++i)
		{
			volume_fraction[i] = m_eps_nodes[i];
		}
		// Σ α_s V_cell 应当等于 Σ V_p（映射器的构造保证），这里把残差记进报告
		const double part_vol = m_coupler.GetParticleVolume();
		if (part_vol > 0.0 && !m_eps_nodes.empty())
		{
			const auto& mapper = m_coupler.GetMapper();
			const double assigned = mapper.GetStats().assigned_volume;
			m_partition_residual = std::fabs(assigned - part_vol) / part_vol;
			m_alpha_max = mapper.GetStats().alpha_max;
		}
		Log::info("DEMCFDSimulation: ε_g 已注入（物理节点 {} 个，含 ghost 共 {} 个），"
			"ε ∈ [{:.6f}, {:.6f}]", m_coupler.GetGrid().PhysicalNi()
			* m_coupler.GetGrid().PhysicalNj() * m_coupler.GetGrid().PhysicalNk(),
			m_eps_nodes.size(), m_eps_min, m_eps_max);
	}

	void DEMCFDSimulation::RemapVolumeFraction()
	{
		// 颗粒运动后 ε 也要跟着变。放在 CouplingPreStep（即 PreSolve → BoundaryCondition
		// 之前）调用：那里求解器还没开始这一步，ghost 会在 BoundaryCondition() 里
		// 按新的物理节点值重新零梯度外推，因此不需要在这里碰 ghost。
		m_coupler.ComputeGasVolumeFraction(m_eps_nodes, m_eps_min, m_eps_max, false);
		InjectVolumeFraction();
	}

	void DEMCFDSimulation::CouplingPreStep(int iter)
	{
		if (!m_two_way_active)
		{
			// 纯单向耦合（阶段 3）：这一步什么都不做，行为与阶段 3 完全一致
			(void)iter;
			return;
		}
		// ---- 1) 颗粒运动过 ⇒ ε 变了，先重建并注入 ----
		if (m_remap_each_step && m_motion.mode != DEMCFDCoupler::MotionMode::Fixed)
		{
			RemapVolumeFraction();
		}
		// ---- 2) 用**步初**气态评估受力（这一份就是下面要散射进流体的反作用力）----
		EvaluateForces();
		m_results_pre = m_results;
		// ---- 3) 反作用力散列 → 求解器动量源项 ----
		if (m_react_on_fluid)
		{
			m_coupler.ScatterReaction(m_results_pre, m_react_pg_resolved != 0,
				m_source_star, m_scatter_summary);
			m_solver->SetInterphaseMomentumSource(m_source_star.data());
			m_solver->SetInterphaseWork(m_react_work);
			m_weight_deficit = std::max(m_weight_deficit,
				m_scatter_summary.weight_deficit);
		}
		else
		{
			m_solver->SetInterphaseMomentumSource(nullptr);
		}
	}

	void DEMCFDSimulation::CouplingPostStep(int iter)
	{
		// 注意：基类的迭代计数在一步里被递增了两次（`SolveField()` 的循环里一次、
		// `PreSolve()` 里再一次），所以**钩子拿到的 iter 与流场输出文件名的编号差 1**。
		// 流场帧名、backup/iter=N 目录用的都是 GlobalData 里的值，这里统一跟随它，
		// 否则 particle_force.csv 的帧号与 result/<iter>.dat 对不上。
		const int iter_out = GlobalData::GetInt("iteration.current_iter");
		(void)iter;
		const double time = GlobalData::GetDouble("iteration.current_time");
		// ⚠ 不能用 iteration.dt：`NSFieldSimulation::CalcTimeStep()` 先把**未削减**的
		//    dt 写进 GlobalData，之后才为"正好落在 end_time"把它削减并单独更新
		//    current_time。于是最后一步的 iteration.dt 与实际推进的时间不一致
		//    （实测末步 3.52e-3 vs 实际 2.80e-3，颗粒被多推 26% 的一步）。
		//    改成取 current_time 的增量 —— 那才是气体真正走过的物理时间。
		const double dt = time - m_last_time_for_dt;
		m_last_time_for_dt = time;

		// ---- 1) 用**步初**受力推进颗粒，并核对反作用的动量账本 ----
		// 推进用的是 m_results_pre（与散射进流体的 F^n 严格同一份），
		// 这是标准的交错（lagged）耦合：气相与颗粒在 [t^n, t^{n+1}] 上都用 t^n 的力。
		if (m_two_way_active && dt > 0.0 && !m_results_pre.empty())
		{
			if (m_react_on_fluid)
			{
				double src_imp[3];
				m_coupler.SumSourceMomentum(m_source_star, src_imp);
				for (int d = 0; d < 3; ++d)
				{
					m_source_impulse[d] += src_imp[d] * dt;
					m_particle_impulse[d] -=
						m_scatter_summary.particle_force_dimensionless[d] * dt;
				}
				// 残差的分母取 Σ_p Σ_d |F*_d|（L1 尺度），**不能**取 ΣF* ——
				// 粒子布置对称时 ΣF* 本身恒为 0（这本身是对称性的好检验），
				// 拿它当分母会得到 0/0 的虚假大残差。
				double denom = m_scatter_summary.particle_force_l1;
				double numer = 0.0;
				for (int d = 0; d < 3; ++d)
				{
					numer = std::max(numer, std::fabs(src_imp[d]
						+ m_scatter_summary.particle_force_dimensionless[d]));
				}
				if (denom > 0.0)
				{
					m_max_scatter_residual = std::max(m_max_scatter_residual, numer / denom);
				}
			}
			if (m_motion.mode == DEMCFDCoupler::MotionMode::Dynamic)
			{
				m_last_substeps = m_coupler.AdvanceParticles(dt, m_motion, m_results_pre);
				++m_steps_advanced;
			}
		}

		// ---- 2) 在**当前**（步末）气态上重算受力，供输出与下游校验 ----
		// 单向耦合与双向耦合都走这一条：CSV 的每一行对应的都是该帧时刻的受力，
		// 因此阶段 3 的 V4/V5/V6/V8（把 CSV 与 result/<iter>.dat 逐值对照）依然成立。
		EvaluateForces();

		// ---- 3) 输出 ----
		int write_interval = m_force_write_interval;
		if (write_interval == 0)
		{
			write_interval = GlobalData::GetInt("iteration.write_interval");
		}
		if (write_interval > 0 && iter_out % write_interval == 0)
		{
			WriteForceHistory(iter_out, time);
		}
		const int flow_interval = GlobalData::GetInt("iteration.write_interval");
		if (m_write_coupling_frames && flow_interval > 0 && iter_out % flow_interval == 0)
		{
			WriteForceFrame(iter_out, time);
		}

		m_last_iter = iter_out;
		m_last_time = time;
	}

	// ==================================================================
	// 受力评估
	// ==================================================================
	void DEMCFDSimulation::UpdateViscosityField()
	{
		const double* density = m_data_manager->GetPrim(ID_DENSITY);
		const double* pressure = m_data_manager->GetPrim(ID_PRESSURE);
		// μ 走求解器同一套 Sutherland（Gas 由参数类持有，是求解器 InitSolver() 里
		// 建的那个对象，所以这里的 μ 与求解器内部完全一致）
		Gas* gas = m_para->GetGas();
		if (gas == nullptr)
		{
			throw ZaranError("DEMCFDSimulation: 求解器的 Gas 尚未建立");
		}
		const double gamma = gas->GetGamma();
		const size_t n = m_eps_nodes.size();
		m_viscosity.assign(n, 0.0);
		for (size_t i = 0; i < n; ++i)
		{
			const double r = density[i];
			if (!(r > 0.0) || !(pressure[i] > 0.0))
			{
				continue;
			}
			// T* = γ p*/ρ*（与 Gas::CalcTemperature 同一口径）；夹一下避免 Sutherland 溢出
			double t = gamma * pressure[i] / r;
			t = std::min(1.0e6, std::max(1.0e-6, t));
			m_viscosity[i] = gas->CalcMu(t);
		}
	}

	void DEMCFDSimulation::EvaluateForces()
	{
		UpdateViscosityField();

		DEMCFDCoupler::PrimitiveView prim;
		prim.density = m_data_manager->GetPrim(ID_DENSITY);
		prim.velocity[0] = m_data_manager->GetPrim(ID_VELOCITY_X);
		prim.velocity[1] = m_data_manager->GetPrim(ID_VELOCITY_Y);
		prim.velocity[2] = m_data_manager->GetPrim(ID_VELOCITY_Z);
		prim.pressure = m_data_manager->GetPrim(ID_PRESSURE);
		prim.viscosity = m_viscosity.data();

		// ∇p 在耦合器内部算（EvaluateParticleForces 会自己调用 ComputePressureGradient），
		// 这里只需要把节点量交给它 —— 插值与体积分配用同一套权重。
		m_coupler.EvaluateParticleForces(prim, m_results);

		m_last_sum_force[0] = m_last_sum_force[1] = m_last_sum_force[2] = 0.0;
		m_sum_abs_force = 0.0;
		for (size_t i = 0; i < m_results.size(); ++i)
		{
			const auto& r = m_results[i];
			const double f = std::sqrt(r.force.force[0] * r.force.force[0]
				+ r.force.force[1] * r.force.force[1]
				+ r.force.force[2] * r.force.force[2]);
			m_max_force = std::max(m_max_force, f);
			m_sum_abs_force += f;
			for (int d = 0; d < 3; ++d)
			{
				m_last_sum_force[d] += r.force.force[d];
			}
			if (i < m_coupler.GetParticles().size())
			{
				const auto& p = m_coupler.GetParticles()[i];
				const double s = std::sqrt(p.velocity[0] * p.velocity[0]
					+ p.velocity[1] * p.velocity[1] + p.velocity[2] * p.velocity[2]);
				m_max_particle_speed = std::max(m_max_particle_speed, s);
			}
		}
	}

	// ==================================================================
	// 输出
	// ==================================================================
	void DEMCFDSimulation::WriteForceHistory(int iter, double time)
	{
		const std::string dir = GlobalData::GetString("work_dir") + "/"
			+ GlobalData::GetString("output.result_folder");
		const std::string path = dir + "/particle_force.csv";

		std::ofstream out;
		if (!m_header_written)
		{
			out.open(path, std::ios::trunc);
			if (!out.is_open())
			{
				Log::warn("DEMCFDSimulation: cannot write '{}'", path);
				return;
			}
			out << kForceHeader;
			m_header_written = true;
		}
		else
		{
			out.open(path, std::ios::app);
			if (!out.is_open())
			{
				Log::warn("DEMCFDSimulation: cannot append '{}'", path);
				return;
			}
		}

		const auto& particles = m_coupler.GetParticles();
		// 全部浮点列用 %.17g：校验脚本要拿这些值做独立复算，9 位有效数字证明不了任何事
		//（与"逐位退化检验必须把输出精度提到 17 位"是同一条理由）
		char line[4096];
		for (size_t i = 0; i < particles.size() && i < m_results.size(); ++i)
		{
			const auto& p = particles[i];
			const auto& r = m_results[i];
			const auto& s = r.gas;
			const auto& f = r.force;
			const double fmag = std::sqrt(f.force[0] * f.force[0] + f.force[1] * f.force[1]
				+ f.force[2] * f.force[2]);
			const double mass = m_coupler.ParticleMass(i);
			// 颗粒加速度 = (曳力 + 压力梯度力)/m + 体积力加速度
			double acc[3];
			for (int d = 0; d < 3; ++d)
			{
				acc[d] = (mass > 0.0)
					? f.force[d] / mass + m_body_acceleration[d] : 0.0;
			}
			// 实际散射进流体的反作用力（= −(曳力 [+ 浮力])，与散列口径一致）
			double freact[3];
			for (int d = 0; d < 3; ++d)
			{
				freact[d] = -(f.drag[d] + (m_react_pg_resolved ? f.pressure_gradient[d] : 0.0));
			}
			std::snprintf(line, sizeof(line),
				// ⚠ 说明符个数必须与 kForceHeader 的列数**逐一对齐**（47 个）。
				//   这里曾经多写了一个 %.17g：前面 24 列类型都还对得上，于是
				//   `branch` 拿到一个 int 却被当成 %.17g 打印（写出 4.9e-324 这种次正规数），
				//   `beta_over_es` 拿到 double 却被 %d 打印成随机整数。
				//   错位后的行**仍然解析得动**（都是数字），只能靠
				//   "表头列数 == 数据列数" 这条断言抓 —— 见 verify.py 的 V1。
				"%d,%.17g,%llu,%d,"
				"%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%lld,"
				"%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
				"%.17g,%.17g,%.17g,"
				"%.17g,%.17g,%d,%.17g,"
				"%.17g,%.17g,%.17g,"
				"%.17g,%.17g,%.17g,"
				"%.17g,%.17g,%.17g,%.17g,"
				"%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%.17g,%.17g,%.17g\n",
				iter, time,
				static_cast<unsigned long long>(p.id), p.group,
				p.x, p.y, p.z, p.radius, p.density,
				s.eps_gas, s.weight_sum, s.cells_touched,
				s.density, s.velocity[0], s.velocity[1], s.velocity[2],
				s.pressure, s.temperature, s.viscosity,
				s.grad_p[0], s.grad_p[1], s.grad_p[2],
				f.re, f.cd, f.branch, f.beta_over_es,
				f.drag[0], f.drag[1], f.drag[2],
				f.pressure_gradient[0], f.pressure_gradient[1], f.pressure_gradient[2],
				f.force[0], f.force[1], f.force[2], fmag,
				mass,
				p.velocity[0], p.velocity[1], p.velocity[2],
				acc[0], acc[1], acc[2],
				m_last_substeps,
				freact[0], freact[1], freact[2]);
			out << line;
		}
	}

	void DEMCFDSimulation::WriteForceFrame(int iter, double time)
	{
		Visual vis;
		const std::string dir = GlobalData::GetString("work_dir") + "/"
			+ GlobalData::GetString("output.result_folder");
		const auto& particles = m_coupler.GetParticles();
		const int n = static_cast<int>(particles.size());

		std::vector<double> x(n), y(n), z(n), radius(n), density(n), eps(n);
		std::vector<double> re(n), cd(n), rho(n), pres(n), ugx(n), ugy(n), ugz(n);
		std::vector<double> fx(n), fy(n), fz(n), fmag(n);
		std::vector<double> upx(n), upy(n), upz(n), urx(n), ury(n), urz(n);
		std::vector<double> rfx(n), rfy(n), rfz(n), mass(n);
		for (int i = 0; i < n; ++i)
		{
			const auto& p = particles[i];
			const auto& r = m_results[i];
			x[i] = p.x; y[i] = p.y; z[i] = p.z;
			radius[i] = p.radius; density[i] = p.density;
			eps[i] = r.gas.eps_gas;
			re[i] = r.force.re; cd[i] = r.force.cd;
			rho[i] = r.gas.density; pres[i] = r.gas.pressure;
			ugx[i] = r.gas.velocity[0]; ugy[i] = r.gas.velocity[1]; ugz[i] = r.gas.velocity[2];
			fx[i] = r.force.force[0]; fy[i] = r.force.force[1]; fz[i] = r.force.force[2];
			fmag[i] = std::sqrt(fx[i] * fx[i] + fy[i] * fy[i] + fz[i] * fz[i]);
			upx[i] = p.velocity[0]; upy[i] = p.velocity[1]; upz[i] = p.velocity[2];
			urx[i] = ugx[i] - upx[i]; ury[i] = ugy[i] - upy[i]; urz[i] = ugz[i] - upz[i];
			rfx[i] = -(r.force.drag[0] + (m_react_pg_resolved ? r.force.pressure_gradient[0] : 0.0));
			rfy[i] = -(r.force.drag[1] + (m_react_pg_resolved ? r.force.pressure_gradient[1] : 0.0));
			rfz[i] = -(r.force.drag[2] + (m_react_pg_resolved ? r.force.pressure_gradient[2] : 0.0));
			mass[i] = m_coupler.ParticleMass(static_cast<size_t>(i));
		}

		std::vector<Visual::PointScalarField> fields = {
			{ "X", &x }, { "Y", &y }, { "Z", &z },
			{ "Radius", &radius }, { "Density", &density }, { "Mass", &mass },
			{ "Eps_gas", &eps },
			{ "Rho_gas", &rho }, { "P_gas", &pres },
			{ "U_gas_x", &ugx }, { "U_gas_y", &ugy }, { "U_gas_z", &ugz },
			{ "U_p_x", &upx }, { "U_p_y", &upy }, { "U_p_z", &upz },
			{ "U_rel_x", &urx }, { "U_rel_y", &ury }, { "U_rel_z", &urz },
			{ "Re", &re }, { "Cd", &cd },
			{ "Fx", &fx }, { "Fy", &fy }, { "Fz", &fz }, { "F_mag", &fmag },
			{ "Freact_x", &rfx }, { "Freact_y", &rfy }, { "Freact_z", &rfz },
		};
		char zone[160];
		std::snprintf(zone, sizeof(zone),
			"particle force iter=%d t=%.4f (%s/%s)", iter, time,
			DragForceModel::ModelName(m_options.model), m_particle_motion_name.c_str());
		vis.WritePointsTecplotASCII(dir + "/particles_force_" + std::to_string(iter) + ".dat",
			"Zaran3 DEM-CFD one-way coupling: particle forces",
			zone, n, fields, time);
	}

	void DEMCFDSimulation::WriteReport() const
	{
		const std::string dir = GlobalData::GetString("work_dir") + "/"
			+ GlobalData::GetString("output.result_folder");
		std::string buf;
		buf += "key,value\n";
		char line[256];
		auto add = [&](const char* k, const char* fmt, double v)
		{
			std::snprintf(line, sizeof(line), "%s,", k);
			buf += line;
			std::snprintf(line, sizeof(line), fmt, v);
			buf += line;
			buf += "\n";
		};
		auto add_int = [&](const char* k, long long v)
		{
			buf += k;
			buf += ",";
			buf += std::to_string(v);
			buf += "\n";
		};
		auto add_str = [&](const char* k, const std::string& v)
		{
			buf += k;
			buf += ",";
			buf += v;
			buf += "\n";
		};

		add_str("drag_model", DragForceModel::ModelName(m_options.model));
		add_str("mapping_method", m_mapping_method);
		add_int("mapping_sub", m_mapping_sub);
		// ---- 双向耦合 ----
		add_str("particle_motion", m_particle_motion_name);
		add_str("particle_integrator", m_integrator_name);
		add_int("react_on_fluid", m_react_on_fluid ? 1 : 0);
		add_int("react_pressure_gradient_force", m_react_pg_resolved);
		add("body_accel_x_m_s2", "%.12e", m_body_acceleration[0]);
		add("body_accel_y_m_s2", "%.12e", m_body_acceleration[1]);
		add("body_accel_z_m_s2", "%.12e", m_body_acceleration[2]);
		add_int("steps_advanced_particles", m_steps_advanced);
		add("max_particle_speed_m_s", "%.12e", m_max_particle_speed);
		// 动量账本：反作用力散列（构造保证 Σ_cells src·V_cell = −Σ_p F*）
		add("source_impulse_x", "%.12e", m_source_impulse[0]);
		add("source_impulse_y", "%.12e", m_source_impulse[1]);
		add("source_impulse_z", "%.12e", m_source_impulse[2]);
		add("particle_impulse_x", "%.12e", m_particle_impulse[0]);
		add("particle_impulse_y", "%.12e", m_particle_impulse[1]);
		add("particle_impulse_z", "%.12e", m_particle_impulse[2]);
		add("max_scatter_residual", "%.6e", m_max_scatter_residual);
		add("weight_deficit", "%.6e", m_weight_deficit);
		add("force_ref_N", "%.12e", m_coupler.GetScales().RefForce());
		// 无量纲参考量：校验脚本只有拿到它们才能把 SI 的受力/气态换算回求解器的无量纲量
		const auto& sc = m_coupler.GetScales();
		add("ref_length_m", "%.12e", sc.ref_length);
		add("ref_density_kg_m3", "%.12e", sc.ref_density);
		add("ref_velocity_m_s", "%.12e", sc.ref_velocity);
		add("ref_temperature_K", "%.12e", sc.ref_temperature);
		add("ref_gamma", "%.12e", sc.ref_gamma);
		add("ref_pressure_Pa", "%.12e", sc.RefPressure());
		add("ref_force_N", "%.12e", sc.RefForce());
		add("grid_dx_m", "%.12e", m_coupler.GetGrid().dx);
		add("grid_dy_m", "%.12e", m_coupler.GetGrid().dy);
		add("grid_dz_m", "%.12e", m_coupler.GetGrid().dz);
		add_int("particle_count", static_cast<long long>(m_coupler.GetParticles().size()));
		add_int("particle_outside_domain", m_outside_count);
		add("particle_volume_m3", "%.12e", m_coupler.GetParticleVolume());
		add("rep_diameter_m", "%.12e", m_coupler.GetRepDiameter());
		add("dx_over_dp", "%.9f", m_dx_over_dp[0]);
		add("dy_over_dp", "%.9f", m_dx_over_dp[1]);
		add("dz_over_dp", "%.9f", m_dx_over_dp[2]);
		add("eps_gas_min", "%.12e", m_eps_min);
		add("eps_gas_max", "%.12e", m_eps_max);
		add("alpha_solid_max", "%.12e", m_alpha_max);
		add("partition_residual", "%.6e", m_partition_residual);
		add("include_pressure_gradient_force", "%.0f",
			m_options.drag.include_pressure_gradient_force ? 1.0 : 0.0);
		add("volume_fraction_from_particles", "%.0f",
			m_options.volume_fraction_from_particles ? 1.0 : 0.0);
		add("last_iter", "%.0f", static_cast<double>(m_last_iter));
		add("last_time", "%.12e", m_last_time);
		add("max_particle_force_N", "%.12e", m_max_force);
		add("sum_abs_particle_force_N", "%.12e", m_sum_abs_force);
		add("last_sum_Fx_N", "%.12e", m_last_sum_force[0]);
		add("last_sum_Fy_N", "%.12e", m_last_sum_force[1]);
		add("last_sum_Fz_N", "%.12e", m_last_sum_force[2]);

		WriteText(dir + "/coupling_report.csv", buf);
		Log::info("DEMCFDSimulation: report written to '{}/coupling_report.csv'", dir);
		Log::info("DEMCFDSimulation: 最大相间力 |F|max = {:E} N，末态 ΣF = ({:E}, {:E}, {:E}) N",
			m_max_force, m_last_sum_force[0], m_last_sum_force[1], m_last_sum_force[2]);
	}
}
