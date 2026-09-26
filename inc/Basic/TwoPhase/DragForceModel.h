/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file DragForceModel.h
 * \brief Inter-phase force models (fluid -> particle) for the one-way coupled DEM-CFD solver.
 * \author Chen Jie.
 *
 *  阶段 3 的相间力模块（见 docs/TWO_PHASE_STAGE3_DEM_CFD.md）
 *  ------------------------------------------------------
 *  本文件只做**物理**：给气相状态、颗粒状态与局部空隙率，返回作用在颗粒上的力。
 *  它不持有网格、不依赖求解器状态、不做时间推进，因此可以脱离 CFD 被验证到解析精度
 *  （`tests/demcfd_shock_tube/verify.py` 的 F 组就是干这件事）。
 *
 *  统一的写法（本文件所有模型都归到同一个式子）
 *  --------------------------------------------
 *      F_drag = V_p · (β/ε_s) · (u_g − u_p)
 *
 *  `β` 是 CFD-DEM 里惯用的"相间动量交换系数"（单位混合物体积），`β/ε_s` 即**单位颗粒
 *  体积**上的曳力系数。它与教科书里"单颗粒曳力 × 空隙率修正"是同一个式子：
 *
 *      β/ε_s = (3/4)·C_D·ρ_g·ε_g^(−n)·|u_rel| / d_p
 *      V_p·(3/4)/d_p = (π d³/6)·(3/4)/d_p = π d²/8 = ½·A_p
 *      ⇒  F_drag = ½·ρ_g·C_D·A_p·|u_rel|·u_rel · ε_g^(−n)
 *
 *  这样五种模型（含稠密支的 Ergun）只需要给出各自的 `β/ε_s`，力的方向与符号处理
 *  只有一份代码，不会出现"某个模型把相对速度写反"这类只差一个符号的 bug。
 *
 *  单位
 *  ----
 *  本类与单位制无关：只要 ρ_g、μ_g、d_p、|u_rel|、∇p 用**同一套自洽的单位**，
 *  输出的力就是该单位制下的力。耦合器（DEMCFDCoupler）负责把求解器的无量纲气态
 *  换算成 SI 再调用，于是输出直接是牛顿。
 */

#pragma once
#include "BasicType.h"

#include <string>

namespace zaran
{
	/// @brief 相间力（流体 → 颗粒，单向耦合）模型集合。
	class DragForceModel
	{
	public:
		// ------------------------------------------------------------
		// 模型
		// ------------------------------------------------------------
		enum class Model
		{
			/// @brief 解析 Stokes 曳力 F = 3πμ d |u_rel|（Re → 0 极限，用于标定）
			/// @details 由 β/ε_s = 18μ/d² 给出，对 |u_rel| = 0 也是良定义的
			Stokes,
			/// @brief 单颗粒 Schiller–Naumann 相关式，**不含**空隙率修正（无界绕流）
			/// @details C_D = 24/Re·(1+0.15Re^0.687)（Re < 1000），Re ≥ 1000 取 0.44
			SchillerNaumann,
			/// @brief Wen–Yu：Schiller–Naumann + 空隙率修正 ε_g^(−2.65)
			/// @details Re 取 ε_g·ρ_g·d·|u_rel|/μ_g（Gidaspow 的原文口径）
			WenYu,
			/// @brief Ergun 压降反推的曳力（稠密支），β/ε_s = 150με_s/(ε_g d²) + 1.75ρ_g|u_rel|/d
			Ergun,
			/// @brief Gidaspow 分段：ε_g ≥ 阈值取 Wen–Yu，否则取 Ergun（默认）
			/// @details 两支在阈值为**分段间断**（Gidaspow 原文即如此，非实现 bug）：
			///          阈值处的跳变幅度会打进日志，便于识别"算例跑到了哪一支"
			Gidaspow,
		};

		/// @brief 气相状态
		struct GasState
		{
			double density = 0.0;                  ///< 密度
			double velocity[3] = { 0.0, 0.0, 0.0 };///< 速度
			double pressure = 0.0;                 ///< 压强（仅用于记录/检查）
			double temperature = 0.0;              ///< 温度（仅用于记录/检查）
			double viscosity = 0.0;                ///< 动力粘性系数 μ
		};

		/// @brief 颗粒状态
		struct ParticleState
		{
			double radius = 0.0;
			double density = 0.0;                  ///< 颗粒密度（单向耦合下不参与受力，仅用于记录）
			double velocity[3] = { 0.0, 0.0, 0.0 };//< 颗粒速度（静止时为 0）
		};

		/// @brief 建模开关
		struct Options
		{
			/// @brief Gidaspow 的分段阈值（空隙率高于它走稀疏支）
			double dense_threshold = 0.8;
			/// @brief 稀疏支求 Re 时是否乘 ε_g（Gidaspow/Wen–Yu 原文为 true）
			bool dilute_re_uses_voidage = true;
			/// @brief 是否叠加压力梯度力 −V_p∇p
			bool include_pressure_gradient_force = true;
			/// @brief 空隙率修正指数 n（Wen–Yu 为 2.65）
			double voidage_exponent = 2.65;
		};

		/// @brief 一次受力评估的全部中间量（全部写出到 particle_force.csv，便于独立复算）
		struct Result
		{
			double drag[3] = { 0.0, 0.0, 0.0 };            ///< 曳力
			double pressure_gradient[3] = { 0.0, 0.0, 0.0 }; ///< 压力梯度力 −V_p∇p
			double force[3] = { 0.0, 0.0, 0.0 };           ///< 合力
			double relative_velocity[3] = { 0.0, 0.0, 0.0 };///< u_g − u_p
			double relative_speed = 0.0;
			double re = 0.0;             ///< 本次评估用的雷诺数（速度为 0 时记为 0）
			double cd = 0.0;             ///< 曳力系数（Ergun 支不定义，记 0）
			double beta_over_es = 0.0;   ///< 实际使用的 β/ε_s（单位颗粒体积的曳力系数）
			double eps_gas = 1.0;        ///< 本次评估用的空隙率
			int    branch = 0;           ///< 0 = 稀疏/单颗粒支，1 = 稠密（Ergun）支
		};

		// ------------------------------------------------------------
		// 接口
		// ------------------------------------------------------------
		/// @brief 解析模型名（大小写不敏感；别名：sn / schiller_naumann / wen_yu / wenyu / gidaspow / ergun / stokes）
		/// @return 是否识别成功
		static bool ParseModel(const std::string& name, Model& model);
		/// @brief 模型名（写进日志与报告）
		static const char* ModelName(Model model);
		/// @brief 该模型是否属于"稠密支"式（Ergun 型，不定义 C_D）
		static bool IsDenseModel(Model model);

		/// @brief 单颗粒 Schiller–Naumann 曳力系数
		/// @details Re < 1000 用 24/Re·(1+0.15Re^0.687)；Re ≥ 1000 用 0.44。
		///          Re = 0 或 Re < 0 返回 0（与"相对速度为零 ⇒ 曳力为零"自洽）。
		static double SchillerNaumannCd(double re);

		/// @brief 单位颗粒体积的曳力系数 β/ε_s（F_drag = V_p·(β/ε_s)·u_rel）
		/// @param model            模型
		/// @param opts             建模开关
		/// @param eps_gas          气相体积分数 ε_g（= 1 − α_s）
		/// @param gas_density      ρ_g
		/// @param gas_viscosity    μ_g
		/// @param diameter         d_p
		/// @param speed            |u_g − u_p|
		/// @param[out] re          本次使用的 Re（速度为 0 时为 0）
		/// @param[out] cd          本次使用的 C_D（Ergun 支为 0）
		/// @param[out] dense_used  是否走了稠密支
		static double DragCoefficientPerVolume(Model model, const Options& opts,
			double eps_gas, double gas_density, double gas_viscosity,
			double diameter, double speed, double& re, double& cd, bool& dense_used);

		/// @brief 评估单个颗粒受到的相间力
		/// @param grad_p 气相压强梯度（长度 3；不需要时传 nullptr）
		static Result Evaluate(Model model, const Options& opts,
			const GasState& gas, const ParticleState& particle,
			double eps_gas, const double grad_p[3]);
	};
}
