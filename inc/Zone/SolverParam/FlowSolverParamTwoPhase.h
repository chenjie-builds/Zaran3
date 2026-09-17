/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file FlowSolverParamTwoPhase.h
 * \brief Solver parameters for the two-phase (volume-fraction) Euler solver.
 * \author Chen Jie.
 *
 *  本参数类描述**气相体积分数场** ε_g(x)（阶段 1 中为冻结的规定场，不含颗粒），
 *  以及体积分数模型的建模开关。继承 FlowSolverParamUniform 以复用 CFL / RK /
 *  限制器 / Riemann 求解器 / 初值类型等全部既有参数解析。
 */

#pragma once
#include "FlowSolverParamUniform.h"

namespace zaran
{
	/// @brief 体积分数的初值分布类型
	enum class VolumeFractionType
	{
		Uniform,   ///< 全场常值
		Step,      ///< 沿 x 的阶跃（两段常值，膜片位置可指定）
		Sine,      ///< 沿 x 的正弦分布，用于检验 ε 梯度项与静水平衡
	};

	/// @brief 两相（体积分数）均匀结构网格 Euler 求解器的参数。
	class FlowSolverParamTwoPhase : public FlowSolverParamUniform
	{
	public:
		FlowSolverParamTwoPhase();
		virtual ~FlowSolverParamTwoPhase();
		/// @brief 从 GlobalData 读取参数
		void Init() override;

	public:
		// ---------------- 体积分数初值 ----------------
		const VolumeFractionType& GetVolumeFractionType() const { return m_volume_fraction_type; }
		const double& GetVolumeFractionValue() const { return m_volume_fraction_value; }
		const double& GetVolumeFractionXStep() const { return m_volume_fraction_x_step; }
		const double& GetVolumeFractionLeft() const { return m_volume_fraction_left; }
		const double& GetVolumeFractionRight() const { return m_volume_fraction_right; }
		const double& GetVolumeFractionMean() const { return m_volume_fraction_mean; }
		const double& GetVolumeFractionAmplitude() const { return m_volume_fraction_amplitude; }
		const double& GetVolumeFractionWavelength() const { return m_volume_fraction_wavelength; }
		const double& GetVolumeFractionPhase() const { return m_volume_fraction_phase; }

		// ---------------- 建模开关 ----------------
		/// @brief 是否加入孔隙率梯度力 +p∇ε（1 = 用 ε∇p 形式，0 = 用守恒通量 ∇(εp) 形式）
		/// @details 两种形式在 ε ≡ 1 与 ε 均匀时完全等价；ε 变化时
		///          - 1：气相动量方程离散为 ε∇p，与 DEM-CFD 惯用的"压力梯度力 -V∇p"配套；
		///          - 0：完全守恒形式 ∂(εU)/∂t + ∇·(εF) = S，动量严格守恒但缺少 p∇ε 力。
		const int& GetPorosityGradientForce() const { return m_porosity_gradient_force; }
		/// @brief 体积分数下限（低于该值视为非物理，仅在检查时报警）
		const double& GetVolumeFractionMin() const { return m_volume_fraction_min; }

	private:
		VolumeFractionType m_volume_fraction_type = VolumeFractionType::Uniform;
		double m_volume_fraction_value = 1.0;
		double m_volume_fraction_x_step = 0.5;
		double m_volume_fraction_left = 1.0;
		double m_volume_fraction_right = 1.0;
		double m_volume_fraction_mean = 0.8;
		double m_volume_fraction_amplitude = 0.2;
		double m_volume_fraction_wavelength = 1.0;
		double m_volume_fraction_phase = 0.0;

		int m_porosity_gradient_force = 1;
		double m_volume_fraction_min = 1.0e-3;
	};
}
