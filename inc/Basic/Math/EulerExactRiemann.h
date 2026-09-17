/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file EulerExactRiemann.h
 * \brief Exact solution of the 1-D Euler Riemann problem (Toro, 3rd ed., Ch. 4).
 *
 *  用于两件事：
 *   1) 给基准测试提供"真值"，从而计算 L1 / L2 / Linf 误差范数；
 *   2) 作为精确解 ghost 边界条件，避免出流边界反射污染内部解。
 *
 *  依赖：仅 C++ 标准库（不依赖 Eigen / 网格 / 框架任何其它模块）。
 */

#pragma once
#include <string>

namespace zaran
{
	//! 一维 Euler 方程 Riemann 问题的精确解。
	/*!
	 *  状态以原始变量 (rho, u, p) 描述，两端比热比取同一常数 gamma。
	 *  求解流程（Toro 4.3 节）：
	 *   - 用 Newton 迭代解压强方程 f(p*) = f_L(p*) + f_R(p*) + (u_R - u_L) = 0；
	 *     初值取 PVRS / TRBS 两段式猜测（Toro 4.53-4.55），保证强激波下也收敛；
	 *   - u* 由 (4.51) 恢复；
	 *   - 采样时按 xi = x/t 判定落在左态 / 左波 / 星区 / 右波 / 右态。
	 */
	class EulerExactRiemann
	{
	public:
		//! 单侧状态（原始变量）
		struct State
		{
			double rho = 1.0;
			double u = 0.0;   //!< 只有 x 方向速度（一维问题）
			double p = 1.0;

			State() = default;
			State(double r, double vel, double pres) : rho(r), u(vel), p(pres) {}
		};

		//! 波型
		enum class WaveType
		{
			Shock,      //!< 激波
			Rarefaction //!< 稀疏波
		};

		EulerExactRiemann(const State& left, const State& right, double gamma = 1.4);

		//! 星区压强 p*
		double GetPStar() const { return m_p_star; }
		//! 星区速度 u*
		double GetUStar() const { return m_u_star; }
		//! 左右两侧的波型
		WaveType GetWaveType(int side) const { return m_shock[side] ? WaveType::Shock : WaveType::Rarefaction; }
		//! Newton 迭代次数与最终残差（用于报告求解器本身的收敛性）
		int GetNewtonIter() const { return m_newton_iter; }
		double GetNewtonResidual() const { return m_newton_res; }

		//! 在相似坐标 xi = x / t 处取样
		State Sample(double xi) const;
		//! 在物理坐标 (x, t) 处取样；t <= 0 时返回左右初始态的阶跃判断
		State Sample(double x, double t) const;

		//! 波型文字描述，用于报告（"shock" / "rarefaction"）
		static const char* WaveTypeName(WaveType type);

		//! 声速
		double SoundSpeed(const State& s) const;

	private:
		//! 压强函数 f_K(p)，Toro (4.39) 激波 / (4.47) 稀疏
		double PressureFunction(double p, const State& s) const;
		//! 压强函数的导数 f'_K(p)，Toro (4.40) / (4.48)
		double PressureFunctionDerivative(double p, const State& s) const;
		//! 星区密度，Toro (4.50) 激波 / (4.55) 稀疏
		double StarDensity(int side) const;
		//! 激波速度 S_K，Toro (4.64)
		double ShockSpeed(int side) const;
		//! 稀疏波头/尾速度，Toro (4.56)
		void RarefactionSpeed(int side, double& head, double& tail, double& c_star) const;
		//! 求解 p* / u*
		void SolveStarState();

		State m_state[2];
		double m_gamma = 1.4;
		bool m_shock[2] = { false, false };
		double m_p_star = 0.0;
		double m_u_star = 0.0;
		int m_newton_iter = 0;
		double m_newton_res = 0.0;
	};
}
