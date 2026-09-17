/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file EulerExactRiemann.cpp
 * \brief Exact solution of the 1-D Euler Riemann problem.
 */

#include "EulerExactRiemann.h"
#include <cmath>
#include <algorithm>

namespace zaran
{
	// 迭代控制
	namespace
	{
		constexpr int    kMaxNewtonIter = 100;
		constexpr double kTolerance = 1.0e-12;
		constexpr double kMinPressure = 1.0e-14;
	}

	EulerExactRiemann::EulerExactRiemann(const State& left, const State& right, double gamma)
		: m_gamma(gamma)
	{
		m_state[0] = left;
		m_state[1] = right;
		if (m_gamma <= 1.0)
		{
			m_gamma = 1.4;
		}
		SolveStarState();
	}

	double EulerExactRiemann::SoundSpeed(const State& s) const
	{
		return std::sqrt(m_gamma * s.p / s.rho);
	}

	const char* EulerExactRiemann::WaveTypeName(WaveType type)
	{
		return (type == WaveType::Shock) ? "shock" : "rarefaction";
	}

	// ------------------------------------------------------------------
	// 压强函数及其导数：Toro (4.39)/(4.40) 激波， (4.47)/(4.48) 稀疏波
	// ------------------------------------------------------------------
	double EulerExactRiemann::PressureFunction(double p, const State& s) const
	{
		const double c = SoundSpeed(s);
		if (p > s.p)
		{
			const double a = 2.0 / ((m_gamma + 1.0) * s.rho);
			const double b = (m_gamma - 1.0) / (m_gamma + 1.0) * s.p;
			return (p - s.p) * std::sqrt(a / (p + b));
		}
		const double exponent = (m_gamma - 1.0) / (2.0 * m_gamma);
		return 2.0 * c / (m_gamma - 1.0) * (std::pow(p / s.p, exponent) - 1.0);
	}

	double EulerExactRiemann::PressureFunctionDerivative(double p, const State& s) const
	{
		const double c = SoundSpeed(s);
		if (p > s.p)
		{
			const double a = 2.0 / ((m_gamma + 1.0) * s.rho);
			const double b = (m_gamma - 1.0) / (m_gamma + 1.0) * s.p;
			return std::sqrt(a / (p + b)) * (1.0 - 0.5 * (p - s.p) / (p + b));
		}
		const double exponent = (m_gamma + 1.0) / (2.0 * m_gamma);
		return 1.0 / (s.rho * c) * std::pow(p / s.p, -exponent);
	}

	double EulerExactRiemann::StarDensity(int side) const
	{
		const State& s = m_state[side];
		if (m_shock[side])
		{
			// Toro (4.50)
			const double ratio = m_p_star / s.p;
			const double nu = (m_gamma - 1.0) / (m_gamma + 1.0);
			return s.rho * (ratio + nu) / (nu * ratio + 1.0);
		}
		// Toro (4.55)
		return s.rho * std::pow(m_p_star / s.p, 1.0 / m_gamma);
	}

	double EulerExactRiemann::ShockSpeed(int side) const
	{
		const State& s = m_state[side];
		const double c = SoundSpeed(s);
		const double factor = std::sqrt(((m_gamma + 1.0) / (2.0 * m_gamma)) * (m_p_star / s.p)
			+ (m_gamma - 1.0) / (2.0 * m_gamma));
		return (side == 0) ? (s.u - c * factor) : (s.u + c * factor);
	}

	void EulerExactRiemann::RarefactionSpeed(int side, double& head, double& tail, double& c_star) const
	{
		const State& s = m_state[side];
		const double c = SoundSpeed(s);
		c_star = c * std::pow(m_p_star / s.p, (m_gamma - 1.0) / (2.0 * m_gamma));
		if (side == 0)
		{
			head = s.u - c;
			tail = m_u_star - c_star;
		}
		else
		{
			head = s.u + c;
			tail = m_u_star + c_star;
		}
	}

	// ------------------------------------------------------------------
	// 求解 p* / u*
	//
	//  压强方程 f(p) = f_L(p) + f_R(p) + (u_R - u_L) 在 p > 0 上严格单调递增
	//  （Toro 定理 4.1），因此这里用"括号保护的安全 Newton"：
	//  每步先收缩括号 [lo, hi]，Newton 步若越出括号（或非有限）就退化为二分。
	//  这样在强双稀疏波（p* -> 0，f' -> ∞，纯 Newton 几乎不动）下也保证收敛。
	// ------------------------------------------------------------------
	void EulerExactRiemann::SolveStarState()
	{
		const State& left = m_state[0];
		const State& right = m_state[1];
		const double c_left = SoundSpeed(left);
		const double c_right = SoundSpeed(right);

		const double rho_bar = 0.5 * (left.rho + right.rho);
		const double c_bar = 0.5 * (c_left + c_right);
		const double p_bar = 0.5 * (left.p + right.p);
		const double du = right.u - left.u;

		auto total = [&](double p) -> double
		{
			return PressureFunction(p, left) + PressureFunction(p, right) + du;
		};
		auto derivative = [&](double p) -> double
		{
			return PressureFunctionDerivative(p, left) + PressureFunctionDerivative(p, right);
		};

		// --- 括号：f(0+) < 0 < f(hi)，必要时放大 hi ---
		double lo = 0.0;
		double hi = std::max(p_bar, std::max(left.p, right.p)) * 2.0 + 1.0;
		int guard = 0;
		while (total(hi) < 0.0 && guard++ < 200)
		{
			hi *= 2.0;
		}

		// --- 初值：PVRS (Toro 4.53)；非物理时退回 TRBS (4.55)；再不行取括号中点 ---
		double p = p_bar - 0.5 * du * rho_bar * c_bar;
		if (!(p > lo && p < hi))
		{
			const double numerator = c_bar - 0.5 * (m_gamma - 1.0) * du;
			const double denominator = c_bar / std::pow(p_bar, (m_gamma - 1.0) / (2.0 * m_gamma));
			if (numerator > 0.0 && denominator > 0.0)
			{
				p = std::pow(numerator / denominator, 2.0 * m_gamma / (m_gamma - 1.0));
			}
			if (!(p > lo && p < hi) || !std::isfinite(p))
			{
				p = 0.5 * (lo + hi);
			}
		}

		m_newton_res = 0.0;
		for (int iter = 0; iter < kMaxNewtonIter; ++iter)
		{
			const double f = total(p);
			const double df = derivative(p);
			m_newton_iter = iter + 1;
			m_newton_res = std::fabs(f);

			if (f < 0.0)
			{
				lo = std::max(lo, p);
			}
			else
			{
				hi = std::min(hi, p);
			}

			double p_new = p - f / df;
			if (!std::isfinite(p_new) || p_new <= lo || p_new >= hi)
			{
				p_new = 0.5 * (lo + hi);
			}
			const double change = std::fabs(p_new - p)
				/ std::max(0.5 * (p_new + p), kMinPressure);
			p = p_new;
			if (change < kTolerance)
			{
				break;
			}
		}

		m_p_star = p;
		m_shock[0] = (m_p_star > left.p);
		m_shock[1] = (m_p_star > right.p);
		// Toro (4.51)
		m_u_star = 0.5 * (left.u + right.u)
			+ 0.5 * (PressureFunction(m_p_star, right) - PressureFunction(m_p_star, left));
	}

	// ------------------------------------------------------------------
	// 采样
	// ------------------------------------------------------------------
	EulerExactRiemann::State EulerExactRiemann::Sample(double xi) const
	{
		const int side = (xi <= m_u_star) ? 0 : 1;
		const State& s = m_state[side];

		if (m_shock[side])
		{
			const double shock_speed = ShockSpeed(side);
			// 左侧：xi < S 时仍是初始态；右侧：xi > S 时仍是初始态
			if ((side == 0 && xi < shock_speed) || (side == 1 && xi > shock_speed))
			{
				return s;
			}
			return State(StarDensity(side), m_u_star, m_p_star);
		}

		double head = 0.0, tail = 0.0, c_star = 0.0;
		RarefactionSpeed(side, head, tail, c_star);

		if (side == 0 && xi < head)
		{
			return s;
		}
		if (side == 1 && xi > head)
		{
			return s;
		}
		if (side == 0 && xi > tail)
		{
			return State(StarDensity(side), m_u_star, m_p_star);
		}
		if (side == 1 && xi < tail)
		{
			return State(StarDensity(side), m_u_star, m_p_star);
		}

		// 稀疏扇内部，Toro (4.56)
		const double c = SoundSpeed(s);
		const double gm = m_gamma - 1.0;
		const double gp = m_gamma + 1.0;
		double u = 0.0;
		double c_local = 0.0;
		if (side == 0)
		{
			u = 2.0 / gp * (c + 0.5 * gm * s.u + xi);
			c_local = 2.0 / gp * (c + 0.5 * gm * (s.u - xi));
		}
		else
		{
			u = 2.0 / gp * (-c + 0.5 * gm * s.u + xi);
			c_local = 2.0 / gp * (c - 0.5 * gm * (s.u - xi));
		}
		State out;
		out.u = u;
		out.rho = s.rho * std::pow(c_local / c, 2.0 / gm);
		out.p = s.p * std::pow(c_local / c, 2.0 * m_gamma / gm);
		return out;
	}

	EulerExactRiemann::State EulerExactRiemann::Sample(double x, double t) const
	{
		if (t <= 0.0)
		{
			return (x < 0.0) ? m_state[0] : m_state[1];
		}
		return Sample(x / t);
	}
}
