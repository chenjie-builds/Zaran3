/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file DragForceModel.cpp
 * \brief Inter-phase force models (fluid -> particle) for the one-way coupled DEM-CFD solver.
 */

#include "DragForceModel.h"
#include "ZaranError.h"

#include <algorithm>
#include <cctype>
#include <cmath>

namespace zaran
{
	namespace
	{
		constexpr double kPi = 3.14159265358979323846;

		/// @brief 把模型名转成小写并去掉 '-' 与 '_'（"Wen-Yu" / "wen_yu" / "wenyu" 等价）
		std::string Normalize(const std::string& name)
		{
			std::string out;
			out.reserve(name.size());
			for (const char c : name)
			{
				if (c == '-' || c == '_' || c == ' ')
				{
					continue;
				}
				out.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
			}
			return out;
		}
	}

	bool DragForceModel::ParseModel(const std::string& name, Model& model)
	{
		const std::string key = Normalize(name);
		if (key == "stokes")
		{
			model = Model::Stokes;
			return true;
		}
		if (key == "schillernaumann" || key == "sn" || key == "singleparticle")
		{
			model = Model::SchillerNaumann;
			return true;
		}
		if (key == "wenyu")
		{
			model = Model::WenYu;
			return true;
		}
		if (key == "ergun")
		{
			model = Model::Ergun;
			return true;
		}
		if (key == "gidaspow")
		{
			model = Model::Gidaspow;
			return true;
		}
		return false;
	}

	const char* DragForceModel::ModelName(Model model)
	{
		switch (model)
		{
		case Model::Stokes:            return "stokes";
		case Model::SchillerNaumann:   return "schiller_naumann";
		case Model::WenYu:             return "wen_yu";
		case Model::Ergun:             return "ergun";
		case Model::Gidaspow:          return "gidaspow";
		}
		return "unknown";
	}

	bool DragForceModel::IsDenseModel(Model model)
	{
		return model == Model::Ergun;
	}

	double DragForceModel::SchillerNaumannCd(double re)
	{
		// Re ≤ 0（相对速度为零或非法）⇒ 曳力为零，C_D 记为 0（若返回 24/Re 会得到 inf，
		// 写进 CSV 后整个文件不可用）。这一约定必须与 verify.py 侧保持一致。
		if (!(re > 0.0))
		{
			return 0.0;
		}
		if (re < 1000.0)
		{
			return (24.0 / re) * (1.0 + 0.15 * std::pow(re, 0.687));
		}
		return 0.44;
	}

	double DragForceModel::DragCoefficientPerVolume(Model model, const Options& opts,
		double eps_gas, double gas_density, double gas_viscosity,
		double diameter, double speed, double& re, double& cd, bool& dense_used)
	{
		re = 0.0;
		cd = 0.0;
		dense_used = false;

		if (!(diameter > 0.0) || !(gas_viscosity > 0.0) || !(gas_density > 0.0))
		{
			return 0.0;
		}
		// 空隙率做一次硬性限制，避免 ε_g → 0 时 ε_g^(−2.65) 溢出
		const double eps = std::min(1.0, std::max(1.0e-6, eps_gas));
		const double eps_s = std::max(0.0, 1.0 - eps);

		// ---------------- 选择实际使用的模型（Gidaspow 分段） ----------------
		Model used = model;
		if (model == Model::Gidaspow)
		{
			used = (eps >= opts.dense_threshold) ? Model::WenYu : Model::Ergun;
		}
		dense_used = (used == Model::Ergun);

		// ---------------- Stokes：β/ε_s = 18μ/d²，对 speed = 0 良定义 ----------------
		if (used == Model::Stokes)
		{
			re = gas_density * diameter * speed / gas_viscosity;
			// 与 Schiller–Naumann 在 Re → 0 的极限一致：C_D = 24/Re
			cd = (re > 0.0) ? 24.0 / re : 0.0;
			return 18.0 * gas_viscosity / (diameter * diameter);
		}

		// ---------------- Ergun（稠密支） ----------------
		if (used == Model::Ergun)
		{
			re = gas_density * diameter * speed / gas_viscosity;
			cd = 0.0;
			const double linear = 150.0 * gas_viscosity * eps_s / (eps * diameter * diameter);
			const double quadratic = 1.75 * gas_density * speed / diameter;
			return linear + quadratic;
		}

		// ---------------- Schiller–Naumann 族（稀疏/单颗粒支） ----------------
		const double re_density = (used == Model::WenYu && opts.dilute_re_uses_voidage)
			? eps * gas_density : gas_density;
		re = re_density * diameter * speed / gas_viscosity;
		cd = SchillerNaumannCd(re);

		double voidage_factor = 1.0;
		if (used == Model::WenYu)
		{
			voidage_factor = std::pow(eps, -opts.voidage_exponent);
		}
		return 0.75 * cd * gas_density * speed * voidage_factor / diameter;
	}

	DragForceModel::Result DragForceModel::Evaluate(Model model, const Options& opts,
		const GasState& gas, const ParticleState& particle,
		double eps_gas, const double grad_p[3])
	{
		Result out;
		out.eps_gas = eps_gas;

		if (!(particle.radius > 0.0))
		{
			return out;
		}
		const double diameter = 2.0 * particle.radius;
		const double volume = (4.0 / 3.0) * kPi * particle.radius
			* particle.radius * particle.radius;

		for (int c = 0; c < 3; ++c)
		{
			out.relative_velocity[c] = gas.velocity[c] - particle.velocity[c];
		}
		out.relative_speed = std::sqrt(out.relative_velocity[0] * out.relative_velocity[0]
			+ out.relative_velocity[1] * out.relative_velocity[1]
			+ out.relative_velocity[2] * out.relative_velocity[2]);

		double re = 0.0, cd = 0.0;
		bool dense = false;
		out.beta_over_es = DragCoefficientPerVolume(model, opts, eps_gas,
			gas.density, gas.viscosity, diameter, out.relative_speed, re, cd, dense);
		out.re = re;
		out.cd = cd;
		out.branch = dense ? 1 : 0;

		const double k = volume * out.beta_over_es;
		for (int c = 0; c < 3; ++c)
		{
			out.drag[c] = k * out.relative_velocity[c];
		}

		if (opts.include_pressure_gradient_force && grad_p != nullptr)
		{
			for (int c = 0; c < 3; ++c)
			{
				out.pressure_gradient[c] = -volume * grad_p[c];
			}
		}

		for (int c = 0; c < 3; ++c)
		{
			out.force[c] = out.drag[c] + out.pressure_gradient[c];
		}
		return out;
	}
}
