/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file FlowSolverParamTwoPhase.cpp
 * \brief Solver parameters for the two-phase (volume-fraction) Euler solver.
 */

#include "FlowSolverParamTwoPhase.h"
#include "GlobalData.h"
#include "Log.h"

namespace zaran
{
	namespace
	{
		/// @brief ε 场来源的可读名字（日志里出现 int 枚举值太难看，也容易看错）
		const char* VolumeFractionTypeName(VolumeFractionType type)
		{
			switch (type)
			{
			case VolumeFractionType::Uniform:  return "uniform";
			case VolumeFractionType::Step:     return "step";
			case VolumeFractionType::Sine:     return "sine";
			case VolumeFractionType::External: return "external(由 DEM-CFD 耦合器注入)";
			}
			return "unknown";
		}
	}

	FlowSolverParamTwoPhase::FlowSolverParamTwoPhase()
	{
	}

	FlowSolverParamTwoPhase::~FlowSolverParamTwoPhase()
	{
	}

	void FlowSolverParamTwoPhase::Init()
	{
		// 基类负责 CFL / RK / 限制器 / Riemann / 初值类型 / 重构阶数 / 网格校验
		FlowSolverParamUniform::Init();

		// ---------------- [init.volume_fraction] ----------------
		if (GlobalData::IsExist("init.volume_fraction.type"))
		{
			const std::string type = GlobalData::GetString("init.volume_fraction.type");
			if (type == "uniform")
			{
				m_volume_fraction_type = VolumeFractionType::Uniform;
			}
			else if (type == "step")
			{
				m_volume_fraction_type = VolumeFractionType::Step;
			}
			else if (type == "sine")
			{
				m_volume_fraction_type = VolumeFractionType::Sine;
			}
			else if (type == "external")
			{
				m_volume_fraction_type = VolumeFractionType::External;
			}
			else
			{
				Log::warn("FlowSolverParamTwoPhase: unknown volume fraction type '{}', "
					"fallback to uniform.", type);
				m_volume_fraction_type = VolumeFractionType::Uniform;
			}
		}
		if (GlobalData::IsExist("init.volume_fraction.value"))
		{
			m_volume_fraction_value = GlobalData::GetDouble("init.volume_fraction.value");
		}
		if (GlobalData::IsExist("init.volume_fraction.x_diaphragm"))
		{
			m_volume_fraction_x_step = GlobalData::GetDouble("init.volume_fraction.x_diaphragm");
		}
		if (GlobalData::IsExist("init.volume_fraction.left_value"))
		{
			m_volume_fraction_left = GlobalData::GetDouble("init.volume_fraction.left_value");
		}
		if (GlobalData::IsExist("init.volume_fraction.right_value"))
		{
			m_volume_fraction_right = GlobalData::GetDouble("init.volume_fraction.right_value");
		}
		if (GlobalData::IsExist("init.volume_fraction.mean"))
		{
			m_volume_fraction_mean = GlobalData::GetDouble("init.volume_fraction.mean");
		}
		if (GlobalData::IsExist("init.volume_fraction.amplitude"))
		{
			m_volume_fraction_amplitude = GlobalData::GetDouble("init.volume_fraction.amplitude");
		}
		if (GlobalData::IsExist("init.volume_fraction.wavelength"))
		{
			m_volume_fraction_wavelength = GlobalData::GetDouble("init.volume_fraction.wavelength");
		}
		if (GlobalData::IsExist("init.volume_fraction.phase"))
		{
			m_volume_fraction_phase = GlobalData::GetDouble("init.volume_fraction.phase");
		}

		// ---------------- [two_phase] ----------------
		if (GlobalData::IsExist("two_phase.porosity_gradient_force"))
		{
			m_porosity_gradient_force = GlobalData::GetInt("two_phase.porosity_gradient_force");
		}
		if (GlobalData::IsExist("two_phase.volume_fraction_min"))
		{
			m_volume_fraction_min = GlobalData::GetDouble("two_phase.volume_fraction_min");
		}

		if (m_volume_fraction_wavelength <= 0.0)
		{
			Log::warn("FlowSolverParamTwoPhase: non-positive wavelength {}, fallback to 1.0",
				m_volume_fraction_wavelength);
			m_volume_fraction_wavelength = 1.0;
		}

		Log::info("Two-phase Euler parameters: volume_fraction_type={}, value={}, "
			"step=({}, {} -> {}), sine=(mean {}, amp {}, lambda {}, phase {}), "
			"porosity_gradient_force={}",
			VolumeFractionTypeName(m_volume_fraction_type), m_volume_fraction_value,
			m_volume_fraction_x_step, m_volume_fraction_left, m_volume_fraction_right,
			m_volume_fraction_mean, m_volume_fraction_amplitude,
			m_volume_fraction_wavelength, m_volume_fraction_phase,
			m_porosity_gradient_force);
	}
}
