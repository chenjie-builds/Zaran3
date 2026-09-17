/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file FlowSolverParamUniform.cpp
 * \brief Solver parameters for the uniform Cartesian structured-grid Euler solver.
 */

#include "FlowSolverParamUniform.h"
#include "GlobalData.h"
#include "Log.h"

namespace zaran
{
	FlowSolverParamUniform::FlowSolverParamUniform()
	{
	}

	FlowSolverParamUniform::~FlowSolverParamUniform()
	{
	}

	void FlowSolverParamUniform::Init()
	{
		// 基类负责 CFL、RK 系数、限制器、Riemann 求解器、来流条件、初值类型、
		// 以及 [structure] 节的度量参数（均匀网格不使用，但保证控制文件结构一致）
		FlowSolverParamStruct::Init();

		if (GlobalData::IsExist("space.order"))
		{
			m_recon_order = GlobalData::GetInt("space.order");
		}
		if (GlobalData::IsExist("space.require_uniform_grid"))
		{
			m_require_uniform = GlobalData::GetInt("space.require_uniform_grid");
		}
		if (GlobalData::IsExist("space.first_order_steps"))
		{
			m_first_order_steps = GlobalData::GetInt("space.first_order_steps");
		}

		if (m_recon_order != 1 && m_recon_order != 2)
		{
			Log::warn("FlowSolverParamUniform: unsupported reconstruction order {}, fallback to 2",
				m_recon_order);
			m_recon_order = 2;
		}

		Log::info("Euler uniform solver parameters: reconstruction_order={}, require_uniform_grid={}, "
			"limiter={}, cfl={}",
			m_recon_order, m_require_uniform,
			static_cast<int>(GetLimiterType()), GetCflNumber());
	}
}
