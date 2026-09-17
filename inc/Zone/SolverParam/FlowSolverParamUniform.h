/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file FlowSolverParamUniform.h
 * \brief Solver parameters for the uniform Cartesian structured-grid Euler solver.
 * \author Chen Jie.
 */

#pragma once
#include "FlowSolverStructPara.h"

namespace zaran
{
	/// @brief 均匀笛卡尔结构网格 Euler 求解器的参数。
	/// @details 继承 FlowSolverParamStruct —— 这样 FieldNS / NSFieldStruct 内部对参数做的
	///          static_pointer_cast<FlowSolverParamStruct> 仍然成立，不会引入类型欺骗；
	///          同时 [structure] 节的度量参数虽然对均匀网格不用，但保证控制文件结构统一。
	///          CFL、RK 系数、限制器、Riemann 求解器、来流条件等全部沿用基类。
	class FlowSolverParamUniform : public FlowSolverParamStruct
	{
	public:
		FlowSolverParamUniform();
		virtual ~FlowSolverParamUniform();
		/// @brief 从 GlobalData 读取参数
		void Init() override;

	public:
		/// @brief 重构阶数：1 = Godunov 一阶，2 = MUSCL 二阶
		const int& GetReconstructionOrder() const { return m_recon_order; }
		/// @brief 是否要求网格严格均匀（读取网格时做校验）
		const int& GetRequireUniformGrid() const { return m_require_uniform; }
		/// @brief 先用一阶格式推进的步数（强间断算例的稳定启动，0 = 不启用）
		const int& GetFirstOrderSteps() const { return m_first_order_steps; }

	private:
		/// @brief 重构阶数（控制文件 [space] order）
		int m_recon_order = 2;
		/// @brief 是否校验网格均匀性（控制文件 [space] require_uniform_grid）
		int m_require_uniform = 1;
		/// @brief 一阶起步步数（控制文件 [space] first_order_steps）
		int m_first_order_steps = 0;
	};
}
