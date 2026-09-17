/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file NSFieldStructTwoPhase.h
 * \brief Two-phase NS field on a uniform Cartesian structured grid.
 * \author Chen Jie.
 *
 *  与 NSFieldStructUniform 的关系
 *  -----------------------------
 *  网格（GridStruct）、残差统计（NSFieldStruct::CalcResidual）、输出通道
 *  （Visual 的 NS_Structured 分支）完全一致，只是把三件东西换成两相版本：
 *    求解器      EulerSolverStructUniform  -> EulerTwoPhaseStructUniform
 *    求解器参数  FlowSolverParamUniform    -> FlowSolverParamTwoPhase
 *    数据管理器  DataManagerNSStruct       -> DataManagerNSTwoPhase（多一个 ε 场）
 */

#pragma once
#include "NSFieldStruct.h"
#include "EulerTwoPhaseStructUniform.h"
#include "FlowSolverParamTwoPhase.h"
#include "DataManagerNSTwoPhase.h"

namespace zaran
{
	/// @brief 两相（体积分数加权）均匀结构网格 Euler 场。
	class NSFieldStructTwoPhase : public NSFieldStruct
	{
	public:
		NSFieldStructTwoPhase(shared_ptr<GridStruct> grid);
		~NSFieldStructTwoPhase();
		/// @brief 取两相求解器
		shared_ptr<EulerTwoPhaseStructUniform> GetTwoPhaseSolver();
		/// @brief 取两相求解器参数
		shared_ptr<FlowSolverParamTwoPhase> GetTwoPhaseSolverPara();
		/// @brief 取两相数据管理器
		shared_ptr<DataManagerNSTwoPhase> GetTwoPhaseDataManager();

	protected:
		void AllocateSolver() override;
		void AllocateSolverPara() override;
		void AllocateDataManager() override;
	};
}
