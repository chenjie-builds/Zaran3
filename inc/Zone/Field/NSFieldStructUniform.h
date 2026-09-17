/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file NSFieldStructUniform.h
 * \brief NS field on a uniform Cartesian structured grid.
 * \author Chen Jie.
 */

#pragma once
#include "NSFieldStruct.h"
#include "EulerSolverStructUniform.h"
#include "FlowSolverParamUniform.h"
#include "DataManagerNSStruct.h"

namespace zaran
{
	/// @brief 均匀结构网格 Euler 场。
	/// @details 与 NSFieldStruct 共用同一套网格（GridStruct）、数据管理（DataManagerNSStruct）、
	///          残差统计（FieldNS::CalcResidual）与输出通道（Visual 的 NS_Structured 分支），
	///          只是把求解器与求解器参数换成均匀网格版本。这样 Tecplot 输出、残差文件、
	///          多块交接面通信等框架能力都可以直接复用。
	class NSFieldStructUniform : public NSFieldStruct
	{
	public:
		NSFieldStructUniform(shared_ptr<GridStruct> grid);
		~NSFieldStructUniform();
		/// @brief 取均匀网格 Euler 求解器
		shared_ptr<EulerSolverStructUniform> GetSolver();
		/// @brief 取均匀网格求解器参数
		shared_ptr<FlowSolverParamUniform> GetSolverPara();

	protected:
		void AllocateSolver() override;
		void AllocateSolverPara() override;
	};
}
