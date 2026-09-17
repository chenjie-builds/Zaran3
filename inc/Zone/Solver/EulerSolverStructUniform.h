/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file EulerSolverStructUniform.h
 * \brief Euler equations solver on a uniform Cartesian structured grid.
 * \author Chen Jie.
 *
 *  设计约定（与项目既有架构保持一致）
 *  ------------------------------------
 *   - 继承自 NSSolver，因此初始化/时间推进/预处理/后处理全部走框架既有流程；
 *   - 流场数据全部由 FieldData / DataManagerNSStruct 管理，变量存于**网格节点**上
 *     （primitive_{0..4}、conservative_{0..4}、residual_{0..4}、dt）；
 *   - 原始变量与守恒变量之间的转换复用 Gas::Prim2Cons / Cons2Prim，
 *     因而无量纲化、状态方程与其它求解器完全一致；
 *   - 界面通量复用框架的 Riemann 求解器（HLLC / Roe / VanLeer / AUSMPW / StegerWarming）；
 *   - 边界条件走 grid->GetBoundMap() 的 inlet / outlet / wall / farfield 约定；
 *   - 结果由 Visual::WriteTecplotBinary 输出 Tecplot，与其它求解器共用同一输出通道。
 *
 *  与 NSSolverStruct 的差别
 *   ------------------------
 *   NSSolverStruct 面向贴体曲线坐标，需要度量系数、雅可比矩阵等；
 *   本求解器限定网格为**均匀笛卡尔网格**，因此：
 *     - CalcCoordTransCoef() 只校验均匀性并记录 dx/dy/dz；
 *     - 雅可比恒为 1，通量差分即普通的中心差；
 *     - 不需要 Metric 相关的任何存储。
 *   这也正是为后续接入 DEM 两相耦合预留的、最简单也最容易验证的流体基座。
 *
 *  体积分数（两相）扩展点
 *   ----------------------
 *   两相版本 EulerTwoPhaseStructUniform 直接继承本类，只通过一个虚钩子接入：
 *
 *       virtual const double* GetVolumeFractionField() const;
 *
 *   单相时返回 nullptr，此时
 *     - 界面通量的体积分数权重恒为 1.0，乘 1.0 在 IEEE 下精确（通量逐位不变）；
 *     - 变量转换的权重恒为 1.0，乘/除 1.0 同样精确（守恒量逐位不变）。
 *   于是"ε ≡ 1 时逐位退化到单相结果"不是靠特判分支，而是由算术本身保证的。
 */

#pragma once
#include "NSSolver.h"
#include "GridStruct.h"
#include "DataManagerNSStruct.h"
#include "FlowSolverParamUniform.h"

namespace zaran
{
	/// @brief 均匀笛卡尔结构网格上的 Euler 方程求解器。
	class EulerSolverStructUniform : public NSSolver
	{
	public:
		EulerSolverStructUniform(index_type index, string name,
			shared_ptr<FlowSolverParam> para, shared_ptr<GridBase> grid,
			shared_ptr<DataManagerNS> data_manager);
		virtual ~EulerSolverStructUniform();

	public:
		DataManagerNSStruct* GetDataManager();
		IdProxyStruct& GetIdxProxy();
		FlowSolverParamUniform* GetPara();
		GridStruct* GetGrid();
		/// @brief 均匀网格步长（由 CalcCoordTransCoef 计算）
		double GetDx() const { return m_dx; }
		double GetDy() const { return m_dy; }
		double GetDz() const { return m_dz; }

	protected:
		// ---------------- NSSolver 接口 ----------------
		/// @brief 初值分派（在基类基础上增加激波管初值）
		void InitField() override;
		void InitFieldFarfield() override;
		void InitFieldFarFieldZeroVel() override;
		void InitFieldBackup() override;
		void InitFieldExplosion() override;
		void InitFieldVortex();
		/// @brief 一维 Riemann（激波管）初值：按 x 方向的膜片位置分两段常值
		void InitFieldRiemann1D();
		/// @brief 均匀网格无需坐标变换，这里只校验均匀性并记录步长
		void CalcCoordTransCoef() override;
		void CalcTimeStepLocal() override;
		void CalcMinTimeStep(double& min_dt) override;
		void RungeKutta() override;
		void Prim2Cons() override;
		void Cons2Prim() override;
		void ZeroResidual() override;
		void CalcConvectionResidual() override;
		void CalcViscousResidual() override {}
		void CalcViscousFlux() override {}
		void CalcViscousFluxGrad() override {}
		void CalcSourceResidual() override {}
		void CheckPrimtive() override;
		void CheckResidual() override;
		void FixPrimtive() override;

		// ---------------- FieldSolver 接口 ----------------
		void BoundaryCondition() override;
		void BackupField(std::string& back_folder) override;

		// ---------------- FlowFieldSolver 接口 ----------------
		void ReduceTimeStep(double& dt) override;

	protected:
		/// @brief 气相体积分数场（节点中心，含 ghost）。
		/// @details 单相求解器返回 nullptr，此时所有体积分数权重恒为 1.0，
		///          乘 1.0 / 除 1.0 在 IEEE 下精确，因而通量与变量转换逐位不变。
		///          两相求解器（EulerTwoPhaseStructUniform）覆盖该函数返回 ε 场。
		virtual const double* GetVolumeFractionField() const { return nullptr; }

	private:
		/// @brief 沿指定方向累加对流通量残差
		void AddFluxResidual(int dir, bool first_order);
		/// @brief 面左右态重构（一阶直接取节点值；二阶用 MUSCL + 限制器）
		void ReconstructFace(const double* const prim[5], index_type idx_m, index_type idx_0,
			index_type idx_p, index_type idx_pp, bool first_order,
			double* prim_left, double* prim_right);
		/// @brief 限制器（复用框架 Limiter 函数库）
		double ApplyLimiter(double slope_left, double slope_right);
		/// @brief 是否需要强制一阶（起步阶段）
		bool UseFirstOrder();

		shared_ptr<GridStruct> m_grid;
		shared_ptr<DataManagerNSStruct> m_data_manager;
		shared_ptr<FlowSolverParamUniform> m_para;
		shared_ptr<IdProxyStruct> m_idx_proxy;
		/// @brief 均匀网格步长
		double m_dx = 0.0, m_dy = 0.0, m_dz = 0.0;
	};
}
