/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file EulerTwoPhaseStructUniform.h
 * \brief Two-phase (volume-fraction weighted) Euler solver on a uniform Cartesian grid.
 * \author Chen Jie.
 *
 *  阶段 1 的目标（见 docs/TWO_PHASE_EULER_STAGE1.md）
 *  ------------------------------------------------
 *  在均匀结构网格上建立**气相体积分数加权**的可压缩 Euler 求解器，
 *  但暂时不含颗粒：ε_g(x) 是给定的冻结场（阶段 2 起由 DEM 的粒子→网格
 *  体积分配写入同一个数组）。
 *
 *  离散方程（ε 冻结 ⇒ ∂ε/∂t = 0）
 *  ------------------------------
 *      ∂(ε ρ)   /∂t + ∇·(ε ρ u)              = 0
 *      ∂(ε ρ u) /∂t + ∇·(ε ρ u u + ε p I)    = p ∇ε  (+ 相间力，阶段 3 起)
 *      ∂(ε E)   /∂t + ∇·(ε (E + p) u)        = 0     (+ 相间功/热，阶段 3 起)
 *
 *  - 守恒量按"每单位总体积的气相量"存：cons = ε·cons_intrinsic；
 *  - 界面通量 = ε_face · F_Godunov(气相原始变量)，ε_face 取两侧节点算术平均；
 *  - 动量方程把 ∇·(ε p I) 展开出一项 p∇ε 作为源项，于是左侧仍是纯 Godunov 通量，
 *    而等效于 DEM-CFD 惯用的 ε∇p（与颗粒反作用力 -V∇p 天然配套）。
 *
 *  为什么"ε ≡ 1 ⇒ 逐位退化"是可证的
 *  -------------------------------
 *    ε ≡ 1 时 ε_face = 0.5(1+1) = 1.0、权重 = 1.0，乘/除 1.0 在 IEEE 下是精确运算；
 *    而 p∇ε 用中心差分得到 (1-1)/(2dx) = 0.0，加到残差上也是精确的空操作。
 *    因此两相求解器在 ε ≡ 1 下的每一步都与单相求解器逐位相同——这不是特判，
 *    而是算术本身的性质。verify_two_phase.py 的第 1 项检验就是对这个论断的回归。
 */

#pragma once
#include "EulerSolverStructUniform.h"
#include "DataManagerNSTwoPhase.h"
#include "FlowSolverParamTwoPhase.h"

namespace zaran
{
	/// @brief 带气相体积分数加权的均匀结构网格 Euler 求解器。
	class EulerTwoPhaseStructUniform : public EulerSolverStructUniform
	{
	public:
		EulerTwoPhaseStructUniform(index_type index, string name,
			shared_ptr<FlowSolverParam> para, shared_ptr<GridBase> grid,
			shared_ptr<DataManagerNS> data_manager);
		virtual ~EulerTwoPhaseStructUniform();

	public:
		DataManagerNSTwoPhase* GetTwoPhaseDataManager();
		FlowSolverParamTwoPhase* GetTwoPhasePara();
		/// @brief 气相体积分数场（节点中心，含 ghost）
		const double* GetVolumeFraction() const;

		/// @brief 注入相间反作用力源项（DEM-CFD 双向耦合，**无量纲**，节点量，长度 3*N）
		/// @details 约定：`src[3*idx + dir]` 是加在节点 idx 的第 dir 个动量方程上的
		///          **单位体积力**（无量纲形式 `F*·L_ref/(ρ_ref a_ref²)`）。
		///          传 nullptr 关闭。每步由耦合器按"颗粒受力的负值"经同一套网格权重散列得到，
		///          于是 `Σ_cell src·V_cell = −Σ_particle F*`（牛顿第三定律）是**构造保证**的。
		///          ⚠ 只反作用**曳力**：颗粒的压力梯度力 −V_p∇p 是气相方程里 −ε_g∇p 的
		///          对偶项（见 docs/TWO_PHASE_STAGE3/4），再反作用一次就是双计。
		void SetInterphaseMomentumSource(const double* src) { m_interphase_source = src; }
		/// @brief 相间力是否对气相**做功能**（把 S·u 加进能量方程）
		/// @details 严谨的写法里相间力在能量方程里应当有 S·u 这一项（相间力对气相做的功）。
		///          打开后能量账本闭合：`d(∭εE)/dt = ∭S·u dV = −Σ_p F_drag·u_g`，
		///          且 `d(∭εE + Σ½m_p|u_p|²)/dt = −耗散 + 浮力功`（可据此做严格的能量收支检验）。
		///          不少 CFD-DEM 代码省略这一项（只把力加进动量方程），那样气相会"白拿动量"。
		///          默认打开。
		void SetInterphaseWork(bool enabled) { m_interphase_work = enabled; }
		/// @brief 当前是否接了相间源项
		bool HasInterphaseMomentumSource() const { return m_interphase_source != nullptr; }

	protected:
		// ---------------- 求解器接入点 ----------------
		/// @brief 体积分数场钩子：本类返回 ε 场，基类据此对通量与变量转换加权
		const double* GetVolumeFractionField() const override;
		/// @brief 先铺 ε 场，再交给基类做流场初值（基类末尾会调用虚函数 Prim2Cons）
		void InitField() override;
		/// @brief 基类边界条件 + ε 的零梯度虚拟节点
		void BoundaryCondition() override;
		/// @brief 基类对流通量残差 + 孔隙率梯度动量源 p∇ε
		void CalcConvectionResidual() override;
		/// @brief 基类物理性检查 + ε 范围检查
		void CheckPrimtive() override;
		/// @brief 相间反作用力源项（双向耦合）；无源项时与基类的空实现等价
		void CalcSourceResidual() override;

	private:
		/// @brief 按参数类型铺 ε 场（全部节点，含 ghost）
		void InitVolumeFraction();
		/// @brief 在参数给定的节点坐标处求 ε
		double EvaluateVolumeFraction(double x, double y, double z) const;
		/// @brief ε 的虚拟节点：零梯度外推（冻结孔隙率的自然取法）
		void FillVolumeFractionGhost();
		/// @brief 动量残差附加 +p∇ε（仅当 two_phase.porosity_gradient_force = 1）
		void AddPorosityGradientSource();

		shared_ptr<DataManagerNSTwoPhase> m_tp_data_manager;
		shared_ptr<FlowSolverParamTwoPhase> m_tp_para;
		/// @brief 相间反作用力源项（3*N，无量纲"单位体积力"），nullptr = 未接
		const double* m_interphase_source = nullptr;
		/// @brief 相间力是否对气相做功能（把 S·u 加进能量方程，默认开）
		bool m_interphase_work = true;
	};
}
