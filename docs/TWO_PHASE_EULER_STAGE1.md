# 两相可压缩 Euler 求解器（阶段 1）—— 实现与验证

> 阶段 1 的目标：在均匀结构网格上建立**气相体积分数加权**的可压缩 Euler 求解器，
> **先不含颗粒**；验收标准是 **ε ≡ 1 时必须退化到单相求解器的结果**。

- 求解器：`EulerTwoPhaseStructUniform`（`inc/Zone/Solver/` + `src/Zone/Solver/`）
- 参数类：`FlowSolverParamTwoPhase`
- 数据管理器：`DataManagerNSTwoPhase`（比单相多注册一个节点量 `volume_fraction`）
- 场类：`NSFieldStructTwoPhase`
- 基准驱动：`benchmarks/euler_uniform/{make_euler_case.py, verify_two_phase.py}`
- 验证算例：`_eulerbench/two_phase/`（脚本自动生成，不入版本库）

---

## 1. 方程形式

气相体积分数 ε_g(x,t) 冻结（阶段 1 为由控制文件规定的场，阶段 2 起由 DEM 的
粒子→网格体积分配写入同一数组），因此 ∂ε/∂t = 0。取体积平均的两相模型：

```
∂(ε ρ)   /∂t + ∇·(ε ρ u)             = 0
∂(ε ρ u) /∂t + ∇·(ε ρ u u + ε p I)   = p ∇ε          (+ 相间力，阶段 3 起)
∂(ε E)   /∂t + ∇·(ε (E + p) u)       = 0             (+ 相间功/热，阶段 3 起)
```

其中 ρ, u, p, E 都是**气相自身**的量。这条写法有三个好处：

1. **守恒量是"每单位总体积的气相量"** `εU`，界面通量是纯 Godunov 通量乘一个标量权重，
   Riemann 求解器原封不动复用；
2. **与 DEM-CFD 的反作用天然配套**：颗粒受到的压力梯度力是 `-V_p ∇p`，气相就应当感受到
   `-ε∇p`；把 `∇·(ε p I)` 展开就得到 `ε∇p + p∇ε`，于是把 `p∇ε` 移到右端正是上式；
3. **能量方程不需要任何修正项**，因为 `ε(ρE + p)u` 本身就是"能量通量 × 面积份额"。

### 1.1 与"完全守恒形式"的差别，以及为什么默认选 ε∇p

若把 `p∇ε` 也留在左端，就得到完全守恒的 `∂(εU)/∂t + ∇·(εF) = S`。两种写法在
**ε 均匀**（特别是 ε ≡ 1）时完全等价，ε 变化时相差一项 `p∇ε`：

| 约定 | 控制文件 | 特点 |
|---|---|---|
| **ε∇p**（默认） | `two_phase.porosity_gradient_force = 1` | 与颗粒的 `-V∇p` 力配套；**有离散精确的静水平衡解** |
| ∇(εp)（守恒） | `two_phase.porosity_gradient_force = 0` | 动量严格守恒，但缺少 `p∇ε`；均匀压强下静止气体也会被"ε 梯度"推起来 |

验证脚本 T3 会把两种约定都跑一遍，用来证明 ε∇p 形式的平衡性是真实的（不是"残差本来就小"）。

---

## 2. 离散

### 2.1 守恒量与原始变量的对应

原始变量**不变**：`primitive_{0..4} = (ρ, u, v, w, p)`，与单相求解器、与状态方程
`Gas::Prim2Cons`、与 Tecplot 输出完全一致。守恒量按总体积加权：

```
cons = ε · Gas::Prim2Cons(prim)          (Prim2Cons)
prim = Gas::Cons2Prim(cons / ε)          (Cons2Prim)
```

### 2.2 界面通量与残差

```
F_face = ε_face · F_Godunov(prim_L^recon, prim_R^recon),   ε_face = ½(ε_L + ε_R)
d(εU)_i/dt = -[ F_face(i+½) - F_face(i-½) ] / Δx  +  δ_{d,动量} · p_i [∇ε]_i
[∇ε]_i = (ε_{i+1} - ε_{i-1}) / (2 Δx)
```

ε 的界面权重取**算术平均**、梯度取**中心差分**，二者配对是有原因的，见下。

### 2.3 静水平衡（well-balanced）的证明

均匀静止态（u = 0、ρ、p 均为常数）下，每个界面上的 Riemann 通量都是
`F = (0, p, 0, 0, 0)`（HLLC 在左右态相同时逐位返回物理通量）。于是第 i 个节点的
动量残差为

```
res_i = -p [ ½(ε_i+ε_{i+1}) - ½(ε_{i-1}+ε_i) ] / Δx + p (ε_{i+1}-ε_{i-1})/(2Δx)
      = -p (ε_{i+1}-ε_{i-1})/(2Δx) + p (ε_{i+1}-ε_{i-1})/(2Δx)
      = 0
```

**与 ε 的具体分布无关** —— 任意 ε(x) 下"均匀压强 + 静止"都是（精确算术下的）定常解。
若把界面权重换成上下游取值、或把梯度换成单侧差分、或把源项符号写反，这个抵消立刻失效，
残差会变成 `O(p|∇ε|)` 量级，速度在几步内涨到 O(1)。这正是 T3 的判据。
（浮点上两者并不逐位相等，残差在 `~1e-16·p/Δx` 的舍入水平。）

### 2.4 为什么 "ε ≡ 1 ⇒ 逐位退化" 是可证的

ε ≡ 1 时

- `ε_face = ½(1.0 + 1.0) = 1.0`，通量乘 1.0 在 IEEE 下**精确**；
- `[∇ε]_i = (1.0 - 1.0)/(2Δx) = 0.0`，`p·0.0 = 0.0`，加到残差上是**精确的空操作**；
- `cons = 1.0 · cons_intrinsic`、`prim = cons / 1.0` 同样**精确**；
- CFL 判据与单相完全相同（ε 加权不改变系统的特征速度 `u, u±c`）。

因此两相求解器在 ε ≡ 1 下的每一步都与单相求解器逐位相同。这一点**不是靠特判分支实现的**，
而是由算术本身保证的：`EulerSolverStructUniform` 里只有一个虚钩子

```cpp
virtual const double* GetVolumeFractionField() const { return nullptr; }   // 单相返回 nullptr
```

返回 `nullptr` 时权重恒为 `1.0`。这与 `NSSolverStruct`→`EulerSolverStructUniform` 的
"复用既有框架"原则一致：单相求解器的数值路径**一行都没有改**。

### 2.5 边界与虚拟节点

- 流场：沿用单相求解器对 `inlet / outlet / wall` 的处理（`grid->GetBoundMap()`）；
- ε：**零梯度外推**（`ε_ghost = ε_边界节点`），冻结孔隙率的自然取法；
- 注意 well-balanced 的抵消**不依赖** ε 虚拟节点的取值（见 2.3），所以边界附近不会额外引入误差。

### 2.6 时间步

`CalcTimeStepLocal` 沿用单相的 `dt = CFL / Σ_d (|u_d| + c)/Δ_d`。ε 加权方程的
特征速度仍是 `u` 与 `u ± c`，与 ε 无关，所以不需要为 ε 缩小步长。
（`p∇ε` 是"几何"源项，不引入新的波速。）

### 2.7 CFL 的适用边界

需要提醒：ε 变化剧烈且 ε 很小时，ε∇p 形式相当于在方程右端引入 `F·∇ε/ε` 量级的源项，
此时源的刚度可能要求更保守的 CFL。阶段 1 的算例取 ε ∈ [0.6, 1.0]，未触及这一区域；
接入颗粒后若局部 ε 降到 0.4 附近，建议重做一次步长敏感性扫描。

---

## 3. 代码结构与挂载点

```
Solver      → FieldSolver → FlowFieldSolver → NSSolver → EulerSolverStructUniform
                                                          └── EulerTwoPhaseStructUniform
Field       → FieldNS → NSFieldStruct → NSFieldStructUniform
                                      └── NSFieldStructTwoPhase
SolverParam → SolverParam → FlowSolverParam → FlowSolverParamStruct → FlowSolverParamUniform
                                                                      └── FlowSolverParamTwoPhase
DataManager → DataManager → DataManagerNS → DataManagerNSStruct → DataManagerNSTwoPhase
```

| 挂载点 | 内容 |
|---|---|
| `FieldSolverType` | 新增 `Euler_TwoPhase_Uniform` |
| `FieldGenerator::Create()` | 新增分支 → `NSFieldStructTwoPhase` |
| `Application::SolveField()` | 新增 `task.solver = "EulerTwoPhase"` 分支（要求 `grid_type = "Structured"`） |
| `Visual` | 场数据的 `FieldData` 含 `volume_fraction` 时自动追加输出；单相输出格式**不变** |
| `InitFieldType` | 复用 `Riemann1D` / `Vortex` 等既有初值类型 |

**新增文件**

```
inc/Zone/Solver/EulerTwoPhaseStructUniform.h        src/Zone/Solver/EulerTwoPhaseStructUniform.cpp
inc/Zone/SolverParam/FlowSolverParamTwoPhase.h      src/Zone/SolverParam/FlowSolverParamTwoPhase.cpp
inc/Zone/DataManager/DataManagerNSTwoPhase.h        src/Zone/DataManager/DataManagerNSTwoPhase.cpp
inc/Zone/Field/NSFieldStructTwoPhase.h              src/Zone/Field/NSFieldStructTwoPhase.cpp
benchmarks/euler_uniform/verify_two_phase.py
```

**修改文件（都是加法式改动）**

```
inc/Zone/Solver/EulerSolverStructUniform.h    新增虚钩子 GetVolumeFractionField()
src/Zone/Solver/EulerSolverStructUniform.cpp  Prim2Cons/Cons2Prim/AddFluxResidual 三处乘 ε 权重
inc/Zone/Solver/FieldSolver.h                 FieldSolverType 新增枚举
src/Generator/Field/FieldGenerator.cpp        新增分支
src/Main/Application.cpp                      新增 "EulerTwoPhase" 分支
inc/Basic/PostProcess/Visual.h                新增带体积分数的重载
src/Basic/PostProcess/Visual.cpp              条件输出 volume_fraction；ASCII 精度提到 17 位
benchmarks/euler_uniform/make_euler_case.py   --two-phase / --vf-* / --porosity-gradient-force / --t-end / --max-iter
```

输出格式：Tecplot ASCII 的 `VARIABLES` 行在两相场下多一项 `Volume_fraction`；
二进制 `.plt` 的变量表头同步增加 `volume_fraction`（只在存在该场时才增加，
单相算例的 `.plt` 表头保持 9 个变量不变）。

---

## 4. 使用方式

```bash
cmake --build build-ninja

# 生成一个两相算例（ε 沿 x 正弦变化）
python benchmarks/euler_uniform/make_euler_case.py \
    --case sod --out tests/euler_two_phase_sod \
    --two-phase --vf-type sine --vf-mean 0.8 --vf-amplitude 0.2 \
    --tecplot-ascii

bin/Release/Zaran3.10.2.exe tests/euler_two_phase_sod

# 四项验收
python benchmarks/euler_uniform/verify_two_phase.py \
    --exe bin/Release/Zaran3.10.2.exe [--quick] [--only t1 t3]
```

---

## 5. 验证结果

四项检验全部通过（`verify_two_phase.py` 默认网格；`--quick` 为粗网格冒烟）。

### T1 ε ≡ 1 逐位退化到单相求解器

同网格、同阶数、同限制器、同 CFL、同结束时间，只把 `task.solver` 从 `"Euler"` 换成
`"EulerTwoPhase"` 并把 ε 设为常值 1，比较最终帧的 Tecplot ASCII 输出（17 位有效数字）：

| 算例 | 阶数 | 网格 | max Δρ/ρ | max Δp/p | max Δu/u |
|---|---|---|---|---|---|
| Sod（含激波/接触/稀疏波） | 2 | 400 | **0.000e+00** | **0.000e+00** | **0.000e+00** |
| Sod | 1 | 400 | **0.000e+00** | **0.000e+00** | **0.000e+00** |
| 等熵涡（二维光滑） | 2 | 81×81 | **0.000e+00** | **0.000e+00** | **0.000e+00** |

五个物理量**逐位相同**（不是"在容差内"，是二进制相同），输出的 `Volume_fraction` 列恒为 1。
这与 2.4 节的论证一致：退化不是靠分支实现的，而是乘/除 1.0 在 IEEE 下精确。

### T2 均匀 ε 不变性

ε 为常值时方程右端的 ε 可整体约去，因此原始变量应与 ε ≡ 1 逐位一致（只剩浮点舍入）：

| 算例 | ε₀ | max Δρ/ρ | max Δu/u | max Δp/p |
|---|---|---|---|---|
| Sod（含间断） | 0.6 | 9.5e-13 | 4.5e-12 | 1.8e-12 |
| 等熵涡（光滑） | 0.6 | 6.0e-15 | 1.8e-15 | 1.3e-15 |

Sod 上稍大（1e-12）是间断处的舍入被非线性放大；光滑算例上是机器精度。
两者都比"真正的建模错误"（会给出 O(1) 偏差）小 10 个量级以上。

### T3 静水平衡

u = 0、ρ = p = 1 全场均匀、ε(x) = 0.8 + 0.2 sin(2πx)（壁面边界，t = 1.0，200 节点）：

| 配置 | max\|u\| | max\|u\|/c | max\|Δρ\|/ρ |
|---|---|---|---|
| `porosity_gradient_force = 1`（ε∇p，默认） | **4.64e-16** | **3.92e-16** | 2.22e-16 |
| `porosity_gradient_force = 0`（∇(εp)，对照） | 1.82e-01 | 1.58e-01 | 4.43e-01 |

默认配置下速度停在浮点舍入水平（机器零），说明 `p∇ε` 源项与界面加权**精确抵消**（2.3 节）。
对照组（去掉 `p∇ε`）速度涨到 0.158 且密度剖面被改掉 44%，说明这个检验**确实有分辨力**：
只要符号、系数、或"用界面值还是节点值"写错一处，残差就会变成 `O(p|∇ε|)` 量级，
速度在几十步内涨到 O(1)。

### T4 全局守恒（壁面 + 非均匀 ε）

41×41×3 网格、六个壁面、ε = 0.9 + 0.1 sin(2πx)、爆炸初值、t = 0.2：

| 积分量 | 初始 | 末态 | 相对漂移 |
|---|---|---|---|
| ∭ ερ dV | 2.8366875000e-01 | 2.8366875000e-01 | **0.000e+00** |
| ∭ εE dV | 1.2650223214e+00 | 1.2650223214e+00 | **0.000e+00** |
| ∭ ερu dV | 0.0000000000e+00 | −5.7425575823e-03 | 绝对变化 5.74e-03（预期） |

质量与能量**逐位守恒**。动量不守恒是物理的：壁面压强与孔隙率梯度力 `p∇ε` 都在改变气相动量。
（初值关于中心对称，所以 ∭ερu 初值为 0，只能报绝对变化。）

---

## 6. 实施中查出并修复的两个真实缺陷

T4 一开始是 FAIL 的（质量 +4.26%、能量 +4.33%）。定位过程如下，两个缺陷**都不在 ε 逻辑里**，
而是本来就存在于单相求解器的边界处理中，被"全局守恒账本"这个新判据逼了出来。

### 6.1 现象与二分定位

| 变体 | Δm/m | ΔE/E |
|---|---|---|
| 单相、二阶 | 4.261e-02 | 4.335e-02 |
| 两相、ε ≡ 1、二阶 | **4.261e-02**（与单相完全相同） | 4.335e-02 |
| 单相、**一阶** | **0.000e+00** | 1.6e-16 |
| 两相、变 ε、一阶 | 0.000e+00 | 0.000e+00 |

两条结论立刻成立：① 单相算例同样泄漏、且数值与两相 ε≡1 **完全相同** ⇒ 与 ε 无关；
② 一阶格式精确守恒、二阶泄漏 ⇒ 问题出在**高阶重构**。

逐帧追踪漂移（`write_interval = 10`）显示：t < 0.074 时漂移恒为 1e-15 机器零，**冲击波一碰到
壁面就开始增长**，呈指数式爬升后饱和。这把范围锁定在壁面。

离线用 Python 复刻 `ReconstructFace + HLLC` 并按"虚拟节点 = 边界节点的镜像"建模，算出的
壁面质量通量**精确为 0**——与观测矛盾。于是给求解器加临时探针，直接在界面处打印原始模板值：

```
DEBUG RAW u: m=-4.740764E-02 g1(虚拟)=-4.740764E-02 bnd=4.740764E-02 pp=2.378483E-01
DEBUG face dir=0 at (2,23,3): W_L=(u -4.740764E-02) W_R=(u -1.589231E-02)
```

**关键一行是 `m` 与 `g1` 相等**：第二层虚拟节点（`m`）本该是 `-pp = -2.378e-01`，
实际却是 `-bnd = -4.7408e-02` —— 所有虚拟层都镜像了**同一个**边界节点。
于是左斜率 `(m - g1) = 0`、右斜率不为零，重构在壁面失去对称性：

```
W_L = g1 + ½·lim(0, bnd-g1)          = g1 = -bnd              （不动）
W_R = bnd - ½·lim(bnd-g1, pp-bnd)    = -0.0159  （被斜率修正拉动）
```

界面两侧速度同号（不镜像），质量通量变成 `O(ρu)`。用求解器累计整个边界面的质量通量，
量到 **0.276**（与实测 d(∭ερ)/dt ≈ 0.26 同量级），与解析模型的矛盾就此解决。

> 一阶格式为什么恰好为 0：一阶重构直接取 `W_L = g1`、`W_R = bnd`，无论虚拟节点取多少层都是
> 精确镜像对，所以完全掩盖了这个缺陷。**这也是"必须用二阶做守恒检验"的一个实例。**

### 6.2 缺陷一：反射界面的虚拟节点没有逐层镜像

**修法**：第 k 层虚拟节点 = 第 k 个内部节点关于壁面的镜像，即
`ghost(k) = mirror(node[is + k - 1])`。这样左斜率 `(m - g1)` 与右斜率 `(g1... )` 成同一对
参数，VanLeer 限制器对其对称 ⇒ 左右重构态精确成镜像 ⇒ 界面法向速度严格为 0。

### 6.3 缺陷二：虚拟节点用了不知道 ε 的变量换算

修完 6.2 后单相与 ε≡1 **精确守恒**，但变 ε 的两相算例仍漏 7.6%（一阶漏得更多）。原因是
两相把守恒量按"每单位总体积的气相量" `cons = ε·cons_intrinsic` 存储，而边界填虚拟节点时
仍调用 `Gas::Cons2Prim`（它假定 cons 是本征量）。于是 ε≠1 时虚拟节点的原始变量
`ρ_ghost = ε·ρ_src`、`p_ghost = ε·p_src` —— 被整体缩放了 ε，界面两侧同样不再对称。

**修法**：先 `/ε` 还原成本征守恒量再做镜像，写回时再 `×ε_ghost`
（新增小工具 `WeightCons`；单相权重恒 1.0，乘/除 1.0 精确，不影响已有的逐位退化）。

这两个缺陷都属于同一类错误：**虚拟节点的状态没有和离散格式"同构"地构造**。
它们不报错、不破坏稳定性、也不会在局部剖面里明显看出，只表现为守恒量缓慢漂移——
必须靠全局守恒账本才查得出来。

### 6.4 修复后的复验

| 变体 | Δm/m | ΔE/E |
|---|---|---|
| 单相、二阶 | **0.000e+00** | 1.6e-16 |
| 两相、ε ≡ 1、二阶 | **0.000e+00** | 1.6e-16 |
| 两相、变 ε、二阶 | **0.000e+00** | **0.000e+00** |
| 两相、变 ε、一阶 | **0.000e+00** | **0.000e+00** |
| 单相、一阶 | **0.000e+00** | 1.6e-16 |

**19/19 项 DEM 回归全部通过**，两相改动与边界修复对既有功能无副作用
（`sod` / `lax` / `vortex` 等用 `outlet` 边界的单相算例数值未变）。

---

## 7. 局限与到阶段 2 的接口

**当前局限**

1. ε 是**冻结的规定场**。颗粒接入后 ε 会随堆积结构变化，届时需处理 ∂ε/∂t 带来的
   几何功项（能量方程的 `p ∂ε/∂t`）。
2. 只支持 `inlet / outlet / wall` 三种边界；`farfield`、`symmetry`、周期性边界尚未实现
   （`EulerSolverStructUniform::BoundaryCondition` 遇到未识别名字会直接抛异常，不静默降级）。
3. 孔隙率梯度力 `p∇ε` 是一阶源项；ε 极陡（如床层自由面）时是整个格式里最粗糙的一环。
4. ε 变化剧烈且数值很小时源项刚度上升，建议对 CFL 做一次敏感性扫描（2.7 节）。
5. 单机 OpenMP，无 MPI。

**阶段 2 的接入点（已预留）**

- **ε 场已经是 `FieldData` 里的一个节点量**（`volume_fraction`），粒子→网格的体积分配
  只需往这个数组写，不用改求解器接口；
- 每次 `BoundaryCondition()` 会自动把 ε 的虚拟节点按零梯度补好；
- 相间力/相间功可以挂在 `CalcSourceResidual()`（目前与 `NSSolverStruct` 一样是空实现），
  与通量残差同一层级、同样参与 RK 子步；
- 时间步协调挂在 `EulerTwoPhaseStructUniform::CalcMinTimeStep` 之后
  （DEM 子循环 + 相间力在 CFD 步内做时间平均）。

**下一步建议**：先把**粒子→网格体积分配**做出来并单独验证（体积守恒精确成立、
网格加密后 ε 场收敛），再动相间力——否则 ε 分配误差和曳力模型误差会混在一起。

