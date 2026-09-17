# 均匀结构网格 Euler 求解器（Zaran3 框架集成版）

> 求解器：`inc/Zone/Solver/EulerSolverStructUniform.h` + `src/Zone/Solver/EulerSolverStructUniform.cpp`
> 参数类：`inc/Zone/SolverParam/FlowSolverParamUniform.h` + `src/Zone/SolverParam/FlowSolverParamUniform.cpp`
> 场类：`inc/Zone/Field/NSFieldStructUniform.h` + `src/Zone/Field/NSFieldStructUniform.cpp`
> 算例与验证：`benchmarks/euler_uniform/{make_euler_case.py, verify_euler_case.py}`

---

## 1. 定位与遵循的框架约定

`EulerSolverStructUniform` 是 Zaran3 中按**项目既有架构规范**实现的均匀笛卡尔结构网格
Euler 方程求解器，作为后续 DEM-CFD 两相耦合的流体基座。它严格遵循以下框架约定：

| 约定 | 落实方式 |
|---|---|
| **继承自 `NSSolver`** | 继承链 `Solver → FieldSolver → FlowFieldSolver → NSSolver → EulerSolverStructUniform`；`Init / Solve / Preprocess / Postprocess / CalcResidual / TimeAdvance` 等主流程全部复用基类实现 |
| **用 `FieldData` + `DataManager` 管理流场数据** | 数据存于网格**节点**，变量名与其它求解器完全一致：`primitive_{0..4}`、`conservative_{0..4}`、`conservative_old_{0..4}`、`residual_{0..4}`、`limiter_{0..4}`、`dt`、`temperture`。半点（面）数据用 `DataManagerNSStruct` 的 `midnode_prim_left/right` 与 `midnode_flux` |
| **变量转换复用 `Gas`** | `Gas::Prim2Cons` / `Cons2Prim`，因此无量纲化、状态方程、`Dimensionless` 参考量与其它求解器完全一致 |
| **复用框架 Riemann 求解器** | `RiemannSolverBuilder` 创建 HLLC / Roe / VanLeer / AUSMPW / StegerWarming |
| **边界条件走 `GetBoundMap()`** | 支持 `inlet` / `outlet` / `wall`（网格文件中分别用编号 5 / 6 / 2） |
| **结果输出 Tecplot** | 场的类型是 `FieldType::NS_Structured`，`Visual::WriteTecplotBinary(shared_ptr<NSFieldStruct>)` 直接可用；输出的变量为 `X, Y, Z, Density, Velocity_x/y/z, Pressure` |

### 与 `NSSolverStruct` 的差别

`NSSolverStruct` 面向**贴体曲线坐标**，需要度量系数、雅可比矩阵、坐标变换系数等一整套机制。
本求解器限定网格为**均匀笛卡尔网格**，因此：

- `CalcCoordTransCoef()` 只做两件事：**校验网格严格均匀**（最大相对间距偏差超过 `1e-6` 直接报错）、
  记录 `dx / dy / dz`；
- 雅可比恒为 1，通量差分即普通中心差；
- 不需要 `Metric` 相关的任何存储（`[structure]` 节的度量参数对本求解器不生效，
  保留它们只是为了让控制文件结构统一）。

---

## 2. 类结构与挂载点

```
Solver
 └── FieldSolver            (FieldSolverType::Euler_Struct_Uniform   ← 新增枚举值)
      └── FlowFieldSolver
           └── NSSolver
                └── EulerSolverStructUniform          ← 本求解器
                     ├── 用 GridStruct（节点 + ghost，与 NS_Struct 共用）
                     ├── 用 DataManagerNSStruct（与 NS_Struct 共用）
                     ├── 用 FlowSolverParamUniform（新增参数类）
                     └── 由 NSFieldStructUniform 创建（新增场类）

Field → FieldNS → NSFieldStruct → NSFieldStructUniform   ← 只覆盖 AllocateSolver / AllocateSolverPara
```

**为什么 `FlowSolverParamUniform` 继承 `FlowSolverParamStruct` 而不是 `FlowSolverParam`**：
`FieldNS` / `NSFieldStruct` 内部会对参数做 `static_pointer_cast<FlowSolverParamStruct>`，
若参数类不继承它就会构成类型欺骗。继承 `FlowSolverParamStruct` 可保证这些既有代码安全，
同时仍能拿到 `[space] order` 等本求解器需要的参数。

**为什么 `NSFieldStructUniform` 不在 `NSFieldStruct` 之上另起炉灶**：
这样 Tecplot 输出、残差统计（`FieldNS::CalcResidual`）、多块交接面通信
（`FieldDataCommInfo`）等框架能力都可以直接复用，无需重复实现。

### 改动清单

| 文件 | 改动 |
|---|---|
| `inc/Zone/Solver/EulerSolverStructUniform.h`、`src/.../*.cpp` | **新增**求解器 |
| `inc/Zone/SolverParam/FlowSolverParamUniform.h`、`src/.../*.cpp` | **新增**参数类 |
| `inc/Zone/Field/NSFieldStructUniform.h`、`src/.../*.cpp` | **新增**场类 |
| `inc/Zone/Solver/FieldSolver.h` | 枚举增加 `FieldSolverType::Euler_Struct_Uniform` |
| `src/Generator/Field/FieldGenerator.cpp` | `Create()` 中增加该类型的场创建分支 |
| `src/Main/Application.cpp` | `task.solver = "Euler"` 的分派分支（要求 `grid_type = "Structured"`） |
| `inc/Zone/SolverParam/FlowSolverPara.h`、`src/.../*.cpp` | `InitFieldType` 增加 `Riemann1D`（一维激波管初值） |
| `src/Main/FieldSimulation.cpp` | 增加 `output.tecplot_ascii` 开关（默认仍输出二进制 `.plt`） |
| `inc/Basic/Math/EulerExactRiemann.h`、`src/.../*.cpp` | 精确 Riemann 解（验证用，从原 `inc/EulerUniform/` 移入标准位置） |

---

## 3. 离散方案

- **方程**：三维无粘可压缩 Euler 方程，理想气体；
- **网格**：均匀笛卡尔结构网格，节点中心有限差分/有限体积；
- **空间**：`dU_i/dt = -(F_{i+1/2} - F_{i-1/2})/Δx`，
  面通量由 Riemann 求解器给出，左右态可选
  一阶（直接取节点值）或二阶（MUSCL + 限制器，重构在**原始变量**上进行）；
- **时间**：SSP-RK，`cons = rk_coef[0]·cons_old + rk_coef[1]·cons + rk_coef[2]·dt·residual`
  （与 `NSSolverStruct::RungeKutta` 使用同一套 RK 系数约定，均匀网格下雅可比为 1）；
- **步长**：`dt = CFL / max_cell [ Σ_d (|v_d| + c)/Δx_d ]`。

**限制器映射**（`LimiterType` → 框架 `Limiter.h` 中的函数）：

| 控制文件 `space.limiter` | 行为 |
|---|---|
| `vk` | `LimiterVanLeer`（默认，平滑、保单调） |
| `barth` | `LimiterMinMod`（更保守） |
| `1st-order` | 斜率取 0，退化为 Godunov 一阶 |
| `noLimiter` | 中心差分斜率，不限制（仅用于对照，会产生振荡） |

---

## 4. 编译与运行

```bash
# 编译（Ninja + MSVC，与主程序共用构建目录）
cmake --build build-ninja

# 生成算例（mesh.dat / mesh.inp / zaran.toml）
python benchmarks/euler_uniform/make_euler_case.py --case sod \
    --out tests/euler_uniform_sod --tecplot-ascii

# 运行（参数为工作目录，控制文件固定为 <work_dir>/zaran.toml）
bin/Release/Zaran3.10.2.exe tests/euler_uniform_sod
```

控制文件中与本求解器相关的键：

```toml
[task]
solver     = "Euler"        # 选择均匀网格 Euler 求解器
grid_type  = "Structured"   # 复用结构网格
rk_stage   = 2              # 2 = SSP-RK2

[space]
order                = 2    # 1 = Godunov 一阶，2 = MUSCL 二阶
limiter              = "vk"
first_order_steps    = 0    # 前 N 步用一阶格式启动
require_uniform_grid = 1    # 校验网格严格均匀（推荐开启）
```

初值类型支持：`FarFlow` / `FarFieldNoVelocity` / `Explosion` / `Vortex`（等熵涡，框架自带）/
**`Riemann1D`（新增，一维激波管）**。

---

## 5. 算例

| 算例 | 初值 | 边界 | 用途 |
|---|---|---|---|
| `sod` | `Riemann1D`：左 (1,0,1)、右 (0.125,0,0.1)，膜片 x=0.5 | 全 `outlet` | 稀疏波/接触/激波三波结构 |
| `lax` | `Riemann1D`：左 (0.445,0.698,3.528)、右 (0.5,0,0.571) | 全 `outlet` | 强激波鲁棒性 |
| `strong` | `Riemann1D`：左 (1,0.75,1)、右 (0.125,0,0.1) | 全 `outlet` | 非零初速强激波 |
| `vortex` | `Vortex`（Hu-Shu 等熵涡，静止远场） | 全 `outlet` | **光滑解，测收敛阶** |
| `conserve` | `Explosion` | 全 `wall` | 封闭系统守恒性 |

> 说明：`Vortex` 初值取自框架既有的 `NSSolverStruct::InitFieldVortex`（远场静止的等熵涡）。
> 采用静止涡而非标准版（远场 u=1）的好处是涡心始终位于域中心，不会靠近周期边界，
> 因而参考解里**不需要叠加周期像**——这一点在实际调试中非常关键（详见第 7 节）。

---

## 6. 验证结果

环境：Windows / MSVC 14.51 / Release / CFL = 0.5 / HLLC / `limiter = vk`。

### 6.1 一维激波管（对 Toro 精确解）

| 阶数 | N | L1(ρ) | L1(u) | L1(p) | L∞(ρ) | 观测阶 |
|---|---|---|---|---|---|---|
| 1 | 100 | 2.0476e-02 | 3.9945e-02 | 1.8849e-02 | 8.3343e-02 | — |
| 1 | 200 | 1.3260e-02 | 2.2993e-02 | 1.1407e-02 | 8.8622e-02 | 0.63 |
| 1 | 400 | 8.4021e-03 | 1.2955e-02 | 6.7920e-03 | 8.6596e-02 | 0.66 |
| 2 | 100 | 5.2046e-03 | 9.6929e-03 | 3.9187e-03 | 7.0884e-02 | — |
| 2 | 200 | 2.8240e-03 | 4.8759e-03 | 1.9468e-03 | 6.7252e-02 | 0.88 |
| 2 | 400 | 1.5736e-03 | 2.5185e-03 | 9.8743e-04 | 6.9221e-02 | 0.84 |

符合预期：含间断算例的 L1 阶只有 0.6–0.9，**L∞ 不收敛**（有 O(1) 平台）；
但二阶的 L1 误差约为同网格一阶的 **1/3–1/4**。

### 6.2 二维等熵涡（光滑解）

| 阶数 | N×N | L1(ρ) | L1(u) = L1(v) | L∞(ρ) | 观测阶 |
|---|---|---|---|---|---|
| 1 | 41 | 5.6405e-03 | 9.9338e-03 | 1.0959e-01 | — |
| 1 | 81 | 3.2040e-03 | 5.5828e-03 | 5.4617e-02 | 0.82 |
| 1 | 161 | 1.7257e-03 | 2.9785e-03 | 2.5890e-02 | 0.89 |
| 2 | 41 | 5.0460e-04 | 9.2692e-04 | 7.4926e-03 | — |
| 2 | 81 | 1.3751e-04 | 2.5367e-04 | 2.2341e-03 | 1.88 |
| 2 | 161 | 3.1915e-05 | 5.6006e-05 | 8.6995e-04 | **2.11** |

**这是二阶精度的硬证据**：一阶格式的观测阶趋向 1.0，MUSCL 二阶格式趋向 2.1。
（`L1(u)` 与 `L1(v)` 完全相等是正确特征——静止涡的速度分量关于涡心对称。）

### 6.3 与改造前自包含实现的对照

改造前的独立实现（已移至 `_eulerbench/stale_standalone/`）在同一算例上的结果：

| 格式 | 自包含版 L1(ρ) @N=400 | 框架版 L1(ρ) @N=400 |
|---|---|---|
| 一阶 | 8.436e-03 | **8.402e-03** |
| 二阶 | 1.638e-03 | **1.574e-03** |

一阶几乎逐位吻合，二阶框架版略优 —— 说明移植过程没有引入数值退化。

### 6.4 网格校验与框架一致性

运行日志中的确认项：

```
Euler uniform solver parameters: reconstruction_order=2, require_uniform_grid=1, limiter=3, cfl=0.5
Riemann 1-D initial field: diaphragm x = 0.5, left = (rho 1, u 0, p 1), right = (rho 0.125, u 0, p 0.1)
Euler uniform grid: dim=3, nodes=(206,9,9), dx=5.025126E-03, dy=2.000000E-02, dz=2.000000E-02,
                    max spacing deviation=1.777831E-14
```

- 节点数 `(206,9,9)` = 物理 `(200,3,3)` + 两侧各 3 层 ghost，与框架一致；
- **最大相对间距偏差 1.78e-14**，确认网格严格均匀；
- Tecplot 输出文件 `result/<iter>.dat`（或 `.plt`），
  变量为 `X, Y, Z, Density, Velocity_x, Velocity_y, Velocity_z, Pressure`。

---

## 7. 实施中踩到的坑

### 7.1 MSVC 内部编译器错误 C1001（OpenMP）

`#pragma omp parallel for collapse(3)` 的循环边界若直接写 `static_cast<int>(nk)`，
MSVC 会报 `fatal error C1001: 内部编译器错误`。

**解决**：把范围提前取成普通 `int` 局部变量（见源文件中的 `LoopBounds` / `NodeCounts` 辅助结构），
循环条件里只出现简单变量。

### 7.2 `DataManagerNS` 的指针重载不完整

`SetPrim(int, const double*)`、`SetCons(int, const double*)`、`SetResidual(int, const double*)`
有指针重载，但 **`SetConsOld` 只有逐分量的三元重载**，写成 `SetConsOld(idx, cons)` 会编译失败。

### 7.3 二维等熵涡的周期像问题（重要经验）

在早前的自包含版本中，曾用「远场 u = 1」的标准 Hu-Shu 涡做验证。t = 2 时涡心移到 x = 7，
其**周期像**位于 x = −3，到 x = 0 处距离仅 3，`exp((1−9)/2) = 0.018`，
对 v 的贡献达 `0.796 × 3 × 0.018 = 4.4e-2` —— 与数值解远场观测到的 `max|v| ≈ 4.2e-2` 完全吻合。

**后果**：远场误差出现一个不随网格变化的平台，二阶观测阶从 2.1 掉到 0.32，看起来像格式失效。

**排查手法**：把误差**按到涡心的距离分环带**统计——若近场正常收敛、远场是不收敛的常数平台，
则问题在**参考解**而不在求解器。

**结论**：周期域基准算例的参考解必须计入周期像；或者直接改用**静止涡**
（框架自带的 `InitFieldType::Vortex` 正是静止涡，因此本版算例不受此问题影响）。

---

## 8. 与 DEM 耦合的衔接点

求解器已按框架规范就位，后续两相扩展只需在既有位置追加：

| 要加的东西 | 挂在哪里 |
|---|---|
| 体积分数 ε_g | 新增一个由 `DataManagerNS` 管理的场（如 `volume_fraction`），参与 `AddFluxResidual` 的面通量与残差 |
| 相间源项 S_mom | `CalcSourceResidual()`（目前为空实现，与 `NSSolverStruct` 一致），在 `m_residual` 上叠加 |
| 曳力 / ∇p 力 | 由 DEM 侧提供，按节点映射为源项；映射时的多线程归约可参考 `AddFluxResidual` 阶段 2 的写法 |
| 时间步协调 | `EulerSolverStructUniform::CalcMinTimeStep` 之后与 DEM 的 `dt` 做子循环匹配 |
| 精确解边界 | `inc/Basic/Math/EulerExactRiemann.h` 已备好（Toro 精确解，含强稀疏波下的括号保护） |
