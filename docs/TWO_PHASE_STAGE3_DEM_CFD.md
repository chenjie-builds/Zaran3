# DEM-CFD 单向耦合求解器（阶段 3）—— 架构、算例与验证

> **阶段目标**：把阶段 1 的「气相体积分数加权 Euler 求解器」与阶段 2 的
> 「粒子 → 网格映射器」接起来，做成一个**单向耦合**的 DEM-CFD 求解器：
> 流体对颗粒施加曳力与压力梯度力，**颗粒不对流体反馈**；
> 本阶段颗粒保持静止、不推进位移。
> 以**激波管**为基础算例，并给出可判定的验证阶梯。
>
> 上游：`docs/TWO_PHASE_EULER_STAGE1.md`（ε 加权 Euler）、
> `docs/TWO_PHASE_STAGE2_MAPPING.md`（粒子 → 网格映射）。

---

## 0. 「单向耦合」到底耦合了什么（先把边界划清）

| 链路 | 本阶段 | 说明 |
|---|---|---|
| **几何**：粒子占据体积 ⇒ 气相体积分数 ε_g < 1 | ✅ 做 | 粒子不动 ⇒ ε_g **只算一次**，是冻结场 |
| **动量**：流体 → 颗粒的曳力 | ✅ 做 | 每个时间步按当地气态 + 局部空隙率算 |
| **动量**：流体 → 颗粒的压力梯度力 −V_p∇p | ✅ 做 | 激波/间断处它与曳力同量级 |
| **动量**：颗粒 → 流体的反作用力 | ❌ **不做** | 本阶段的定义；双向耦合见第 9 节 |
| **位移**：颗粒受力后运动、ε_g 随时间变化 | ❌ **不做** | `coupling.freeze_particles = true` |
| 能量/质量交换（对流换热、相变、燃烧） | ❌ 不做 | 预留 |

**一个必须说清的点**：把 ε_g 放进流体方程（体积排斥）本身是一种几何耦合，
它**会**改变流场。所以「单向耦合不扰动流体」这句话的正确表述是：

> **ε_g ≡ 1 时**，DEM-CFD 的整条 CFD 链路必须与单相 Euler 求解器**逐位相同**。

这不是特判，是算术保证的（`ε_face = ½(1+1) = 1.0`，乘/除 1.0 在 IEEE 下精确；
`p∇ε` 的中心差分减出 `0.0`，加 0.0 也是精确的空操作）。V3 就是这句话的回归。

---

## 1. 架构：模块与数据接口

```
                    ┌─────────────────────────── DEM 侧（粒子几何） ───────────────────────────┐
 particles.csv ──► ReadDEMParticle ──► std::vector<DEMCFDCoupler::Particle>
 （与 DEM 算例同一格式）                          │  (x,y,z,r,ρ_p,u_p)
                                                 │
                              ┌──────────────────┴──────────────────┐
                              │           DEMCFDCoupler             │ ← 唯一的相间接口
                              │  · 持有 ParticleGridMapper          │
                              │  · ComputeGasVolumeFraction()  DEM→CFD
                              │  · EvaluateParticleForces()    CFD→DEM
                              └──────┬───────────────────┬──────────┘
                                     │                   │
                    DEM → CFD：α_s = Σ w_pc V_p/V_cell   │  CFD → DEM：Σ w_pc·(ρ,u,p,μ)
                    ε_g = 1 − α_s（CFD 物理节点）        │            Σ w_pc·∇p
                                     │                   │
                                     ▼                   ▼
        DataManagerNSTwoPhase::volume_fraction     DragForceModel::Evaluate()
                                     │              F = V_p·(β/ε_s)·u_rel − V_p∇p
                                     │                   │
                                     ▼                   ▼
        EulerTwoPhaseStructUniform                 particle_force.csv
        （ε 加权通量 + p∇ε 源项）                  particles_force_<iter>.dat
   ┌────────────────────────────────────┐
   │ NSFieldSimulation 的时间推进骨架   │  ← 三个虚钩子，DEMCFDSimulation 只实现它们
   │  PrepareCoupling()  → 注入 ε       │
   │  CouplingPreStep()  → （双向耦合在这里散射反作用力）
   │  CouplingPostStep() → 算力 + 输出  │
   └────────────────────────────────────┘
```

### 1.1 新增模块

| 文件 | 作用 | 可否独立验证 |
|---|---|---|
| `inc/Basic/TwoPhase/DragForceModel.h` + `src/.../DragForceModel.cpp` | 相间力模型（Stokes / Schiller–Naumann / Wen–Yu / Ergun / Gidaspow）+ 压力梯度力。**纯函数**，不含网格与求解器状态 | ✅ F1–F6 |
| `inc/Basic/TwoPhase/DEMCFDCoupler.h` + `src/.../DEMCFDCoupler.cpp` | 相间数据接口：粒子 → ε_g、网格 → 粒子插值、逐粒子受力；含无量纲 ↔ SI 换算 | ✅ V2/V4/V8 |
| `inc/Main/DEMCFDSimulation.h` + `src/.../DEMCFDSimulation.cpp` | 算例控制器：读参数/粒子、注入 ε、每个写盘步输出受力 | ✅ V1/V6/V7 |

### 1.2 改动（尽量小、尽量不动既有行为）

| 文件 | 改动 |
|---|---|
| `inc/Zone/SolverParam/FlowSolverParamTwoPhase.h/.cpp` | `VolumeFractionType` 增加 `External`（ε 由外部注入）；日志改成打印可读的类型名 |
| `src/Zone/Solver/EulerTwoPhaseStructUniform.cpp` | `InitField()`：`External` 分支跳过解析 ε，改为**只补 ghost**（零梯度外推）+ 范围检查 |
| `inc/Basic/TwoPhase/ParticleGridMapper.h` | 暴露 `ComputeWeights()` 与 `GetGridSpec()`——插值必须与体积分配共用**同一套**权重 |
| `inc/Main/FieldSimulation.h` + `src/Main/FieldSimulation.cpp` | `SolveField()` 改为虚函数；新增三个耦合钩子（默认空实现 ⇒ 既有算例行为不变）与 `GetFieldManager()` |
| `src/Main/Application.cpp` | 增加 `task.solver = "EulerDEM"` 分支（流体侧仍复用两相 Euler 场，只是换控制器） |
| `inc/Basic/PostProcess/Visual.h/.cpp` | 新增通用的 `WritePointsTecplotASCII()`（任意列的点云 → Tecplot ASCII） |
| `src/Basic/Physics/PerfectGas.cpp` | **修 bug**：`CalcMul` 的 Sutherland 调用参数错位（见第 8 节） |

### 1.3 交互流程（一个时间步，按实际调用顺序）

```
PrepareCoupling()                       ← 只在最开始调用一次（在 solver->Init() 之前！）
   ReadParameters / ReadParticles
   DEMCFDCoupler::Init(CFD 网格视图, 单位换算, 建模选项)
   ComputeGasVolumeFraction()  →  ε_g 写进 DataManagerNSTwoPhase
   （之后 solver->Init() → InitField() 只补 ghost，不再覆盖物理节点）

while (ContinueSolve()):
   CouplingPreStep(iter)                ← 单向耦合：空
   PreSolve()   → solver->Preprocess() = BoundaryCondition()
                                        （先填 ε 的 ghost，再按 ε 加权填流场 ghost）
   SolveOneStep() → solver->Solve() = TimeAdvance() → RungeKutta()
                                        （通量已按 ε_face 加权，动量含 +p∇ε）
   PostSolve()  → solver->Postprocess() = Cons2Prim → CheckPrimtive → FixPrimtive
                  + 按 write_interval 输出流场
   CouplingPostStep(iter)               ← 单向耦合的主场：
        UpdateViscosityField()          （用求解器的 Sutherland 算 μ*）
        ComputePressureGradient()       （中心差分，边界退化单侧）
        EvaluateParticleForces()        （同一套权重插值 → DragForceModel）
        → particle_force.csv / Tecplot 点云
```

---

## 2. 相间力模型（`DragForceModel`）

### 2.1 统一的写法（所有模型都归到同一个式子）

```
F_drag = V_p · (β/ε_s) · (u_g − u_p)
```

`β` 是 CFD-DEM 惯用的"相间动量交换系数"（单位**混合物体积**），`β/ε_s` 是
**单位颗粒体积**上的曳力系数。它与教科书的"单颗粒曳力 × 空隙率修正"是同一个式子：

```
β/ε_s = (3/4)·C_D·ρ_g·ε_g^(−n)·|u_rel| / d_p
V_p·(3/4)/d_p = (π d³/6)(3/4)/d_p = π d²/8 = ½·A_p
⇒  F_drag = ½·ρ_g·C_D·A_p·|u_rel|·u_rel · ε_g^(−n)
```

把五种模型（含稠密支的 Ergun）都归到 `β/ε_s` 上，力的方向、符号、零速度处理
就只有一份代码，不会出现"某个模型把相对速度写反"这种只差一个符号的 bug。

⚠️ **β 与 β/ε_s 差一个 ε_s，极易搞混**。本模块返回的是 `β/ε_s`；
校验脚本 F4 特意把两者都打印出来对照（第一版就写错在这里）。

### 2.2 五个模型

| 模型 | β/ε_s | C_D / Re 口径 |
|---|---|---|
| `stokes` | `18μ/d²` | 解析极限 `F = 3πμd·u_rel`；对 `|u|=0` 也良定义 |
| `schiller_naumann` | `0.75·C_D·ρ·|u|/d` | `Re = ρd|u|/μ`（无空隙率修正，无界绕流） |
| `wen_yu` | `0.75·C_D·ρ·|u|·ε_g^(−2.65)/d` | `Re = ε_g·ρd|u|/μ`（Gidaspow 原文口径） |
| `ergun` | `150μ·ε_s/(ε_g d²) + 1.75ρ|u|/d` | 不定义 C_D（记 0） |
| `gidaspow`（默认） | ε_g ≥ 阈值 → `wen_yu`；否则 → `ergun` | 阈值默认 0.8 |

`C_D(Re)`：`Re < 1000` 用 `24/Re·(1+0.15Re^0.687)`，`Re ≥ 1000` 取 `0.44`；
`Re ≤ 0` 时 `C_D ≡ 0`（与"相对速度为零 ⇒ 曳力为零"自洽，且避免 `24/0` 写出 `inf`
把整个 CSV/输出文件毁掉）。

**Gidaspow 的分段在阈值处是间断的**（实测 ε_g 从 0.80 降到 0.7999 时
`β/ε_s` 从 2.87e2 跳到 4.80e2，×1.67）。这是原文的形式，不是实现缺陷；
耦合器把每个粒子的 `branch` 列写进 CSV，用来识别算例落在哪一支。

### 2.3 压力梯度力

```
F_∇p = −V_p·∇p          ∇p 用 CFD 节点的中心差分（物理范围边界退化为单侧差分）
```

它与气相的 `porosity_gradient_force = 1`（ε∇p 形式）配套：把两相动量方程写成
`∇·(ερuu + εpI)` 再加源项 `+p∇ε`，等效于 `ε∇p`，而颗粒侧取 `−V_p∇p` 正是
"浮力"的标准写法。F4 验证了这条恒等式：

```
F_drag + F_∇p ≡ (−dp/dz)·V_p/ε_s      （压降反推的单颗粒总力）
```

---

## 3. 单位约定（最容易埋雷的地方）

求解器内部**全部无量纲**，颗粒参数是 **SI**。换算只在耦合器里做一次：

| 量 | 无量纲 → SI | 说明 |
|---|---|---|
| 长度 | `x = x*·L_ref` | `L_ref = [freestream] ref_length` |
| 密度 | `ρ = ρ*·ρ_ref` | `ρ_ref = [freestream] density` |
| 速度 | `u = u*·a_ref` | `a_ref = √(γ R T_ref)`（**声速**，不是 `L_ref/t_ref`） |
| 压强 | `p = p*·ρ_ref a_ref²` | 注意 `p_ref = ρ_ref a_ref²`，而物理上 `p = ρRT = ρa²/γ` |
| 温度 | `T = T*·T_ref`，`T* = γp*/ρ*` | 与 `Gas::CalcTemperature` 同一口径 |
| 粘性 | `μ = μ*·ρ_ref a_ref L_ref` | 与 `PerfectGas::m_mu0` 的定义严格互逆 |
| 力 | `F = F*·ρ_ref a_ref² L_ref²` | `[freestream]` 的 `ref_length=1`、`density=1`、`temperature=273` 给出 `F_ref = 1.097e5 N` |

参考量一律从**求解器自己的 `Dimensionless` 对象**取（`FlowSolverParam::GetDimensionless()`），
不另抄一份 `[freestream]` 的解析式 —— 避免"另一处手写的 a_ref"与求解器漂移。
`coupling_report.csv` 会把 6 个参考量写进去，校验脚本因此也不需要重复这些常数。

---

## 4. 算例：`tests/demcfd_shock_tube`

### 4.1 三步跑完

```bash
python tests/demcfd_shock_tube/generate_case.py     # 生成网格/粒子/控制文件 + 参考算例
bin/Release/Zaran3.10.2.exe tests/demcfd_shock_tube # 跑
python tests/demcfd_shock_tube/verify.py            # 独立校验（自己重跑并比对，共 89 项断言）
```

`verify.py` 默认会用当前二进制**重跑**主算例、`ref_sod` 与 ε≡1 变体，
保证校验的是同一份二进制；`--no-run` 只读现有输出，`--only F` 只跑模型组。

### 4.2 物理设置

| 项 | 值 |
|---|---|
| 计算域 | x ∈ [0,1]，y,z ∈ [0,0.04]（1 维问题放进 3 维薄网格，y/z 间距 ≥ x 间距） |
| 网格 | 401 × 3 × 3 节点 ⇒ `dx = 2.5e-3`，`dy = dz = 2.0e-2` |
| 初始条件 | 标准 Sod：左 (ρ,u,p) = (1, 0, 1)，右 (0.125, 0, 0.1)，膜片 x = 0.5 |
| 离散 | MUSCL 二阶 + HLLC + SSP-RK2，CFL = 0.5；物理时间跑到 `t = 0.2`（~401 步，~4 s） |
| 边界 | 六面 outlet |
| 无量纲参考量 | ρ_ref = 1 kg/m³、T_ref = 273 K、γ = 1.4、空气 ⇒ a_ref = 331.25 m/s，p_ref = 1.0973e5 Pa |

由此得到的物理图像：左室 1.08 bar / 382 K，右室 0.11 bar / 306 K；
激波速度 `U_shock = 1.7522`（无量纲）= 580 m/s，星区 `p* = 0.30313`、`u* = 0.92745`。

### 4.3 粒子怎么设（`particles.csv`，与 DEM 算例同一个格式）

```
id,group,radius,mass,px,py,pz,vx,vy,vz,ox,oy,oz,motion_type
1,0,4.0e-4,2.6808e-7,0.35,0.02,0.02,0,0,0,0,0,0,0
...
8,1,4.0e-4,2.6808e-7,0.45125,0.02,0.02,0,0,0,0,0,0,0
```

* **半径 / 密度与网格绑定**：`r = 0.4 mm`（`d_p = 0.8 mm`, `ρ_p = 1000 kg/m³`）
  ⇒ `Δx/d_p = 3.125`，`Δy/d_p = Δz/d_p = 25`。
  阶段 2 实测的体积分数法适用下限是 `Δ/d_p ≳ 3`，生成器对三个方向分别算并给判定；
  改半径必须同时改 `--nx`（或改 `--mapping-method exact`）。
* **密度**由 `mass/V_p` 给出（与 DEM 输入口径一致）；`mass = 0` 时退回
  `[coupling] particle_density`。单向耦合下密度不参与受力，但报告里要用。
* **位置**：
  * 5 个粒子放在**膜片下游** x = 0.55…0.75（每 0.05 一个），截面中心
    （正好落在 y/z 的中间节点上）。下游初速为 0 ⇒ **激波到达之前曳力恒等于 0**，
    力-时间历程上的台阶就是激波到达时刻，可直接与解析激波速度对照（V7）。
  * 2 个放在**上游** x = 0.35 / 0.45：它们只经历**稀疏扇**（光滑加速），
    用来检验曳力在光滑区的精度，不经过任何间断（V5 里这一侧的光滑区偏差 3.8e-8）。
  * 1 个放在 `x = 0.45 + 0.5Δx = 0.45125`（`group = 1`）：SubCell(2) 的 8 个子立方体
    中心会分到相邻两个单元（权重 0.5/0.5），**插值不再退化成"取最近节点"**。
    其余粒子都恰好落在节点上，只有它真正验证了插值算子（V4）。
* **速度**：`freeze_particles = true` 会把输入速度强制为 0 并告警 —— 本阶段颗粒静止。
* 粒子必须落在 CFD 的**控制体范围**里（比节点跨度多出半格：
  `[x0−Δx/2, x0+(n−½)Δx]`）；质心在域外的粒子会被计数并告警。

### 4.4 输出哪些量

| 文件 | 内容 |
|---|---|
| `result/<iter>.dat` | 流场（Tecplot ORDERED，9 个变量：`X,Y,Z,Density,Velocity_x/y/z,Pressure,Volume_fraction`）|
| `result/particle_force.csv` | **逐粒子、逐步**的 36 列：气态（ρ_g,u,v,w,p,T,μ）、`ε_g`、`Σw`、触及单元数、`∇p`、`Re`、`C_d`、`branch`、`β/ε_s`、曳力三分量、压力梯度力三分量、合力三分量与模 |
| `result/particles_force_<iter>.dat` | 粒子点云 + 力矢量（`Fx,Fy,Fz,F_mag,Re,Cd,ε_g,Rho_gas,P_gas,U_gas_*`），可在 Tecplot 里与流场云图叠加、做时间动画 |
| `result/coupling_report.csv` | 机器可读的关键量：参考量、`Δ/d_p` 三个方向、`ε_g` 范围、`α_max`、配分残差、末态 ΣF、最大 |F| |
| `result/residual.dat` | 流场残差历史（基类输出） |
| `ref_sod/` | **同网格、同离散参数**的单相 Euler 参考算例（逐位退化检验的对照） |

### 4.5 关键参数（`[coupling]` 段）

```toml
particle_file    = "particles.csv"
drag_model       = "gidaspow"     # stokes | schiller_naumann | wen_yu | ergun | gidaspow
mapping_method   = "subcell"      # Δ/d_p ≳ 3 时够用；1–3 用 exact
mapping_sub      = 2              # SubCell 每轴细分数（2 = 经典八分体法）
volume_fraction_from_particles = true   # false ⇒ ε_g ≡ 1（严格"无体积排斥"）
freeze_particles = true           # 本阶段颗粒静止
include_pressure_gradient_force = true
dilute_re_uses_voidage = true     # Wen-Yu 的 Re 含 ε_g
dense_threshold  = 0.8
force_write_interval   = 1        # 每 N 步写一行 CSV（0 = 跟随 write_interval）
write_coupling_frames  = true     # 每个写盘步额外写粒子点云
```

同时控制文件里必须有 `init.volume_fraction.type = "external"` ——
否则会被解析场覆盖，表现为"粒子完全不起作用"；`DEMCFDSimulation::BindField()`
会直接拦下来报错。

---

## 5. 验证阶梯与实测结果

`verify.py` 一共 **89 项断言**，全部通过（F 组 32 项 + V 组 57 项）。

### 5.1 F 组：相间力模型（不需要 CFD）

| 检验 | 判据 | 实测 |
|---|---|---|
| F1 Stokes 极限：`V_p·(β/ε_s)·u` vs 解析 `3πμd·u` | 相对偏差 ≤ 1e-14 | **≤ 3.05e-16** |
| F1 `C_D·Re/24 → 1`（Re→0） | ≤ 5e-3 | Re=1e-6 → 1.000011；Re=1e-3 → 1.001303 |
| F2 `C_D` 与标准单球阻力曲线对照（Re=0.1…1e4） | ≤ 10%（SN 的标称精度） | **3.08% / 4.15% / 1.25% / 2.03% / 4.35% / 7.32%** |
| F2 `C_D` 单调不增、Re≥1000 取 0.44 平台 | 精确 | ✅ |
| F3 `V_p(β/ε_s)u ≡ ½ρC_DA|u|u·ε_g^(−2.65)` | ≤ 1e-15 | **≤ 2.7e-16** |
| F3 Gidaspow（ε≥0.8）与 Wen–Yu 逐位一致 | 逐位 | ✅ |
| F4 Ergun 压降 + 气相动量平衡反推 β | ≤ 1e-13 | **0.0 / 2.1e-16 / 2.2e-16** |
| F4 `F_drag + F_∇p ≡ (−dp/dz)V_p/ε_s` | ≤ 1e-14 | **0.0 / 2.0e-16** |
| F5 零相对速度 ⇒ Re=C_D=0 且力恒为 0 | 精确 | ✅（5 个模型） |
| F5 `cos(F_drag, u_rel) = 1` | ≤ 1e-14 | **1.0000000000000000** |
| F5 压力梯度力 = −V_p∇p | 逐位 | **0.000e+00** |
| F6 Gidaspow 分段判定（ε=0.9…0.5） | 分支正确 | ✅，阈值跳变 ×1.67 |

### 5.2 V 组：耦合算例

| 检验 | 判据 | 实测 |
|---|---|---|
| **V1 结构** | 表头列数 == 数据列数；整型列原文为整数；行数 = 帧数×粒子数；时间严格递增 | 36 列 ×3216 行 / 402 帧 / 8 粒子 ✅ |
| V1 点云与 CSV 互证 | ≤ 1e-12 | **0.000e+00** |
| **V2 体积配分（端到端）** | \|∭(1−ε)dV − ΣV_p\|/ΣV_p ≤ 1e-9 | **1.032e-13**（C++ 内部 9.642e-16）|
| V2 `1−ε_min == α_max`（报告与流场互证） | ≤ 1e-12 | ✅ |
| **V3 ε≡1 与单相 Euler 逐位退化** | 逐位（== 0） | **0.000e+00**（82 帧 × 3609 点 × 8 变量）|
| **V4 插值**：从 `.dat` 用同一套 SubCell 权重重算 ρ_g/u_g/ε_g | ≤ 1e-11 | **≤ 2.22e-16**；`Σw = 1.000000000`；错位粒子 `cells = 2` |
| **V5 精确 Riemann 解对照**（气态） | 光滑区（±5Δx 内近似均匀）≤ 2% | 稀疏扇侧 **3.84e-8**；激波侧 **1.78e-2** |
| V5 全体记录的相对 L2（含间断） | 仅报告 | ρ 8.1e-3 / 4.6e-2；u 2.4e-2 / 8.4e-2；p 1.0e-2 / 4.7e-2 |
| **V6 力-时间历程 vs 精确解预测** | 光滑子集 rms/峰值 ≤ 10% | **≤ 0.19%**；峰值比 0.973–1.002 |
| **V7 激波到达时刻** vs `(x_p−x_d)/U_shock` | ≤ 5% | **4.81% / 2.40% / 1.58% / 1.21% / 0.96%** |
| V7 激波到达前曳力恒为 0 | 逐位（≤1e-12·F_max） | **712 条记录全为 0** |
| **V8 压力梯度力复核**（从 `.dat` 独立算 ∇p） | ≤ 1e-9 | **0.000e+00** |
| **V9 对称性** | ΣFy、ΣFz 机器零 | **1.46e-16 / 1.32e-16 N**（尺度 14.3 N）|
| V9 寄生 v/w 相对量 | ≤ 1e-3 | **2.46e-5**（薄网格固有，见下）|
| V10 体积排斥**确实**改变流场 | 差非零且 < 5% | **max\|Δρ\| = 2.10e-4** ✅ |
| **V11 粘性系数 = 标准 Sutherland 律** | ≤ 0.5% | **最大 1.067e-4** |

### 5.3 几条"为什么这样卡"的说明

* **V3 是整条阶梯里最强的一条**：它同时覆盖了 ε 的注入路径、ghost 填充、
  通量加权、`p∇ε` 源项、时间推进与输出。它一旦为 0，就说明"单向耦合对流体零扰动"
  这句话在本实现里没有例外。
* **V5/V6/V7 的严格判据只加在"精确解局部均匀"的记录上**（判据：±5Δx 内相对变化
  < 1e-3）。含激波/接触间断的记录数值解必然被抹开，那里的偏差是**离散误差**
  而不是模型误差；把它们混进来只会得到一个随网格变化的、没有信息量的数字。
  稀疏扇侧的光滑区偏差 3.8e-8 说明二阶格式在光滑区几乎精确。
* **V6 的归一化尺度**取该粒子**精确解的峰值力**。用"光滑子集的 RMS"归一化时，
  上游粒子的"零区"（曳力恒为 0）会让分母趋零、指标虚高到 1.0 —— 这是度量的错，
  不是解的错。
* **V9 的 v/w = 2.46e-5 不是耦合引入的**：1 维初值放进 3 维薄网格时，y/z 方向的
  通量差本应恒为 0，但 MUSCL 在这两个方向的斜率是舍入噪声，`1/Δy = 50` 的放大
  让它在几百步后长到 ~1e-5。**同样的量级在 `ref_sod` 里也存在**（V3 的逐位退化
  证明两边完全一致）。要压掉它得用真正的 1 维网格或对称边界。

---

## 6. 稠密变体：让 ε 真正起作用

默认算例是**稀疏**布置（全域平均 α_s = 1.3e-6），这是刻意的 ——
只有 ε_g ≈ 1，流场才近似自由场 Sod，V5–V7 才有解析对照可言。

要让体积排斥真的改变流场，需要一个**颗粒堆积带**：

```bash
python tests/demcfd_shock_tube/generate_case.py --out tests/demcfd_shock_tube_plug \
    --packing plug --nx 41 --radius 0.004 --plug-width 0.20
bin/Release/Zaran3.10.2.exe tests/demcfd_shock_tube_plug
```

实测（**100 个粒子**，r = 4 mm，`nx = 41` ⇒ Δx = 25 mm，堆积带宽 0.20 m = 8 个单元）：

| 量 | 值 |
|---|---|
| ε_g 范围（粒子处） | **[0.6247, 0.8392]**，全域 α_max = 0.3753 |
| 体积配分残差 | **2.27e-15** |
| 曳力分支 | **Ergun（稠密）8160 条 + Wen–Yu（稀疏）340 条** —— 两支都被真正跑到 |
| ΣF_x / max\|F\| | 204.5 N / 2.59 N |
| ΣF_y / ΣF_z | 1.78e-14 / 1.97e-14 N（机器零） |
| 与单相参考的最大密度差 | **max\|Δρ\| = 1.73e-2**（稀疏算例是 2.10e-4，**大 80 倍**）|

末帧（t = 0.2）密度剖面对照，耦合算例 vs 单相同网格参考：

| x | 0.20 | 0.40 | 0.475 | 0.55 | 0.60 | 0.70 | 0.85 | 0.95 |
|---|---|---|---|---|---|---|---|---|
| 耦合 | 0.99996 | 0.64698 | 0.49435 | 0.42285 | 0.42941 | 0.32866 | 0.23662 | 0.12500 |
| 单相 | 0.99996 | 0.64697 | 0.49490 | 0.42486 | 0.43063 | 0.32918 | 0.23655 | 0.12500 |

（`nx = 41` 很粗：稀疏波扇只跨约 5 个单元、接触间断被抹到 2–3 个单元，
所以 0.475/0.85 处是过渡值而非平台值。**这个变体的价值在于让 ε 真正起作用、
并跑到 Ergun 支**，不在定量精度。）

⚠️ 该变体**不适用** `verify.py` 的 V5–V7：堆积带会反射/透射激波，
"把精确自由场解在粒子位置取样"这个前提不再成立。它是**演示**算例，不是验证算例。
生成器也会提示 `Δy/d_p = 2.5 < 3`（此时 SubCell 的 α 场偏抹平，要更严格就减小半径
或加密网格）。**改参数后建议直接删掉该目录重跑** —— 帧名是迭代号，
两次不同参数的运行可能留下帧号错位的旧文件（本工程踩过同类陷阱）。

---

## 7. 与阶段 1/2 的验收是怎么拼起来的

| 阶段 | 验收对象 | 关键判据 | 状态 |
|---|---|---|---|
| 1 | ε 加权 Euler 求解器 | ε≡1 逐位退化；均匀 ε 不变性；静水平衡；壁面 + 非均匀 ε 守恒 | ✅ |
| 2 | 粒子 → 网格映射器 | 体积精确配分；α 有界；力散射守恒；对偶插值 | ✅ |
| 3 | 单向耦合 | **ε≡1 逐位退化（V3）**；插值同源（V4）；力-时间历程对精确解（V5–V7）；ε 的影响非零且有界（V10） | ✅ |
| 3.5 | 双向耦合（未做） | 相间反作用力守恒 + 动量账本 + 最小流化速度 | — |

---

## 8. 过程中发现并修掉的两个 bug

### 8.1 `PerfectGas::CalcMul` 的 Sutherland 调用参数错位（既有代码，已修）

```cpp
// 原来
return m_mu0 * Southerland(T, m_T0, m_Ts);
// 4 参数 Southerland 的签名是 Southerland(T, mu0, T0, Ts) —— 第 2 个是 mu0！
// ⇒ m_T0 = 273.16/T_ref 被当成 μ0、m_Ts = 110.4/T_ref 被当成 T0
// ⇒ 空气在 439 K 下给出 μ = 1.35e-4 Pa·s 而不是 2.43e-5（偏大 5.6 倍）
// 修正
return Southerland(T, m_mu0, m_T0, m_Ts);
```

该函数此前**没有任何调用点**（粘性残差在现有算例里都没开），所以一直没暴露；
阶段 3 的曳力需要 `Re = ρ|u|d/μ`，才把它翻出来。修好后
`μ(T_ref) = 1.7151e-5 Pa·s`，与标准 Sutherland 律的相对偏差 ≤ 1.07e-4（V11）。

### 8.2 粒子力 CSV 的 printf 说明符比表头多一个（本次新写的代码，已修）

格式串里有 37 个说明符而表头只有 36 个列名。前 24 列类型恰好都对得上，
于是 `branch`（int）被 `%.17g` 打印成 `4.94e-324` 这种**次正规数**、
`beta_over_es`（double）被 `%d` 打印成一个随机整数。

阴险之处在于：**错位后的行仍然解析得动**（两列都还是"数字"），
`int(float("4.94e-324")) = 0` 更会把 `branch` 悄悄变成合法的 0。
只能靠两条断言抓住：

```
表头列数 == 数据行字段数
整型列（iter/group/cells/branch）的**原始文本**必须匹配 ^[+-]?\d+$
```

（后者不能写成 `int(float(...))` —— 那正是会静默通过的那条路。）

---

## 9. 局限与向双向耦合的扩展点

**本阶段的已知局限**

1. **无相间反作用力**：流场完全按自行演化；`CalcSourceResidual()` 仍是空实现。
2. **颗粒静止**：`freeze_particles = true`；ε_g 是冻结场，不随时间变化。
3. **曳力用插值气态**：`(ρ,u,p)` 与 ∇p 都按粒子→网格的同一套权重插值。
   Δ/d_p ≳ 3 时粒子小于单元，这是"点源法"的标准做法；
   但分辨率不足时（Δ/d_p < 3）插值到的气态并不能代表球体周围的真实绕流。
4. **颗粒自身也被算进自己的 ε_g**：粒子的曳力用的是插值到粒子中心的 ε_g，
   其中包含该粒子自己贡献的 α_s（本算例 ε_g = 0.99973 而不是 1，
   让曳力偏大 0.07%）。严格做法是"挖掉自贡献"或用单元值，这里保留并在 CSV 中
   显式输出 `ε_g`，便于识别。
5. **无湍流、无对流换热、无相变/燃烧**；粘性项未启用（`is_viscous = false`）。
6. **网格必须均匀**：ε 的映射器只支持均匀笛卡尔网格。
7. **单机串行散点累加**：粒子多、网格细时需要 cell list 与并行归约
   （阶段 2 已记录同一件事）。

**双向耦合要做的三件事（结构已经留好）**

1. **反作用力**：在 `CouplingPreStep()` 里把 `βV_p(u_g−u_p)`（颗粒受力取负）
   经 `ParticleGridMapper::ScatterParticleVector()` 散列到网格，
   存成 `DataManagerNS` 的一个源项数组，在 `CalcSourceResidual()` 里参与 RK 子步。
   `Σ_cell F = Σ_particle F` 这条守恒性质阶段 2 已验证到机器精度。
2. **颗粒运动**：同一个钩子里按 DEM 推进颗粒（子循环：DEM 的 dt 比 CFD 小 1–2 个量级），
   并把 `freeze_particles` 关掉、每步重建 ε_g（此时 ε 不再是冻结场，
   `VolumeFractionType::External` 的注入路径每步重跑一次即可）。
3. **动量账本**：用 `∭(ερu + Σm_p u_p)` 的总动量输运与边界通量对照，
   作为双向耦合的第一条验收判据（阶段 1 的壁面守恒检验可直接复用）。

**下一步的自然验证阶梯**

| 级别 | 算例 | 判据 |
|---|---|---|
| 1 | 单球定常绕流（给定解析或文献 `C_D–Re`） | 与 SN/Ergun 相关式一致 |
| 2 | 固定床压降 | 与 Ergun 方程一致（判别 Re 与 ε 的口径） |
| 3 | 最小流化速度 | 与 `u_mf = …` 理论值一致 |
| 4 | 单个颗粒在激波后的加速（双向） | 与点源法文献结果一致 |
| 5 | 激波穿过颗粒床（本节的 plug 变体 + 双向） | 激波衰减/反射系数 |

---

## 10. 文件清单

```
新增
  inc/Basic/TwoPhase/DragForceModel.h        src/Basic/TwoPhase/DragForceModel.cpp
  inc/Basic/TwoPhase/DEMCFDCoupler.h         src/Basic/TwoPhase/DEMCFDCoupler.cpp
  inc/Main/DEMCFDSimulation.h                src/Main/DEMCFDSimulation.cpp
  tests/demcfd_shock_tube/{generate_case.py, verify.py, zaran.toml, mesh.*, particles.csv}
  tests/demcfd_shock_tube/ref_sod/（单相参考算例，由生成器产出）
  tests/demcfd_shock_tube_plug/（稠密变体演示，由生成器产出）
  docs/TWO_PHASE_STAGE3_DEM_CFD.md（本文）

改动
  inc/Zone/SolverParam/FlowSolverParamTwoPhase.h      + VolumeFractionType::External
  src/Zone/SolverParam/FlowSolverParamTwoPhase.cpp    + external 解析 / 可读日志
  src/Zone/Solver/EulerTwoPhaseStructUniform.cpp      + External 分支（只补 ghost）
  inc/Basic/TwoPhase/ParticleGridMapper.h             + ComputeWeights / GetGridSpec
  inc/Main/FieldSimulation.h  src/Main/FieldSimulation.cpp
                                                      + 三个耦合钩子 + 虚 SolveField
  src/Main/Application.cpp                            + task.solver = "EulerDEM"
  inc/Basic/PostProcess/Visual.h  src/.../Visual.cpp  + WritePointsTecplotASCII
  src/Basic/Physics/PerfectGas.cpp                    ！CalcMul 的 Sutherland 调用修正
```
